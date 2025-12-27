/*
 * B3D is freely redistributable under the MIT License. See the file
 * "LICENSE" for information on usage and redistribution of this file.
 */

/* Mesh voxelizer implementation
 * Based on voxelizer by Karim Naaji (MIT License).
 *
 * Fixed-point design:
 * - Internal math uses b3d_scalar_t and B3D_FP_* macros
 * - Public API remains float for user convenience
 * - Conversion happens at API boundaries via B3D_FLOAT_TO_FP/B3D_FP_TO_FLOAT
 * - All multiplications use B3D_FP_MUL, divisions use B3D_FP_DIV
 */

#include <limits.h>
#include <stdbool.h>

#include "b3d-voxel.h"
#include "b3d.h"
#include "math-toolkit.h"

/* Precision factor to reduce "holes" artifact (relative to voxel_size) */
#define B3D_VOXEL_PRECISION_F 0.1f

/* Hash table size for deduplication (must be power of 2) */
#define B3D_VOXEL_HASH_SIZE 16384

/* Hash entry for open-addressing deduplication */
typedef struct {
    int ix, iy, iz; /* Integer grid coordinates */
    bool occupied;  /* Slot in use */
} b3d_voxel_hash_entry_t;

/* Hash table state for deduplication context */
typedef struct {
    b3d_voxel_hash_entry_t *entries;
    size_t capacity;
    size_t count;
} b3d_voxel_hash_t;

/* Fixed-point epsilon for comparisons */
#ifdef B3D_FLOAT_POINT
#define B3D_VOXEL_EPSILON 1e-7f
#define B3D_VOXEL_EPSILON_SQ4 (4.0f * B3D_VOXEL_EPSILON * B3D_VOXEL_EPSILON)
#define B3D_VOXEL_ZERO B3D_FLOAT_TO_FP(0.0f)
#else
/* In Q15.16: use ~0.005 (328 in Q15.16) for epsilon.
 * For squared comparison, must account for shift in b3d_voxel_mul:
 * epsilon_sq4 = (4 * 328 * 328) >> 16 = 6 in Q15.16
 */
#define B3D_VOXEL_EPSILON 328
#define B3D_VOXEL_EPSILON_SQ4 \
    ((4LL * B3D_VOXEL_EPSILON * B3D_VOXEL_EPSILON) >> B3D_FP_BITS)
#define B3D_VOXEL_ZERO 0
#endif

/* Internal 3D position using b3d_scalar_t for fixed-point compatibility */
typedef struct {
    b3d_scalar_t x, y, z;
} b3d_voxel_spos_t;

/* Internal triangle representation */
typedef struct {
    b3d_voxel_spos_t p[3];
} b3d_voxel_stri_t;

/*
 * Scalar operations using B3D_FP_* macros for fixed-point compatibility.
 * These work correctly in both float and fixed-point modes.
 */
static inline b3d_scalar_t b3d_voxel_mul(b3d_scalar_t a, b3d_scalar_t b)
{
    return B3D_FP_MUL(a, b);
}

static inline b3d_scalar_t b3d_voxel_div(b3d_scalar_t a, b3d_scalar_t b)
{
    return B3D_FP_DIV(a, b);
}

static inline b3d_scalar_t b3d_voxel_add(b3d_scalar_t a, b3d_scalar_t b)
{
    return B3D_FP_ADD(a, b);
}

static inline b3d_scalar_t b3d_voxel_sub_s(b3d_scalar_t a, b3d_scalar_t b)
{
    return B3D_FP_SUB(a, b);
}

static inline b3d_scalar_t b3d_voxel_neg(b3d_scalar_t a)
{
    return B3D_FP_SUB(B3D_VOXEL_ZERO, a);
}

static inline b3d_scalar_t b3d_voxel_min_s(b3d_scalar_t a, b3d_scalar_t b)
{
    return a < b ? a : b;
}

static inline b3d_scalar_t b3d_voxel_max_s(b3d_scalar_t a, b3d_scalar_t b)
{
    return a > b ? a : b;
}

static inline b3d_scalar_t b3d_voxel_abs_s(b3d_scalar_t x)
{
    return b3d_fp_abs(x);
}

/* Vector operations on internal scalar positions */
static inline b3d_scalar_t b3d_voxel_dot(b3d_voxel_spos_t a, b3d_voxel_spos_t b)
{
    return b3d_voxel_add(
        b3d_voxel_add(b3d_voxel_mul(a.x, b.x), b3d_voxel_mul(a.y, b.y)),
        b3d_voxel_mul(a.z, b.z));
}

static inline b3d_voxel_spos_t b3d_voxel_cross(b3d_voxel_spos_t a,
                                               b3d_voxel_spos_t b)
{
    return (b3d_voxel_spos_t) {
        b3d_voxel_sub_s(b3d_voxel_mul(a.y, b.z), b3d_voxel_mul(a.z, b.y)),
        b3d_voxel_sub_s(b3d_voxel_mul(a.z, b.x), b3d_voxel_mul(a.x, b.z)),
        b3d_voxel_sub_s(b3d_voxel_mul(a.x, b.y), b3d_voxel_mul(a.y, b.x)),
    };
}

static inline b3d_voxel_spos_t b3d_voxel_sub(b3d_voxel_spos_t a,
                                             b3d_voxel_spos_t b)
{
    return (b3d_voxel_spos_t) {
        b3d_voxel_sub_s(a.x, b.x),
        b3d_voxel_sub_s(a.y, b.y),
        b3d_voxel_sub_s(a.z, b.z),
    };
}

/* Convert internal scalar position to public float position */
static inline b3d_voxel_pos_t b3d_voxel_to_fpos(b3d_voxel_spos_t p)
{
    return (b3d_voxel_pos_t) {
        B3D_FP_TO_FLOAT(p.x),
        B3D_FP_TO_FLOAT(p.y),
        B3D_FP_TO_FLOAT(p.z),
    };
}

/*
 * Integer floor with overflow protection.
 * Clamps input to safe range before conversion to prevent undefined behavior.
 */
static inline int b3d_voxel_floori(b3d_scalar_t x)
{
#ifdef B3D_FLOAT_POINT
    /* Clamp to int range to prevent undefined behavior */
    if (x >= (float) INT_MAX)
        return INT_MAX;
    if (x <= (float) INT_MIN)
        return INT_MIN;
    int i = (int) x;
    return (x < 0.0f && x != (float) i) ? i - 1 : i;
#else
    /* Fixed-point: extract integer part, adjust for negative */
    int i = B3D_FP_TO_INT(x);
    /* If there's a fractional part and value is negative, subtract 1 */
    if (x < 0 && (x & ((1 << B3D_FP_BITS) - 1)) != 0)
        return i - 1;
    return i;
#endif
}

/*
 * Integer ceil with overflow protection.
 * Clamps input to safe range before conversion to prevent undefined behavior.
 */
static inline int b3d_voxel_ceili(b3d_scalar_t x)
{
#ifdef B3D_FLOAT_POINT
    /* Clamp to int range to prevent undefined behavior */
    if (x >= (float) INT_MAX)
        return INT_MAX;
    if (x <= (float) INT_MIN)
        return INT_MIN;
    int i = (int) x;
    return (x > 0.0f && x != (float) i) ? i + 1 : i;
#else
    /* Fixed-point: extract integer part, adjust for positive with fraction */
    int i = B3D_FP_TO_INT(x);
    /* If there's a fractional part and value is positive, add 1 */
    if (x > 0 && (x & ((1 << B3D_FP_BITS) - 1)) != 0)
        return i + 1;
    return i;
#endif
}

/* Plane-box overlap test */
static bool b3d_voxel_plane_box_overlap(b3d_voxel_spos_t normal,
                                        b3d_scalar_t d,
                                        b3d_voxel_spos_t halfbox)
{
    b3d_voxel_spos_t vmin, vmax;

    if (normal.x > B3D_VOXEL_ZERO) {
        vmin.x = b3d_voxel_neg(halfbox.x);
        vmax.x = halfbox.x;
    } else {
        vmin.x = halfbox.x;
        vmax.x = b3d_voxel_neg(halfbox.x);
    }

    if (normal.y > B3D_VOXEL_ZERO) {
        vmin.y = b3d_voxel_neg(halfbox.y);
        vmax.y = halfbox.y;
    } else {
        vmin.y = halfbox.y;
        vmax.y = b3d_voxel_neg(halfbox.y);
    }

    if (normal.z > B3D_VOXEL_ZERO) {
        vmin.z = b3d_voxel_neg(halfbox.z);
        vmax.z = halfbox.z;
    } else {
        vmin.z = halfbox.z;
        vmax.z = b3d_voxel_neg(halfbox.z);
    }

    if (b3d_voxel_add(b3d_voxel_dot(normal, vmin), d) > B3D_VOXEL_ZERO)
        return false;

    if (b3d_voxel_add(b3d_voxel_dot(normal, vmax), d) >= B3D_VOXEL_ZERO)
        return true;

    return false;
}

/* Triangle-box overlap test (Moller's SAT algorithm) */
static bool b3d_voxel_tri_box_overlap(b3d_voxel_spos_t boxcenter,
                                      b3d_voxel_spos_t halfbox,
                                      b3d_voxel_stri_t tri)
{
    b3d_voxel_spos_t v0, v1, v2;
    b3d_voxel_spos_t e0, e1, e2;
    b3d_scalar_t min, max, p0, p1, p2, rad, fex, fey, fez;

    /* Translate triangle to box center */
    v0 = b3d_voxel_sub(tri.p[0], boxcenter);
    v1 = b3d_voxel_sub(tri.p[1], boxcenter);
    v2 = b3d_voxel_sub(tri.p[2], boxcenter);

    /* Compute edge vectors */
    e0 = b3d_voxel_sub(v1, v0);
    e1 = b3d_voxel_sub(v2, v1);
    e2 = b3d_voxel_sub(v0, v2);

    /* AXISTEST_X01(e0.z, e0.y, fez, fey) */
    fex = b3d_voxel_abs_s(e0.x);
    fey = b3d_voxel_abs_s(e0.y);
    fez = b3d_voxel_abs_s(e0.z);

    p0 = b3d_voxel_sub_s(b3d_voxel_mul(e0.z, v0.y), b3d_voxel_mul(e0.y, v0.z));
    p2 = b3d_voxel_sub_s(b3d_voxel_mul(e0.z, v2.y), b3d_voxel_mul(e0.y, v2.z));
    if (p0 < p2) {
        min = p0;
        max = p2;
    } else {
        min = p2;
        max = p0;
    }
    rad = b3d_voxel_add(b3d_voxel_mul(fez, halfbox.y),
                        b3d_voxel_mul(fey, halfbox.z));
    if (min > rad || max < b3d_voxel_neg(rad))
        return false;

    /* AXISTEST_Y02(e0.z, e0.x, fez, fex) */
    p0 = b3d_voxel_add(b3d_voxel_neg(b3d_voxel_mul(e0.z, v0.x)),
                       b3d_voxel_mul(e0.x, v0.z));
    p2 = b3d_voxel_add(b3d_voxel_neg(b3d_voxel_mul(e0.z, v2.x)),
                       b3d_voxel_mul(e0.x, v2.z));
    if (p0 < p2) {
        min = p0;
        max = p2;
    } else {
        min = p2;
        max = p0;
    }
    rad = b3d_voxel_add(b3d_voxel_mul(fez, halfbox.x),
                        b3d_voxel_mul(fex, halfbox.z));
    if (min > rad || max < b3d_voxel_neg(rad))
        return false;

    /* AXISTEST_Z12(e0.y, e0.x, fey, fex) */
    p1 = b3d_voxel_sub_s(b3d_voxel_mul(e0.y, v1.x), b3d_voxel_mul(e0.x, v1.y));
    p2 = b3d_voxel_sub_s(b3d_voxel_mul(e0.y, v2.x), b3d_voxel_mul(e0.x, v2.y));
    if (p2 < p1) {
        min = p2;
        max = p1;
    } else {
        min = p1;
        max = p2;
    }
    rad = b3d_voxel_add(b3d_voxel_mul(fey, halfbox.x),
                        b3d_voxel_mul(fex, halfbox.y));
    if (min > rad || max < b3d_voxel_neg(rad))
        return false;

    /* Edge e1 tests */
    fex = b3d_voxel_abs_s(e1.x);
    fey = b3d_voxel_abs_s(e1.y);
    fez = b3d_voxel_abs_s(e1.z);

    /* AXISTEST_X01(e1.z, e1.y, fez, fey) */
    p0 = b3d_voxel_sub_s(b3d_voxel_mul(e1.z, v0.y), b3d_voxel_mul(e1.y, v0.z));
    p2 = b3d_voxel_sub_s(b3d_voxel_mul(e1.z, v2.y), b3d_voxel_mul(e1.y, v2.z));
    if (p0 < p2) {
        min = p0;
        max = p2;
    } else {
        min = p2;
        max = p0;
    }
    rad = b3d_voxel_add(b3d_voxel_mul(fez, halfbox.y),
                        b3d_voxel_mul(fey, halfbox.z));
    if (min > rad || max < b3d_voxel_neg(rad))
        return false;

    /* AXISTEST_Y02(e1.z, e1.x, fez, fex) */
    p0 = b3d_voxel_add(b3d_voxel_neg(b3d_voxel_mul(e1.z, v0.x)),
                       b3d_voxel_mul(e1.x, v0.z));
    p2 = b3d_voxel_add(b3d_voxel_neg(b3d_voxel_mul(e1.z, v2.x)),
                       b3d_voxel_mul(e1.x, v2.z));
    if (p0 < p2) {
        min = p0;
        max = p2;
    } else {
        min = p2;
        max = p0;
    }
    rad = b3d_voxel_add(b3d_voxel_mul(fez, halfbox.x),
                        b3d_voxel_mul(fex, halfbox.z));
    if (min > rad || max < b3d_voxel_neg(rad))
        return false;

    /* AXISTEST_Z0(e1.y, e1.x, fey, fex) */
    p0 = b3d_voxel_sub_s(b3d_voxel_mul(e1.y, v0.x), b3d_voxel_mul(e1.x, v0.y));
    p1 = b3d_voxel_sub_s(b3d_voxel_mul(e1.y, v1.x), b3d_voxel_mul(e1.x, v1.y));
    if (p0 < p1) {
        min = p0;
        max = p1;
    } else {
        min = p1;
        max = p0;
    }
    rad = b3d_voxel_add(b3d_voxel_mul(fey, halfbox.x),
                        b3d_voxel_mul(fex, halfbox.y));
    if (min > rad || max < b3d_voxel_neg(rad))
        return false;

    /* Edge e2 tests */
    fex = b3d_voxel_abs_s(e2.x);
    fey = b3d_voxel_abs_s(e2.y);
    fez = b3d_voxel_abs_s(e2.z);

    /* AXISTEST_X2(e2.z, e2.y, fez, fey) */
    p0 = b3d_voxel_sub_s(b3d_voxel_mul(e2.z, v0.y), b3d_voxel_mul(e2.y, v0.z));
    p1 = b3d_voxel_sub_s(b3d_voxel_mul(e2.z, v1.y), b3d_voxel_mul(e2.y, v1.z));
    if (p0 < p1) {
        min = p0;
        max = p1;
    } else {
        min = p1;
        max = p0;
    }
    rad = b3d_voxel_add(b3d_voxel_mul(fez, halfbox.y),
                        b3d_voxel_mul(fey, halfbox.z));
    if (min > rad || max < b3d_voxel_neg(rad))
        return false;

    /* AXISTEST_Y1(e2.z, e2.x, fez, fex) */
    p0 = b3d_voxel_add(b3d_voxel_neg(b3d_voxel_mul(e2.z, v0.x)),
                       b3d_voxel_mul(e2.x, v0.z));
    p1 = b3d_voxel_add(b3d_voxel_neg(b3d_voxel_mul(e2.z, v1.x)),
                       b3d_voxel_mul(e2.x, v1.z));
    if (p0 < p1) {
        min = p0;
        max = p1;
    } else {
        min = p1;
        max = p0;
    }
    rad = b3d_voxel_add(b3d_voxel_mul(fez, halfbox.x),
                        b3d_voxel_mul(fex, halfbox.z));
    if (min > rad || max < b3d_voxel_neg(rad))
        return false;

    /* AXISTEST_Z12(e2.y, e2.x, fey, fex) */
    p1 = b3d_voxel_sub_s(b3d_voxel_mul(e2.y, v1.x), b3d_voxel_mul(e2.x, v1.y));
    p2 = b3d_voxel_sub_s(b3d_voxel_mul(e2.y, v2.x), b3d_voxel_mul(e2.x, v2.y));
    if (p2 < p1) {
        min = p2;
        max = p1;
    } else {
        min = p1;
        max = p2;
    }
    rad = b3d_voxel_add(b3d_voxel_mul(fey, halfbox.x),
                        b3d_voxel_mul(fex, halfbox.y));
    if (min > rad || max < b3d_voxel_neg(rad))
        return false;

    /* Test AABB overlap in X, Y, Z */
    min = b3d_voxel_min_s(v0.x, b3d_voxel_min_s(v1.x, v2.x));
    max = b3d_voxel_max_s(v0.x, b3d_voxel_max_s(v1.x, v2.x));
    if (min > halfbox.x || max < b3d_voxel_neg(halfbox.x))
        return false;

    min = b3d_voxel_min_s(v0.y, b3d_voxel_min_s(v1.y, v2.y));
    max = b3d_voxel_max_s(v0.y, b3d_voxel_max_s(v1.y, v2.y));
    if (min > halfbox.y || max < b3d_voxel_neg(halfbox.y))
        return false;

    min = b3d_voxel_min_s(v0.z, b3d_voxel_min_s(v1.z, v2.z));
    max = b3d_voxel_max_s(v0.z, b3d_voxel_max_s(v1.z, v2.z));
    if (min > halfbox.z || max < b3d_voxel_neg(halfbox.z))
        return false;

    /* Test triangle plane vs box */
    b3d_voxel_spos_t normal = b3d_voxel_cross(e0, e1);
    b3d_scalar_t d = b3d_voxel_neg(b3d_voxel_dot(normal, v0));
    return b3d_voxel_plane_box_overlap(normal, d, halfbox);
}

/* Squared triangle area x4 (for degenerate triangle detection, avoids sqrt).
 * Returns |cross(ab,ac)|^2 = (2*area)^2, compare against B3D_VOXEL_EPSILON_SQ4.
 */
static inline b3d_scalar_t b3d_voxel_tri_area_sq4(b3d_voxel_stri_t *tri)
{
    b3d_voxel_spos_t ab = b3d_voxel_sub(tri->p[1], tri->p[0]);
    b3d_voxel_spos_t ac = b3d_voxel_sub(tri->p[2], tri->p[0]);
    b3d_voxel_spos_t cross = b3d_voxel_cross(ab, ac);

    /* |cross|^2 = cross.x^2 + cross.y^2 + cross.z^2 */
    return b3d_voxel_add(b3d_voxel_add(b3d_voxel_mul(cross.x, cross.x),
                                       b3d_voxel_mul(cross.y, cross.y)),
                         b3d_voxel_mul(cross.z, cross.z));
}

/* Integer-based spatial hash for O(1) voxel deduplication.
 * Uses prime multipliers for good distribution across grid coordinates.
 */
static inline size_t b3d_voxel_hash_int(int ix, int iy, int iz)
{
    /* Convert to unsigned to handle negative indices correctly */
    uint32_t ux = (uint32_t) ix;
    uint32_t uy = (uint32_t) iy;
    uint32_t uz = (uint32_t) iz;

    size_t h = ux * 73856093u;
    h ^= uy * 19349663u;
    h ^= uz * 83492791u;
    return h;
}

/* Initialize hash table for deduplication.
 * Returns true on success, false if memory not available.
 */
static inline bool b3d_voxel_hash_init(b3d_voxel_hash_t *ht,
                                       b3d_voxel_hash_entry_t *storage,
                                       size_t capacity)
{
    if (!ht || !storage || capacity == 0)
        return false;

    ht->entries = storage;
    ht->capacity = capacity;
    ht->count = 0;

    /* Clear all entries */
    for (size_t i = 0; i < capacity; i++)
        storage[i].occupied = false;

    return true;
}

/* Insert grid coordinates into hash table.
 * Returns true if inserted (new entry), false if already exists or table full.
 */
static inline bool b3d_voxel_hash_insert(b3d_voxel_hash_t *ht,
                                         int ix,
                                         int iy,
                                         int iz)
{
    if (!ht || ht->count >= ht->capacity)
        return false;

    size_t hash = b3d_voxel_hash_int(ix, iy, iz);
    size_t idx = hash % ht->capacity;

    /* Quadratic probing for collision resolution */
    for (size_t probe = 0; probe < ht->capacity; probe++) {
        size_t i = (idx + probe * probe) % ht->capacity;
        b3d_voxel_hash_entry_t *e = &ht->entries[i];

        if (!e->occupied) {
            /* Empty slot: insert new entry */
            e->ix = ix;
            e->iy = iy;
            e->iz = iz;
            e->occupied = true;
            ht->count++;
            return true;
        }

        /* Check if this coordinate already exists */
        if (e->ix == ix && e->iy == iy && e->iz == iz)
            return false; /* Duplicate */
    }

    return false; /* Table full, couldn't find slot */
}

/* Calculate mesh bounding box (public API uses float) */
void b3d_voxel_mesh_aabb(const float *triangles,
                         size_t tri_count,
                         b3d_aabb_t *aabb)
{
    if (!triangles || tri_count == 0 || !aabb)
        return;

    aabb->min.x = aabb->min.y = aabb->min.z = 1e30f;
    aabb->max.x = aabb->max.y = aabb->max.z = -1e30f;

    for (size_t i = 0; i < tri_count * 9; i += 3) {
        float x = triangles[i + 0];
        float y = triangles[i + 1];
        float z = triangles[i + 2];

        if (x < aabb->min.x)
            aabb->min.x = x;
        if (y < aabb->min.y)
            aabb->min.y = y;
        if (z < aabb->min.z)
            aabb->min.z = z;
        if (x > aabb->max.x)
            aabb->max.x = x;
        if (y > aabb->max.y)
            aabb->max.y = y;
        if (z > aabb->max.z)
            aabb->max.z = z;
    }
}

/* Main voxelization function.
 *
 * Uses hash-based O(1) deduplication when out_voxels is provided.
 * Early-exits when output buffer is full (returns max_voxels + 1 to indicate
 * truncation).
 *
 * When out_voxels is NULL (estimation mode), returns an upper bound that may
 * exceed actual unique voxel count due to overlap counting without dedup.
 */
size_t b3d_voxelize(const float *triangles,
                    size_t tri_count,
                    float voxel_size,
                    uint32_t color,
                    b3d_voxel_t *out_voxels,
                    size_t max_voxels)
{
    if (!triangles || tri_count == 0 || voxel_size <= 0.0f)
        return 0;

    /* Convert parameters to internal scalar type */
    b3d_scalar_t vs = B3D_FLOAT_TO_FP(voxel_size);
    b3d_scalar_t half_size = b3d_voxel_div(vs, B3D_FLOAT_TO_FP(2.0f));
    b3d_scalar_t precision =
        b3d_voxel_mul(vs, B3D_FLOAT_TO_FP(B3D_VOXEL_PRECISION_F));
    size_t voxel_count = 0;
    bool buffer_full = false;

    /* Static hash table for O(1) deduplication (~1MB footprint) */
    static b3d_voxel_hash_entry_t hash_storage[B3D_VOXEL_HASH_SIZE];
    b3d_voxel_hash_t ht = {0};

    if (out_voxels && max_voxels > 0)
        b3d_voxel_hash_init(&ht, hash_storage, B3D_VOXEL_HASH_SIZE);

    /* Process each triangle */
    for (size_t ti = 0; ti < tri_count && !buffer_full; ti++) {
        b3d_voxel_stri_t tri;
        const float *tp = triangles + ti * 9;

        /* Convert triangle vertices to internal scalar type */
        tri.p[0] = (b3d_voxel_spos_t) {
            B3D_FLOAT_TO_FP(tp[0]),
            B3D_FLOAT_TO_FP(tp[1]),
            B3D_FLOAT_TO_FP(tp[2]),
        };
        tri.p[1] = (b3d_voxel_spos_t) {
            B3D_FLOAT_TO_FP(tp[3]),
            B3D_FLOAT_TO_FP(tp[4]),
            B3D_FLOAT_TO_FP(tp[5]),
        };
        tri.p[2] = (b3d_voxel_spos_t) {
            B3D_FLOAT_TO_FP(tp[6]),
            B3D_FLOAT_TO_FP(tp[7]),
            B3D_FLOAT_TO_FP(tp[8]),
        };

        /* Skip degenerate triangles */
        if (b3d_voxel_tri_area_sq4(&tri) < B3D_VOXEL_EPSILON_SQ4)
            continue;

        /* Compute triangle AABB in scalar space */
        b3d_scalar_t min_x = b3d_voxel_min_s(
            tri.p[0].x, b3d_voxel_min_s(tri.p[1].x, tri.p[2].x));
        b3d_scalar_t min_y = b3d_voxel_min_s(
            tri.p[0].y, b3d_voxel_min_s(tri.p[1].y, tri.p[2].y));
        b3d_scalar_t min_z = b3d_voxel_min_s(
            tri.p[0].z, b3d_voxel_min_s(tri.p[1].z, tri.p[2].z));
        b3d_scalar_t max_x = b3d_voxel_max_s(
            tri.p[0].x, b3d_voxel_max_s(tri.p[1].x, tri.p[2].x));
        b3d_scalar_t max_y = b3d_voxel_max_s(
            tri.p[0].y, b3d_voxel_max_s(tri.p[1].y, tri.p[2].y));
        b3d_scalar_t max_z = b3d_voxel_max_s(
            tri.p[0].z, b3d_voxel_max_s(tri.p[1].z, tri.p[2].z));

        /* Convert AABB to integer grid indices */
        int ix_min = b3d_voxel_floori(b3d_voxel_div(min_x, vs));
        int ix_max = b3d_voxel_ceili(b3d_voxel_div(max_x, vs));
        int iy_min = b3d_voxel_floori(b3d_voxel_div(min_y, vs));
        int iy_max = b3d_voxel_ceili(b3d_voxel_div(max_y, vs));
        int iz_min = b3d_voxel_floori(b3d_voxel_div(min_z, vs));
        int iz_max = b3d_voxel_ceili(b3d_voxel_div(max_z, vs));

        b3d_voxel_spos_t halfbox = {
            b3d_voxel_add(half_size, precision),
            b3d_voxel_add(half_size, precision),
            b3d_voxel_add(half_size, precision),
        };

        /* Test each voxel in AABB using integer indices */
        for (int ix = ix_min; ix <= ix_max && !buffer_full; ix++) {
            b3d_scalar_t x = b3d_voxel_mul(B3D_INT_TO_FP(ix), vs);
            for (int iy = iy_min; iy <= iy_max && !buffer_full; iy++) {
                b3d_scalar_t y = b3d_voxel_mul(B3D_INT_TO_FP(iy), vs);
                for (int iz = iz_min; iz <= iz_max && !buffer_full; iz++) {
                    b3d_scalar_t z = b3d_voxel_mul(B3D_INT_TO_FP(iz), vs);
                    b3d_voxel_spos_t center = {x, y, z};

                    if (b3d_voxel_tri_box_overlap(center, halfbox, tri)) {
                        if (!out_voxels) {
                            /* Estimation mode: count all overlaps */
                            voxel_count++;
                        } else {
                            /* Use hash for O(1) deduplication */
                            if (b3d_voxel_hash_insert(&ht, ix, iy, iz)) {
                                /* New unique voxel */
                                if (voxel_count < max_voxels) {
                                    /* Convert back to float for output */
                                    out_voxels[voxel_count].pos =
                                        b3d_voxel_to_fpos(center);
                                    out_voxels[voxel_count].color = color;
                                    voxel_count++;
                                } else {
                                    /* Buffer full: early exit */
                                    buffer_full = true;
                                    voxel_count = max_voxels + 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return voxel_count;
}

/* Render voxels as cubes (lit version) */
void b3d_voxel_render(const b3d_voxel_t *voxels, size_t count, float voxel_size)
{
    if (!voxels || count == 0 || voxel_size <= 0.0f)
        return;

    float h = voxel_size * 0.5f;

    /* Cube vertices relative to center */
    static const float cube_verts[8][3] = {
        {-1, -1, -1}, {+1, -1, -1}, {+1, +1, -1}, {-1, +1, -1},
        {-1, -1, +1}, {+1, -1, +1}, {+1, +1, +1}, {-1, +1, +1},
    };

    /* Cube face definitions (2 triangles per face, 6 faces) */
    static const int cube_faces[12][3] = {
        /* Front face */
        {4, 5, 6},
        {4, 6, 7},
        /* Back face */
        {1, 0, 3},
        {1, 3, 2},
        /* Top face */
        {7, 6, 2},
        {7, 2, 3},
        /* Bottom face */
        {0, 1, 5},
        {0, 5, 4},
        /* Right face */
        {5, 1, 2},
        {5, 2, 6},
        /* Left face */
        {0, 4, 7},
        {0, 7, 3},
    };

    /* Face normals for lighting */
    static const float face_normals[6][3] = {
        {0, 0, +1}, /* Front */
        {0, 0, -1}, /* Back */
        {0, +1, 0}, /* Top */
        {0, -1, 0}, /* Bottom */
        {+1, 0, 0}, /* Right */
        {-1, 0, 0}, /* Left */
    };

    for (size_t i = 0; i < count; i++) {
        float cx = voxels[i].pos.x;
        float cy = voxels[i].pos.y;
        float cz = voxels[i].pos.z;
        uint32_t color = voxels[i].color;

        /* Render 12 triangles per cube (2 per face) */
        for (int f = 0; f < 12; f++) {
            b3d_tri_t tri;
            int face_idx = f / 2;

            for (int v = 0; v < 3; v++) {
                int vi = cube_faces[f][v];
                tri.v[v].x = cx + cube_verts[vi][0] * h;
                tri.v[v].y = cy + cube_verts[vi][1] * h;
                tri.v[v].z = cz + cube_verts[vi][2] * h;
            }

            b3d_triangle_lit(&tri, face_normals[face_idx][0],
                             face_normals[face_idx][1],
                             face_normals[face_idx][2], color);
        }
    }
}

/* Render voxels as cubes (unlit version) */
void b3d_voxel_render_flat(const b3d_voxel_t *voxels,
                           size_t count,
                           float voxel_size)
{
    if (!voxels || count == 0 || voxel_size <= 0.0f)
        return;

    float h = voxel_size * 0.5f;

    static const float cube_verts[8][3] = {
        {-1, -1, -1}, {+1, -1, -1}, {+1, +1, -1}, {-1, +1, -1},
        {-1, -1, +1}, {+1, -1, +1}, {+1, +1, +1}, {-1, +1, +1},
    };

    static const int cube_faces[12][3] = {
        {4, 5, 6}, {4, 6, 7}, {1, 0, 3}, {1, 3, 2}, {7, 6, 2}, {7, 2, 3},
        {0, 1, 5}, {0, 5, 4}, {5, 1, 2}, {5, 2, 6}, {0, 4, 7}, {0, 7, 3},
    };

    for (size_t i = 0; i < count; i++) {
        float cx = voxels[i].pos.x;
        float cy = voxels[i].pos.y;
        float cz = voxels[i].pos.z;
        uint32_t color = voxels[i].color;

        for (int f = 0; f < 12; f++) {
            b3d_tri_t tri;

            for (int v = 0; v < 3; v++) {
                int vi = cube_faces[f][v];
                tri.v[v].x = cx + cube_verts[vi][0] * h;
                tri.v[v].y = cy + cube_verts[vi][1] * h;
                tri.v[v].z = cz + cube_verts[vi][2] * h;
            }

            b3d_triangle(&tri, color);
        }
    }
}
