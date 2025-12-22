/*
 * B3D is freely redistributable under the MIT License. See the file
 * "LICENSE" for information on usage and redistribution of this file.
 */

/* Internal math types and functions. */

#ifndef B3D_MATH_H
#define B3D_MATH_H

#include <limits.h>
#include <math.h>
#include <stdint.h>

/*
 * Fixed-point / scalar arithmetic
 *
 * Supports two modes:
 *   - B3D_FLOAT_POINT: Use native float (faster on FPU-equipped hardware)
 *   - Default: Q15.16 fixed-point (better for embedded without FPU)
 */

#ifndef B3D_DEPTH_16BIT
#ifndef B3D_DEPTH_32BIT
#define B3D_DEPTH_32BIT 1
#endif
#endif

#if defined(B3D_FLOAT_POINT) && defined(B3D_DEPTH_16BIT)
#error "B3D_FLOAT_POINT and B3D_DEPTH_16BIT cannot be combined"
#endif

/* Floating-point path for comparison builds */
#ifdef B3D_FLOAT_POINT

typedef float b3d_scalar_t;
#define B3D_FP_BITS 0
#define B3D_FP_ONE 1.0f
#define B3D_FP_HALF 0.5f

#define B3D_INT_TO_FP(i) ((b3d_scalar_t) (i))
#define B3D_FLOAT_TO_FP(f) ((b3d_scalar_t) (f))
#define B3D_FP_TO_INT(f) ((int) (f))
#define B3D_FP_TO_FLOAT(f) ((float) (f))

#define B3D_FP_MUL(a, b) ((a) * (b))
#define B3D_FP_DIV(a, b) ((b) == 0 ? 0 : (a) / (b))
#define B3D_FP_FLOOR(f) floorf(f)
#define B3D_FP_ADD(a, b) ((a) + (b))
#define B3D_FP_SUB(a, b) ((a) - (b))

/* Fixed-point style constants expressed in float for compatibility */
#define B3D_FP_PI 3.1415926536f
#define B3D_FP_PI_HALF 1.5707963268f
#define B3D_FP_2PI 6.2831853072f

static inline b3d_scalar_t b3d_fp_sin(b3d_scalar_t x)
{
    return sinf(x);
}

static inline b3d_scalar_t b3d_fp_cos(b3d_scalar_t x)
{
    return cosf(x);
}

static inline b3d_scalar_t b3d_fp_sqrt(b3d_scalar_t a)
{
    return a <= 0 ? 0 : sqrtf(a);
}

static inline b3d_scalar_t b3d_fp_abs(b3d_scalar_t x)
{
    return fabsf(x);
}

#else /* Fixed-point path */

/* Q15.16 format in int32_t: range ±32,768, precision 1/65536 */
typedef int32_t b3d_fixed_t;
#define B3D_FP_BITS 16
#define B3D_FP_ONE (1LL << B3D_FP_BITS) /* 1LL avoids UB */
#define B3D_FP_HALF (1LL << (B3D_FP_BITS - 1))

/* Use multiplication (not left-shift) to avoid UB with negative values */
#define B3D_INT_TO_FP(i) ((b3d_fixed_t) ((int64_t) (i) * B3D_FP_ONE))
#define B3D_FLOAT_TO_FP(f) ((b3d_fixed_t) ((f) * B3D_FP_ONE))
#define B3D_FP_TO_INT(f) ((f) >> B3D_FP_BITS)
#define B3D_FP_TO_FLOAT(f) ((float) (f) / B3D_FP_ONE)

/* int64_t intermediates required for precision.
 * Use multiplication instead of left-shift to avoid UB with negative values. */
#define B3D_FP_MUL(a, b) ((b3d_fixed_t) (((int64_t) (a) * (b)) >> B3D_FP_BITS))
#define B3D_FP_DIV(a, b) \
    ((b) == 0 ? 0 : (b3d_fixed_t) (((int64_t) (a) * B3D_FP_ONE) / (b)))
#define B3D_FP_FLOOR(f)                                  \
    ((f) & ~((1 << B3D_FP_BITS) - 1)) /* replaces floorf \
                                       */
#define B3D_FP_ADD(a, b) ((a) + (b))
#define B3D_FP_SUB(a, b) ((a) - (b))

/* Define scalar type for fixed-point */
typedef b3d_fixed_t b3d_scalar_t;

/* Fixed-point constants */
#define B3D_FP_PI ((b3d_fixed_t) ((355LL << B3D_FP_BITS) / 113))
#define B3D_FP_PI_HALF (B3D_FP_PI >> 1)
#define B3D_FP_2PI (B3D_FP_PI << 1)
#define B3D_FP_PI_SQ \
    ((b3d_fixed_t) (((int64_t) B3D_FP_PI * B3D_FP_PI) >> B3D_FP_BITS))

/* Bhaskara I sine: sin(x) ≈ 16x(π-x) / (5π² - 4x(π-x)), ~0.3% max error */
static inline b3d_fixed_t b3d_fp_sin(b3d_fixed_t x)
{
    int sign = 1;

    /* Handle negative angles - guard against INT32_MIN overflow */
    if (x < 0) {
        if (x == INT32_MIN)
            x = INT32_MAX; /* Clamp to avoid overflow on negation */
        else
            x = -x;
        sign = -sign;
    }

    /* Fast modulo reduction using int64_t (avoids slow loop for large angles)
     */
    if (x >= B3D_FP_2PI) {
        int64_t x64 = (int64_t) x % (int64_t) B3D_FP_2PI;
        x = (b3d_fixed_t) x64;
    }

    /* Map (π, 2π) to (0, π) with sign flip */
    if (x > B3D_FP_PI) {
        x -= B3D_FP_PI;
        sign = -sign;
    }

    /* Bhaskara I approximation with precomputed pi² */
    b3d_fixed_t xp = B3D_FP_MUL(x, B3D_FP_PI - x);
    b3d_fixed_t denom = 5 * B3D_FP_PI_SQ - 4 * xp;

    if (denom == 0)
        return 0;

    return sign * B3D_FP_DIV(16 * xp, denom);
}

static inline b3d_fixed_t b3d_fp_cos(b3d_fixed_t x)
{
    /* Do addition in int64_t to avoid overflow for large angles near INT32_MAX.
     * b3d_fp_sin handles modulo reduction internally. */
    int64_t x64 = (int64_t) x + (int64_t) B3D_FP_PI_HALF;
    /* Reduce to int32_t range before calling sin */
    if (x64 > INT32_MAX)
        x64 = x64 % (int64_t) B3D_FP_2PI;
    else if (x64 < INT32_MIN)
        x64 = -((-x64) % (int64_t) B3D_FP_2PI);
    return b3d_fp_sin((b3d_fixed_t) x64);
}

/* Integer sqrt on Q16.16: computes floor(sqrt(a)) in fixed-point */
static inline b3d_fixed_t b3d_fp_sqrt(b3d_fixed_t a)
{
    if (a <= 0)
        return 0;

    /* Scale by 2^16 so sqrt preserves fixed-point fraction bits */
    uint64_t n = ((uint64_t) (uint32_t) a) << B3D_FP_BITS;
    uint64_t res = 0;
    uint64_t bit = 1ULL << 62; /* Largest power-of-four starting point */

    while (bit > n)
        bit >>= 2;

    while (bit != 0) {
        if (n >= res + bit) {
            n -= res + bit;
            res = (res >> 1) + bit;
        } else {
            res >>= 1;
        }
        bit >>= 2;
    }

    return res > INT32_MAX ? INT32_MAX : (b3d_fixed_t) res;
}

/* Fixed-point absolute value - guards against INT32_MIN overflow */
static inline b3d_fixed_t b3d_fp_abs(b3d_fixed_t x)
{
    if (x == INT32_MIN)
        return INT32_MAX; /* -INT32_MIN overflows; clamp to max */
    return x < 0 ? -x : x;
}

#endif /* B3D_FLOAT_POINT */

/*
 * Constants
 */
#define B3D_NEAR_DISTANCE 0.1f
#define B3D_FAR_DISTANCE 100.0f
#define B3D_EPSILON 1e-8f /* Near-zero threshold for division guards */
#define B3D_DEGEN_THRESHOLD 0.0001f /* Degenerate trig/scanline threshold */
#define B3D_CULL_THRESHOLD 0.01f /* Back-face culling dot product threshold */
#define B3D_DEPTH_FAR 1e30f      /* Depth buffer clear value (far plane) */
#define B3D_PI 3.1415926536f     /* Pi constant for angle conversions */
#define B3D_CLIP_BUFFER_SIZE 32  /* Maximum triangles in clipping buffer */

/*
 * Depth buffer type (shared with public header)
 */
#ifndef B3D_DEPTH_T_DEFINED
#ifdef B3D_FLOAT_POINT
typedef float b3d_depth_t;
#elif defined(B3D_DEPTH_16BIT)
typedef uint16_t b3d_depth_t;
#else
typedef int32_t b3d_depth_t;
#endif
#define B3D_DEPTH_T_DEFINED
#endif

#if defined(B3D_FLOAT_POINT)
#define B3D_DEPTH_CLEAR B3D_DEPTH_FAR
static inline b3d_depth_t b3d_depth_from_float(float d)
{
    return d;
}
static inline float b3d_depth_to_float(b3d_depth_t d)
{
    return d;
}
#elif defined(B3D_DEPTH_16BIT)
#define B3D_DEPTH_CLEAR 0xFFFF
static inline b3d_depth_t b3d_depth_from_float(float d)
{
    if (d < 0.0f)
        d = 0.0f;
    if (d > 1.0f)
        d = 1.0f;
    return (b3d_depth_t) (d * 65535.0f + 0.5f);
}
static inline float b3d_depth_to_float(b3d_depth_t d)
{
    return (float) d * (1.0f / 65535.0f);
}
#else
#define B3D_DEPTH_CLEAR INT32_MAX
static inline b3d_depth_t b3d_depth_from_float(float d)
{
    return B3D_FLOAT_TO_FP(d);
}
static inline float b3d_depth_to_float(b3d_depth_t d)
{
    return B3D_FP_TO_FLOAT(d);
}
#endif

/*
 * Internal vector/matrix types
 */
typedef struct {
    float x, y, z, w;
} b3d_vec_t;

typedef struct {
    float m[4][4];
} b3d_mat_t;

typedef struct {
    b3d_vec_t p[3];
} b3d_triangle_t;

/*
 * Vector operations
 */
static inline float b3d_vec_dot(b3d_vec_t a, b3d_vec_t b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline float b3d_vec_length(b3d_vec_t v)
{
    return sqrtf(b3d_vec_dot(v, v));
}

static inline b3d_vec_t b3d_vec_add(b3d_vec_t a, b3d_vec_t b)
{
    return (b3d_vec_t) {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}

static inline b3d_vec_t b3d_vec_sub(b3d_vec_t a, b3d_vec_t b)
{
    return (b3d_vec_t) {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}

static inline b3d_vec_t b3d_vec_mul(b3d_vec_t a, float b)
{
    return (b3d_vec_t) {a.x * b, a.y * b, a.z * b, a.w * b};
}

static inline b3d_vec_t b3d_vec_div(b3d_vec_t a, float b)
{
    return (b3d_vec_t) {a.x / b, a.y / b, a.z / b, 1};
}

static inline b3d_vec_t b3d_vec_cross(b3d_vec_t a, b3d_vec_t b)
{
    return (b3d_vec_t) {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                        a.x * b.y - a.y * b.x, 1};
}

static inline b3d_vec_t b3d_vec_norm(b3d_vec_t v)
{
    float l = b3d_vec_length(v);
    if (l < B3D_EPSILON)
        return (b3d_vec_t) {0, 0, 0, 1};
    return (b3d_vec_t) {v.x / l, v.y / l, v.z / l, 1};
}

/*
 * Matrix operations
 */

/* Build matrix from 4 row vectors */
static inline b3d_mat_t b3d_mat_from_rows(b3d_vec_t r0,
                                          b3d_vec_t r1,
                                          b3d_vec_t r2,
                                          b3d_vec_t r3)
{
    return (b3d_mat_t) {{
        {r0.x, r0.y, r0.z, r0.w},
        {r1.x, r1.y, r1.z, r1.w},
        {r2.x, r2.y, r2.z, r2.w},
        {r3.x, r3.y, r3.z, r3.w},
    }};
}

/* Extract row as vec3 (w=0) */
static inline b3d_vec_t b3d_mat_row3(b3d_mat_t m, int r)
{
    return (b3d_vec_t) {m.m[r][0], m.m[r][1], m.m[r][2], 0};
}

static inline b3d_mat_t b3d_mat_ident(void)
{
    return (b3d_mat_t) {{
        [0][0] = 1,
        [1][1] = 1,
        [2][2] = 1,
        [3][3] = 1,
    }};
}

static inline b3d_mat_t b3d_mat_rot_x(float a)
{
    return (b3d_mat_t) {{
        [0][0] = 1,
        [1][1] = cosf(a),
        [1][2] = sinf(a),
        [2][1] = -sinf(a),
        [2][2] = cosf(a),
        [3][3] = 1,
    }};
}

static inline b3d_mat_t b3d_mat_rot_y(float a)
{
    return (b3d_mat_t) {{
        [0][0] = cosf(a),
        [0][2] = sinf(a),
        [2][0] = -sinf(a),
        [1][1] = 1,
        [2][2] = cosf(a),
        [3][3] = 1,
    }};
}

static inline b3d_mat_t b3d_mat_rot_z(float a)
{
    return (b3d_mat_t) {{
        [0][0] = cosf(a),
        [0][1] = sinf(a),
        [1][0] = -sinf(a),
        [1][1] = cosf(a),
        [2][2] = 1,
        [3][3] = 1,
    }};
}

static inline b3d_mat_t b3d_mat_trans(float x, float y, float z)
{
    return (b3d_mat_t) {{
        [0][0] = 1,
        [1][1] = 1,
        [2][2] = 1,
        [3][3] = 1,
        [3][0] = x,
        [3][1] = y,
        [3][2] = z,
    }};
}

static inline b3d_mat_t b3d_mat_scale(float x, float y, float z)
{
    return (b3d_mat_t) {{
        [0][0] = x,
        [1][1] = y,
        [2][2] = z,
        [3][3] = 1,
    }};
}

static inline b3d_mat_t b3d_mat_proj(float fov,
                                     float aspect,
                                     float near,
                                     float far)
{
    fov = 1.0f / tanf(fov * 0.5f / 180.0f * B3D_PI);
    return (b3d_mat_t) {{
        [0][0] = aspect * fov,
        [1][1] = fov,
        [2][2] = far / (far - near),
        [3][2] = (-far * near) / (far - near),
        [2][3] = 1,
        [3][3] = 0,
    }};
}

static inline int b3d_mat_is_identity(b3d_mat_t m)
{
    return m.m[0][0] == 1.0f && m.m[1][1] == 1.0f && m.m[2][2] == 1.0f &&
           m.m[3][3] == 1.0f && m.m[0][1] == 0.0f && m.m[0][2] == 0.0f &&
           m.m[0][3] == 0.0f && m.m[1][0] == 0.0f && m.m[1][2] == 0.0f &&
           m.m[1][3] == 0.0f && m.m[2][0] == 0.0f && m.m[2][1] == 0.0f &&
           m.m[2][3] == 0.0f && m.m[3][0] == 0.0f && m.m[3][1] == 0.0f &&
           m.m[3][2] == 0.0f;
}

static inline b3d_mat_t b3d_mat_mul(b3d_mat_t a, b3d_mat_t b)
{
    if (b3d_mat_is_identity(b))
        return a;
    if (b3d_mat_is_identity(a))
        return b;

    b3d_mat_t matrix = {{{0}}};
    for (int c = 0; c < 4; c++) {
        for (int r = 0; r < 4; r++) {
            matrix.m[r][c] = a.m[r][0] * b.m[0][c] + a.m[r][1] * b.m[1][c] +
                             a.m[r][2] * b.m[2][c] + a.m[r][3] * b.m[3][c];
        }
    }
    return matrix;
}

static inline b3d_vec_t b3d_mat_mul_vec(b3d_mat_t m, b3d_vec_t v)
{
    return (b3d_vec_t) {
        v.x * m.m[0][0] + v.y * m.m[1][0] + v.z * m.m[2][0] + v.w * m.m[3][0],
        v.x * m.m[0][1] + v.y * m.m[1][1] + v.z * m.m[2][1] + v.w * m.m[3][1],
        v.x * m.m[0][2] + v.y * m.m[1][2] + v.z * m.m[2][2] + v.w * m.m[3][2],
        v.x * m.m[0][3] + v.y * m.m[1][3] + v.z * m.m[2][3] + v.w * m.m[3][3],
    };
}

/* Quick inverse for affine transforms with orthonormal rotation.
 * Only valid when the 3x3 rotation part is orthonormal (no scaling/shearing).
 * For general matrices, use full matrix inversion instead.
 */
static inline b3d_mat_t b3d_mat_qinv(b3d_mat_t m)
{
    /* Extract original rows for transpose and translation computation */
    b3d_vec_t m0 = b3d_mat_row3(m, 0);
    b3d_vec_t m1 = b3d_mat_row3(m, 1);
    b3d_vec_t m2 = b3d_mat_row3(m, 2);
    b3d_vec_t t = b3d_mat_row3(m, 3);

    /* Transpose 3x3 rotation, compute inverse translation: -R^T * t */
    return b3d_mat_from_rows(
        (b3d_vec_t) {m0.x, m1.x, m2.x, 0}, (b3d_vec_t) {m0.y, m1.y, m2.y, 0},
        (b3d_vec_t) {m0.z, m1.z, m2.z, 0},
        (b3d_vec_t) {-b3d_vec_dot(t, m0), -b3d_vec_dot(t, m1),
                     -b3d_vec_dot(t, m2), 1});
}

static inline b3d_mat_t b3d_mat_point_at(b3d_vec_t pos,
                                         b3d_vec_t target,
                                         b3d_vec_t up)
{
    b3d_vec_t fwd = b3d_vec_sub(target, pos);

    /* Guard: pos == target produces zero forward vector */
    if (b3d_vec_dot(fwd, fwd) < B3D_EPSILON)
        return b3d_mat_ident();

    fwd = b3d_vec_norm(fwd);

    /* Gram-Schmidt: remove forward component from up */
    up = b3d_vec_sub(up, b3d_vec_mul(fwd, b3d_vec_dot(up, fwd)));

    /* Guard: up parallel to forward produces zero up vector */
    if (b3d_vec_dot(up, up) < B3D_EPSILON)
        return b3d_mat_ident();

    up = b3d_vec_norm(up);
    b3d_vec_t right = b3d_vec_cross(up, fwd);

    return b3d_mat_from_rows((b3d_vec_t) {right.x, right.y, right.z, 0},
                             (b3d_vec_t) {up.x, up.y, up.z, 0},
                             (b3d_vec_t) {fwd.x, fwd.y, fwd.z, 0},
                             (b3d_vec_t) {pos.x, pos.y, pos.z, 1});
}

/*
 * Geometry operations
 */
static inline b3d_vec_t b3d_intersect_plane(b3d_vec_t norm,
                                            float plane_d,
                                            b3d_vec_t start,
                                            b3d_vec_t end)
{
    float ad = b3d_vec_dot(start, norm);
    float bd = b3d_vec_dot(end, norm);
    float denom = bd - ad;
    if (fabsf(denom) < B3D_EPSILON)
        return start;
    float t = (plane_d - ad) / denom;
    if (t < 0.0f)
        t = 0.0f;
    if (t > 1.0f)
        t = 1.0f;
    b3d_vec_t start_to_end = b3d_vec_sub(end, start);
    b3d_vec_t segment = b3d_vec_mul(start_to_end, t);
    return b3d_vec_add(start, segment);
}

static inline int b3d_clip_against_plane(b3d_vec_t plane,
                                         b3d_vec_t norm,
                                         b3d_triangle_t in,
                                         b3d_triangle_t out[2])
{
    norm = b3d_vec_norm(norm);
    float plane_d = b3d_vec_dot(norm, plane);
    b3d_vec_t *inside[3];
    int inside_count = 0;
    b3d_vec_t *outside[3];
    int outside_count = 0;
    float d0 = b3d_vec_dot(in.p[0], norm) - plane_d;
    float d1 = b3d_vec_dot(in.p[1], norm) - plane_d;
    float d2 = b3d_vec_dot(in.p[2], norm) - plane_d;
    if (d0 >= 0)
        inside[inside_count++] = &in.p[0];
    else
        outside[outside_count++] = &in.p[0];
    if (d1 >= 0)
        inside[inside_count++] = &in.p[1];
    else
        outside[outside_count++] = &in.p[1];
    if (d2 >= 0)
        inside[inside_count++] = &in.p[2];
    else
        outside[outside_count++] = &in.p[2];
    if (inside_count == 3) {
        out[0] = in;
        return 1;
    } else if (inside_count == 1 && outside_count == 2) {
        out[0].p[0] = *inside[0];
        out[0].p[1] =
            b3d_intersect_plane(norm, plane_d, *inside[0], *outside[0]);
        out[0].p[2] =
            b3d_intersect_plane(norm, plane_d, *inside[0], *outside[1]);
        return 1;
    } else if (inside_count == 2 && outside_count == 1) {
        out[0].p[0] = *inside[0];
        out[0].p[1] = *inside[1];
        out[0].p[2] =
            b3d_intersect_plane(norm, plane_d, *inside[0], *outside[0]);
        out[1].p[0] = *inside[1];
        out[1].p[1] = out[0].p[2];
        out[1].p[2] =
            b3d_intersect_plane(norm, plane_d, *inside[1], *outside[0]);
        return 2;
    }
    return 0;
}

#endif /* B3D_MATH_H */
