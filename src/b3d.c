/*
 * B3D is freely redistributable under the MIT License. See the file
 * "LICENSE" for information on usage and redistribution of this file.
 */

#include <stdbool.h>
#include <string.h>

#include "b3d.h"
#include "math-toolkit.h"

/* Global state */
int b3d_width, b3d_height;
uint32_t *b3d_pixels;
b3d_depth_t *b3d_depth;

static b3d_mat_t b3d_model, b3d_view, b3d_proj;
static b3d_vec_t b3d_camera;
static b3d_camera_t b3d_camera_params; /* Full camera state for queries */
static float b3d_fov_degrees;          /* Current FOV for queries */
static bool b3d_ortho_mode = false;    /* Orthographic vs perspective */

/* Matrix stack for push/pop operations */
static b3d_mat_t b3d_matrix_stack[B3D_MATRIX_STACK_SIZE];
static int b3d_matrix_stack_top = 0;

/* Lazy matrix computation state */
#ifdef B3D_NO_CULLING
static b3d_mat_t b3d_model_view; /* Cached model*view matrix */
#endif

/* Dirty flag for matrix recomputation */
static bool b3d_model_view_dirty = true;

/* Debug counters */
static size_t b3d_clip_drop_count = 0;

/* Lighting state */
static b3d_vec_t b3d_light_dir = {0.0f, 0.0f, 1.0f, 0.0f}; /* Default: +Z */
static float b3d_ambient = 0.2f;                           /* Default: 20% */

/* Perspective-correct depth conversion constants.
 * Converts interpolated 1/w to normalized depth [0, 1].
 * Formula: depth = B3D_DEPTH_OFFSET - w_inv * B3D_DEPTH_SCALE
 *
 * Fixed-point considerations:
 * - w_inv ranges from 1/far to 1/near (e.g., 0.01 to 10 for default planes)
 * - In Q15.16: 10 = 655360, 0.01 = 655 - both fit comfortably
 * - B3D_FP_MUL uses int64_t intermediate to prevent overflow
 * - Result depth is in [0, B3D_FP_ONE] range
 *
 * Precision note: In Q15.16, near-plane depth may be ~0.0001 instead of
 * exact 0 due to rounding in constant conversion. This is expected and
 * does not affect z-buffer correctness since depth ordering is preserved.
 */
#define B3D_DEPTH_OFFSET \
    (B3D_FAR_DISTANCE / (B3D_FAR_DISTANCE - B3D_NEAR_DISTANCE))
#define B3D_DEPTH_SCALE                       \
    ((B3D_NEAR_DISTANCE * B3D_FAR_DISTANCE) / \
     (B3D_FAR_DISTANCE - B3D_NEAR_DISTANCE))

/* Precomputed depth constants in scalar format.
 * Perspective mode: depth = offset - w_inv * scale (standard formula)
 * Orthographic mode: depth = z_ndc directly (offset=0, scale=-1)
 */
static b3d_scalar_t b3d_depth_offset_fp = 0;
static b3d_scalar_t b3d_depth_scale_fp = 0;

static inline void b3d_update_depth_constants(bool ortho)
{
    if (ortho) {
        /* Orthographic: d = 0 - z_ndc * (-1) = z_ndc */
        b3d_depth_offset_fp = 0;
        b3d_depth_scale_fp = B3D_FLOAT_TO_FP(-1.0f);
    } else {
        /* Perspective: d = offset - (1/w) * scale */
        b3d_depth_offset_fp = B3D_FLOAT_TO_FP((float) B3D_DEPTH_OFFSET);
        b3d_depth_scale_fp = B3D_FLOAT_TO_FP((float) B3D_DEPTH_SCALE);
    }
}

static inline void b3d_init_depth_constants(void)
{
    b3d_update_depth_constants(b3d_ortho_mode);
}

/* Cached screen-space clipping planes (updated when resolution changes) */
static b3d_vec_t b3d_screen_planes[4][2];
static int b3d_planes_cached_w = 0, b3d_planes_cached_h = 0;

#ifdef B3D_NO_CULLING
/* Update cached model*view matrix if dirty (lazy matrix computation) */
static void b3d_update_model_view(void)
{
    if (b3d_model_view_dirty) {
        b3d_model_view = b3d_mat_mul(b3d_view, b3d_model);
        b3d_model_view_dirty = false;
    }
}
#endif

static void b3d_update_screen_planes(void)
{
    if (b3d_planes_cached_w == b3d_width && b3d_planes_cached_h == b3d_height)
        return;

    /* Top edge */
    b3d_screen_planes[0][0] = (b3d_vec_t) {0, 0.5f, 0, 1};
    b3d_screen_planes[0][1] = (b3d_vec_t) {0, 1, 0, 1};
    /* Bottom edge */
    b3d_screen_planes[1][0] = (b3d_vec_t) {0, (float) b3d_height, 0, 1};
    b3d_screen_planes[1][1] = (b3d_vec_t) {0, -1, 0, 1};
    /* Left edge */
    b3d_screen_planes[2][0] = (b3d_vec_t) {0.5f, 0, 0, 1};
    b3d_screen_planes[2][1] = (b3d_vec_t) {1, 0, 0, 1};
    /* Right edge */
    b3d_screen_planes[3][0] = (b3d_vec_t) {(float) b3d_width, 0, 0, 1};
    b3d_screen_planes[3][1] = (b3d_vec_t) {-1, 0, 0, 1};

    b3d_planes_cached_w = b3d_width;
    b3d_planes_cached_h = b3d_height;
}

static inline b3d_scalar_t b3d_depth_load(b3d_depth_t v)
{
#ifdef B3D_DEPTH_16BIT
    /* Direct: uint16 [0, 65535] → fixed-point [0, B3D_FP_ONE]
     * Special-case max value to ensure 0xFFFF → B3D_FP_ONE exactly,
     * so fragments at depth 1.0 pass the far-plane depth test.
     */
    if (v == 0xFFFF)
        return (b3d_scalar_t) B3D_FP_ONE;
    /* Reciprocal multiply for other values: v/65535 ≈ (v * 65537) >> 16 */
    return (b3d_scalar_t) (((uint32_t) v * 65537U) >> 16);
#else
    return (b3d_scalar_t) v;
#endif
}

static inline b3d_depth_t b3d_depth_store(b3d_scalar_t v)
{
#ifdef B3D_DEPTH_16BIT
    /* Direct: fixed-point [0, B3D_FP_ONE] → uint16 [0, 65535]
     * Add 0.5 ulp for unbiased rounding before shift
     */
    if (v <= 0)
        return 0;
    if (v >= B3D_FP_ONE)
        return 0xFFFF;
    return (b3d_depth_t) ((((int64_t) v * 65535) + (1 << (B3D_FP_BITS - 1))) >>
                          B3D_FP_BITS);
#else
    return (b3d_depth_t) v;
#endif
}

static inline b3d_scalar_t b3d_fp_min(b3d_scalar_t a, b3d_scalar_t b)
{
    return a < b ? a : b;
}

static inline b3d_scalar_t b3d_fp_max(b3d_scalar_t a, b3d_scalar_t b)
{
    return a > b ? a : b;
}

#define B3D_FP_DEGEN_THRESHOLD B3D_FLOAT_TO_FP(B3D_DEGEN_THRESHOLD)

/* Swap two vertices (x, y, z components) using token pasting */
#define SWAP_VERTICES(a, b, tmp) \
    do {                         \
        tmp = a##x;              \
        a##x = b##x;             \
        b##x = tmp;              \
        tmp = a##y;              \
        a##y = b##y;             \
        b##y = tmp;              \
        tmp = a##z;              \
        a##z = b##z;             \
        b##z = tmp;              \
    } while (0)

/* Swap span endpoints (x and z) when start > end */
#define SWAP_SPAN(sx, sz, ex, ez, tmp) \
    do {                               \
        tmp = sx;                      \
        sx = ex;                       \
        ex = tmp;                      \
        tmp = sz;                      \
        sz = ez;                       \
        ez = tmp;                      \
    } while (0)

/* Transform all 3 vertices of a triangle by matrix */
#define TRANSFORM_TRI(tri, mat)                        \
    do {                                               \
        (tri).p[0] = b3d_mat_mul_vec(mat, (tri).p[0]); \
        (tri).p[1] = b3d_mat_mul_vec(mat, (tri).p[1]); \
        (tri).p[2] = b3d_mat_mul_vec(mat, (tri).p[2]); \
    } while (0)

/* Perspective divide for all 3 vertices */
#define PERSPECTIVE_DIV(tri)                                \
    do {                                                    \
        (tri).p[0] = b3d_vec_div((tri).p[0], (tri).p[0].w); \
        (tri).p[1] = b3d_vec_div((tri).p[1], (tri).p[1].w); \
        (tri).p[2] = b3d_vec_div((tri).p[2], (tri).p[2].w); \
    } while (0)

/* Convert NDC vertex to screen coordinates */
#define NDC_TO_SCREEN(p, xs, ys)     \
    do {                             \
        (p).x = ((p).x + 1) * (xs);  \
        (p).y = (-(p).y + 1) * (ys); \
    } while (0)

/* Clamp integer to [lo, hi] range - inline function avoids double evaluation */
static inline int b3d_clamp_int(int v, int lo, int hi)
{
    if (v < lo)
        return lo;
    if (v > hi)
        return hi;
    return v;
}

/* Pixel write macro for scanline unrolling.
 * Uses perspective-correct depth: depth = offset - w_inv * scale.
 * Expects: w = current 1/w, w_step = delta 1/w per pixel.
 * Uses cached b3d_depth_offset_fp and b3d_depth_scale_fp.
 * Clamps depth to [0, B3D_FP_ONE] for numerical stability.
 */
#define PUT_PIXEL(i)                                                 \
    do {                                                             \
        b3d_scalar_t d =                                             \
            b3d_depth_offset_fp - B3D_FP_MUL(w, b3d_depth_scale_fp); \
        if (d < 0)                                                   \
            d = 0;                                                   \
        else if (d > B3D_FP_ONE)                                     \
            d = B3D_FP_ONE;                                          \
        if (d < b3d_depth_load(dp[i])) {                             \
            dp[i] = b3d_depth_store(d);                              \
            pp[i] = c;                                               \
        }                                                            \
        w = B3D_FP_ADD(w, w_step);                                   \
    } while (0)

/* Edge interpolation state for rasterizer.
 * Uses 1/w for perspective-correct depth interpolation.
 * t and t_step use higher precision (Q12.20) for smoother edges.
 * Parametric interpolation: value = start + total_delta * t
 */
typedef struct {
    b3d_scalar_t x, w_inv;   /* start position: x coord and 1/w */
    b3d_scalar_t dx, dw_inv; /* total edge delta (end - start) */
    b3d_tri_scalar_t t;      /* interpolation parameter [0, 1] (Q12.20) */
    b3d_tri_scalar_t t_step; /* step per scanline (Q12.20) */
} raster_edge_t;

/* Rasterize one half of a triangle (top or bottom).
 * Left/right edges interpolate x and 1/w with parameter t.
 * Depth is computed per-pixel from interpolated 1/w for perspective
 * correctness. Uses cached fixed-point depth constants.
 */
static void raster_half(int y_start,
                        int y_end,
                        raster_edge_t *left,
                        raster_edge_t *right,
                        uint32_t c)
{
    b3d_scalar_t tmp = 0;

    /* Clamp y range to screen bounds, fast-forward edge parameters.
     * Uses saturating multiply-add to prevent overflow on large skips.
     */
    if (y_start < 0) {
        int skip = -y_start;
        left->t = b3d_tri_mul_add_sat(left->t, left->t_step, skip);
        right->t = b3d_tri_mul_add_sat(right->t, right->t_step, skip);
        y_start = 0;
    }
    if (y_end > b3d_height)
        y_end = b3d_height;
    if (y_start >= y_end)
        return;

    /* Initialize row base for iterative update (addition vs multiplication) */
    size_t row_base = (size_t) y_start * (size_t) b3d_width;

    for (int y = y_start; y < y_end; y++) {
        /* Interpolate x and 1/w along edges using higher-precision t (Q12.20)
         */
        b3d_scalar_t sx = left->x + B3D_FP_MUL_TRI(left->dx, left->t);
        b3d_scalar_t sw = left->w_inv + B3D_FP_MUL_TRI(left->dw_inv, left->t);
        b3d_scalar_t ex = right->x + B3D_FP_MUL_TRI(right->dx, right->t);
        b3d_scalar_t ew =
            right->w_inv + B3D_FP_MUL_TRI(right->dw_inv, right->t);
        if (sx > ex)
            SWAP_SPAN(sx, sw, ex, ew, tmp);
        b3d_scalar_t dx = ex - sx;
        if (dx < B3D_FP_DEGEN_THRESHOLD) {
            left->t += left->t_step;
            right->t += right->t_step;
            row_base += b3d_width;
            continue;
        }

        /* Compute 1/w step per pixel for perspective-correct interpolation.
         * Guard against overflow: if |dw|/dx ratio is too large, the result
         * won't fit in int32_t. This happens on extremely thin spans with
         * large depth range - skip these degenerate cases.
         * Use int64_t to avoid overflow in the comparison itself.
         */
        b3d_scalar_t dw = ew - sw;
        b3d_scalar_t dw_abs = dw < 0 ? -dw : dw;
        /* Max ratio ~32 ensures w_step fits in int32_t with margin */
        if ((int64_t) dw_abs > (int64_t) dx * 32) {
            left->t += left->t_step;
            right->t += right->t_step;
            row_base += b3d_width;
            continue;
        }
        b3d_scalar_t w_step = B3D_FP_DIV(dw, dx);

        int start = B3D_FP_TO_INT(sx), end = B3D_FP_TO_INT(ex);
        start = b3d_clamp_int(start, 0, b3d_width);
        end = b3d_clamp_int(end, 0, b3d_width);
        if (start >= end) {
            left->t += left->t_step;
            right->t += right->t_step;
            row_base += b3d_width;
            continue;
        }

        /* Compute starting 1/w at first visible pixel */
        b3d_scalar_t w = sw + B3D_FP_MUL(w_step, B3D_INT_TO_FP(start) - sx);

        b3d_depth_t *dp = b3d_depth + row_base + start;
        uint32_t *pp = b3d_pixels + row_base + start;
        int n = end - start;
        while (n >= 4) {
            PUT_PIXEL(0);
            PUT_PIXEL(1);
            PUT_PIXEL(2);
            PUT_PIXEL(3);
            dp += 4, pp += 4;
            n -= 4;
        }
        while (n-- > 0) {
            PUT_PIXEL(0);
            dp++, pp++;
        }

        left->t += left->t_step;
        right->t += right->t_step;
        row_base += b3d_width;
    }
}

/* Screen-space vertex for rasterization.
 * @x: screen-space x coordinate
 * @y: screen-space y coordinate
 * @w_inv: 1/w for perspective-correct depth interpolation (w = z_view)
 */
typedef struct {
    b3d_scalar_t x, y, w_inv;
} raster_vertex_t;

/* Internal rasterization function */
static void b3d_rasterize(const raster_vertex_t v[3], uint32_t c)
{
    /* Copy and floor vertices (w_inv is not floored, only screen coords) */
    raster_vertex_t a = {B3D_FP_FLOOR(v[0].x), B3D_FP_FLOOR(v[0].y),
                         v[0].w_inv};
    raster_vertex_t b = {B3D_FP_FLOOR(v[1].x), B3D_FP_FLOOR(v[1].y),
                         v[1].w_inv};
    raster_vertex_t cv = {B3D_FP_FLOOR(v[2].x), B3D_FP_FLOOR(v[2].y),
                          v[2].w_inv};

    /* Screen-space AABB early-out */
    b3d_scalar_t min_x = b3d_fp_min(b3d_fp_min(a.x, b.x), cv.x);
    b3d_scalar_t max_x = b3d_fp_max(b3d_fp_max(a.x, b.x), cv.x);
    b3d_scalar_t min_y = b3d_fp_min(b3d_fp_min(a.y, b.y), cv.y);
    b3d_scalar_t max_y = b3d_fp_max(b3d_fp_max(a.y, b.y), cv.y);
    if (max_x < 0 || min_x >= B3D_INT_TO_FP(b3d_width) || max_y < 0 ||
        min_y >= B3D_INT_TO_FP(b3d_height)) {
        return;
    }

    /* Sort vertices by Y coordinate (bubble sort for 3 elements) */
    raster_vertex_t tmp;
    if (a.y > b.y) {
        tmp = a;
        a = b;
        b = tmp;
    }
    if (a.y > cv.y) {
        tmp = a;
        a = cv;
        cv = tmp;
    }
    if (b.y > cv.y) {
        tmp = b;
        b = cv;
        cv = tmp;
    }

    /* Guard against degenerate triangles (division by zero) */
    b3d_scalar_t dy_total = cv.y - a.y;
    b3d_scalar_t dy_top = b.y - a.y;
    if (dy_total < B3D_FP_DEGEN_THRESHOLD)
        return;

    /* Setup left edge (A to C, spans entire triangle).
     * t_step uses Q12.20 for higher precision interpolation.
     */
    raster_edge_t left = {
        .x = a.x,
        .w_inv = a.w_inv,
        .dx = cv.x - a.x,
        .dw_inv = cv.w_inv - a.w_inv,
        .t = 0,
        .t_step = B3D_TRI_FP_DIV(B3D_TRI_FP_ONE, B3D_FP_TO_TRI(dy_total)),
    };

    /* Setup right edge for top half (A to B) */
    raster_edge_t right = {
        .x = a.x,
        .w_inv = a.w_inv,
        .dx = b.x - a.x,
        .dw_inv = b.w_inv - a.w_inv,
        .t = 0,
        .t_step = (dy_top > B3D_FP_DEGEN_THRESHOLD)
                      ? B3D_TRI_FP_DIV(B3D_TRI_FP_ONE, B3D_FP_TO_TRI(dy_top))
                      : 0,
    };

    /* Rasterize top half: right edge from A toward B */
    raster_half(B3D_FP_TO_INT(a.y), B3D_FP_TO_INT(b.y), &left, &right, c);

    /* Setup right edge for bottom half (B to C) */
    b3d_scalar_t dy_bot = cv.y - b.y;
    right = (raster_edge_t) {
        .x = b.x,
        .w_inv = b.w_inv,
        .dx = cv.x - b.x,
        .dw_inv = cv.w_inv - b.w_inv,
        .t = 0,
        .t_step = (dy_bot > B3D_FP_DEGEN_THRESHOLD)
                      ? B3D_TRI_FP_DIV(B3D_TRI_FP_ONE, B3D_FP_TO_TRI(dy_bot))
                      : 0,
    };

    /* Rasterize bottom half: right edge from B toward C */
    raster_half(B3D_FP_TO_INT(b.y), B3D_FP_TO_INT(cv.y), &left, &right, c);
}

#undef PUT_PIXEL

/* Input validation: reject triangles containing NaN or infinity values.
 * @tri: triangle to validate
 *
 * Returns true if all vertex coordinates are finite, false otherwise.
 */
static inline bool b3d_is_finite_tri(const b3d_tri_t *tri)
{
    for (int i = 0; i < 3; i++) {
        if (!isfinite(tri->v[i].x) || !isfinite(tri->v[i].y) ||
            !isfinite(tri->v[i].z))
            return false;
    }
    return true;
}

/* Public API */

bool b3d_triangle(const b3d_tri_t *tri, uint32_t c)
{
    if (!tri || !b3d_pixels || !b3d_depth)
        return false;
    if (!b3d_is_finite_tri(tri))
        return false;

    b3d_triangle_t t =
        (b3d_triangle_t) {{{tri->v[0].x, tri->v[0].y, tri->v[0].z, 1},
                           {tri->v[1].x, tri->v[1].y, tri->v[1].z, 1},
                           {tri->v[2].x, tri->v[2].y, tri->v[2].z, 1}}};
#ifdef B3D_NO_CULLING
    /* Lazy matrix computation: use cached model*view when culling disabled */
    b3d_update_model_view();
    TRANSFORM_TRI(t, b3d_model_view);
#else
    /* Standard path: separate transforms for world-space culling */
    TRANSFORM_TRI(t, b3d_model);
    b3d_vec_t line_a = b3d_vec_sub(t.p[1], t.p[0]);
    b3d_vec_t line_b = b3d_vec_sub(t.p[2], t.p[0]);
    b3d_vec_t normal = b3d_vec_cross(line_a, line_b);
    b3d_vec_t cam_ray = b3d_vec_sub(t.p[0], b3d_camera);
    bool cull = !(c & B3D_DRAW_BACKFACE);
    if (cull && b3d_vec_dot(normal, cam_ray) > B3D_CULL_THRESHOLD)
        return false;
    TRANSFORM_TRI(t, b3d_view);
#endif

    b3d_triangle_t clipped[2];
    int count = b3d_clip_against_plane((b3d_vec_t) {0, 0, B3D_NEAR_DISTANCE, 1},
                                       (b3d_vec_t) {0, 0, 1, 1}, t, clipped);
    if (count == 0)
        return false;

    /* Far-plane clipping: keep vertices where z <= B3D_FAR_DISTANCE */
    b3d_triangle_t view_clipped[4];
    int view_count = 0;
    for (int i = 0; i < count; ++i) {
        b3d_triangle_t fc[2];
        int n =
            b3d_clip_against_plane((b3d_vec_t) {0, 0, B3D_FAR_DISTANCE, 1},
                                   (b3d_vec_t) {0, 0, -1, 1}, clipped[i], fc);
        for (int j = 0; j < n && view_count < 4; ++j)
            view_clipped[view_count++] = fc[j];
    }
    if (view_count == 0)
        return false;

    b3d_triangle_t buf_a[B3D_CLIP_BUFFER_SIZE], buf_b[B3D_CLIP_BUFFER_SIZE];
    b3d_triangle_t *src = buf_a, *dst = buf_b;
    int src_count = 0;
    for (int n = 0; n < view_count; ++n) {
        t = view_clipped[n];
        TRANSFORM_TRI(t, b3d_proj);
        if (fabsf(t.p[0].w) < B3D_EPSILON || fabsf(t.p[1].w) < B3D_EPSILON ||
            fabsf(t.p[2].w) < B3D_EPSILON)
            continue;

        /* Depth value to interpolate in screen space.
         * Perspective: 1/w for perspective-correct interpolation
         * Orthographic: z_ndc directly (w=1, no perspective correction needed)
         */
        float depth0, depth1, depth2;
        if (b3d_ortho_mode) {
            /* In ortho, z is already normalized depth [0,1] from projection */
            depth0 = t.p[0].z;
            depth1 = t.p[1].z;
            depth2 = t.p[2].z;
        } else {
            /* In perspective, w = z_view; store 1/w for correct interpolation
             */
            depth0 = 1.0f / t.p[0].w;
            depth1 = 1.0f / t.p[1].w;
            depth2 = 1.0f / t.p[2].w;
        }

        PERSPECTIVE_DIV(t);

        /* Store depth value in z for screen-space clipping interpolation */
        t.p[0].z = depth0;
        t.p[1].z = depth1;
        t.p[2].z = depth2;

        float xs = b3d_width * 0.5f;
        float ys = b3d_height * 0.5f;
        NDC_TO_SCREEN(t.p[0], xs, ys);
        NDC_TO_SCREEN(t.p[1], xs, ys);
        NDC_TO_SCREEN(t.p[2], xs, ys);
        if (src_count < B3D_CLIP_BUFFER_SIZE)
            src[src_count++] = t;
        else
            ++b3d_clip_drop_count;
    }

    for (int p = 0; p < 4; ++p) {
        int dst_count = 0;
        for (int i = 0; i < src_count; ++i) {
            int n = b3d_clip_against_plane(b3d_screen_planes[p][0],
                                           b3d_screen_planes[p][1], src[i],
                                           clipped);
            for (int w = 0; w < n; ++w) {
                if (dst_count < B3D_CLIP_BUFFER_SIZE)
                    dst[dst_count++] = clipped[w];
                else
                    ++b3d_clip_drop_count;
            }
        }

        b3d_triangle_t *tmp = src;
        src = dst;
        dst = tmp;
        src_count = dst_count;
    }
    if (src_count == 0)
        return false;

    /* z component now contains 1/w for perspective-correct depth */
    for (int i = 0; i < src_count; ++i) {
        raster_vertex_t rv[3] = {
            {B3D_FLOAT_TO_FP(src[i].p[0].x), B3D_FLOAT_TO_FP(src[i].p[0].y),
             B3D_FLOAT_TO_FP(src[i].p[0].z)},
            {B3D_FLOAT_TO_FP(src[i].p[1].x), B3D_FLOAT_TO_FP(src[i].p[1].y),
             B3D_FLOAT_TO_FP(src[i].p[1].z)},
            {B3D_FLOAT_TO_FP(src[i].p[2].x), B3D_FLOAT_TO_FP(src[i].p[2].y),
             B3D_FLOAT_TO_FP(src[i].p[2].z)},
        };
        b3d_rasterize(rv, c & 0x00FFFFFFU);
    }
    return true;
}

void b3d_reset(void)
{
    b3d_model = b3d_mat_ident();
    b3d_model_view_dirty = true;
}

void b3d_rotate_x(float angle)
{
    b3d_model = b3d_mat_mul(b3d_model, b3d_mat_rot_x(angle));
    b3d_model_view_dirty = true;
}

void b3d_rotate_y(float angle)
{
    b3d_model = b3d_mat_mul(b3d_model, b3d_mat_rot_y(angle));
    b3d_model_view_dirty = true;
}

void b3d_rotate_z(float angle)
{
    b3d_model = b3d_mat_mul(b3d_model, b3d_mat_rot_z(angle));
    b3d_model_view_dirty = true;
}

void b3d_translate(float x, float y, float z)
{
    b3d_model = b3d_mat_mul(b3d_model, b3d_mat_trans(x, y, z));
    b3d_model_view_dirty = true;
}

void b3d_scale(float x, float y, float z)
{
    b3d_model = b3d_mat_mul(b3d_model, b3d_mat_scale(x, y, z));
    b3d_model_view_dirty = true;
}

void b3d_set_fov(float fov_in_degrees)
{
    if (b3d_width <= 0 || b3d_height <= 0)
        return;

    b3d_fov_degrees = fov_in_degrees;
    b3d_proj = b3d_mat_proj(fov_in_degrees, b3d_height / (float) b3d_width,
                            B3D_NEAR_DISTANCE, B3D_FAR_DISTANCE);

    /* Switch to perspective mode if currently orthographic */
    if (b3d_ortho_mode) {
        b3d_ortho_mode = false;
        b3d_update_depth_constants(false);
    }
}

void b3d_ortho(float left,
               float right,
               float bottom,
               float top,
               float near,
               float far)
{
    /* Validate parameters */
    if (!isfinite(left) || !isfinite(right) || !isfinite(bottom) ||
        !isfinite(top) || !isfinite(near) || !isfinite(far))
        return;
    if (fabsf(right - left) < B3D_EPSILON ||
        fabsf(top - bottom) < B3D_EPSILON || fabsf(far - near) < B3D_EPSILON)
        return;

    b3d_proj = b3d_mat_ortho(left, right, bottom, top, near, far);
    b3d_ortho_mode = true;
    b3d_update_depth_constants(true);
}

bool b3d_is_ortho(void)
{
    return b3d_ortho_mode;
}

void b3d_set_camera(const b3d_camera_t *cam)
{
    if (!cam)
        return;

    b3d_camera_params = *cam;
    b3d_camera = (b3d_vec_t) {cam->x, cam->y, cam->z, 1};
    b3d_vec_t up = {0, 1, 0, 1};
    b3d_vec_t target = {0, 0, 1, 1};
    up = b3d_mat_mul_vec(b3d_mat_rot_z(cam->roll), up);
    target = b3d_mat_mul_vec(b3d_mat_rot_x(cam->pitch), target);
    target = b3d_mat_mul_vec(b3d_mat_rot_y(cam->yaw), target);
    target = b3d_vec_add(b3d_camera, target);
    b3d_view = b3d_mat_qinv(b3d_mat_point_at(b3d_camera, target, up));
    b3d_model_view_dirty = true;
}

void b3d_look_at(float x, float y, float z)
{
    b3d_vec_t up = {0, 1, 0, 1};
    b3d_view = b3d_mat_qinv(
        b3d_mat_point_at(b3d_camera, (b3d_vec_t) {x, y, z, 1}, up));
    b3d_model_view_dirty = true;
}

bool b3d_to_screen(float x, float y, float z, int *sx, int *sy)
{
    if (!sx || !sy)
        return false;

    b3d_vec_t p = {x, y, z, 1};
    p = b3d_mat_mul_vec(b3d_model, p);
    p = b3d_mat_mul_vec(b3d_view, p);
    p = b3d_mat_mul_vec(b3d_proj, p);
    if (p.w < B3D_EPSILON)
        return false;

    p = b3d_vec_div(p, p.w);
    float mid_x = b3d_width / 2.0f;
    float mid_y = b3d_height / 2.0f;
    p.x = (p.x + 1.0f) * mid_x;
    p.y = (-p.y + 1.0f) * mid_y;
    *sx = (int) (p.x + 0.5f);
    *sy = (int) (p.y + 0.5f);
    return true;
}

bool b3d_init(uint32_t *pixel_buffer,
              b3d_depth_t *depth_buffer,
              int w,
              int h,
              float fov)
{
    size_t depth_bytes = b3d_buffer_size(w, h, sizeof(b3d_depth_t));
    size_t pixel_bytes = b3d_buffer_size(w, h, sizeof(uint32_t));
    if (!pixel_buffer || !depth_buffer || w <= 0 || h <= 0 || fov <= 0 ||
        depth_bytes == 0 || pixel_bytes == 0) {
        b3d_width = 0;
        b3d_height = 0;
        b3d_pixels = NULL;
        b3d_depth = NULL;
        b3d_matrix_stack_top = 0;
        b3d_clip_drop_count = 0;
        b3d_model_view_dirty = true;
        return false;
    }

    b3d_width = w;
    b3d_height = h;
    b3d_pixels = pixel_buffer;
    b3d_depth = depth_buffer;
    b3d_matrix_stack_top = 0;
    b3d_model_view_dirty = true;
    b3d_fov_degrees = fov;
    b3d_ortho_mode = false; /* Start in perspective mode */
    b3d_init_depth_constants();
    b3d_update_screen_planes();
    b3d_clear();
    b3d_reset();
    b3d_proj = b3d_mat_proj(fov, b3d_height / (float) b3d_width,
                            B3D_NEAR_DISTANCE, B3D_FAR_DISTANCE);
    b3d_set_camera(&(b3d_camera_t) {0, 0, 0, 0, 0, 0});
    return true;
}

void b3d_clear(void)
{
    if (!b3d_depth || !b3d_pixels || b3d_width <= 0 || b3d_height <= 0)
        return;

    b3d_clip_drop_count = 0;
    /* Check for integer overflow when calculating buffer size */
    if ((size_t) b3d_width > SIZE_MAX / (size_t) b3d_height)
        return; /* Prevent overflow */

    size_t count = (size_t) b3d_width * (size_t) b3d_height;
    /* Also check for overflow in the pixel buffer size calculation */
    if (count > SIZE_MAX / sizeof(b3d_pixels[0]))
        return; /* Prevent overflow */

    for (size_t i = 0; i < count; ++i)
        b3d_depth[i] = B3D_DEPTH_CLEAR;
    memset(b3d_pixels, 0, count * sizeof(b3d_pixels[0]));
}

size_t b3d_buffer_size(int w, int h, size_t elem_size)
{
    if (w <= 0 || h <= 0 || elem_size == 0)
        return 0;

    /* Check for overflow: w * h * elem_size */
    size_t sw = (size_t) w;
    size_t sh = (size_t) h;
    if (sw > SIZE_MAX / sh)
        return 0;

    size_t count = sw * sh;
    if (count > SIZE_MAX / elem_size)
        return 0;
    return count * elem_size;
}

bool b3d_push_matrix(void)
{
    if (b3d_matrix_stack_top >= B3D_MATRIX_STACK_SIZE)
        return false;
    b3d_matrix_stack[b3d_matrix_stack_top++] = b3d_model;
    return true;
}

bool b3d_pop_matrix(void)
{
    if (b3d_matrix_stack_top <= 0)
        return false;
    b3d_model = b3d_matrix_stack[--b3d_matrix_stack_top];
    b3d_model_view_dirty = true;
    return true;
}

void b3d_get_model_matrix(float out[16])
{
    if (!out)
        return;

    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c)
            out[r * 4 + c] = b3d_model.m[r][c];
    }
}

void b3d_set_model_matrix(const float m[16])
{
    if (!m)
        return;

    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c)
            b3d_model.m[r][c] = m[r * 4 + c];
    }
    b3d_model_view_dirty = true;
}

size_t b3d_get_clip_drop_count(void)
{
    return b3d_clip_drop_count;
}

void b3d_get_camera(b3d_camera_t *out)
{
    if (!out)
        return;
    *out = b3d_camera_params;
}

float b3d_get_fov(void)
{
    return b3d_fov_degrees;
}

void b3d_get_view_matrix(float out[16])
{
    if (!out)
        return;

    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c)
            out[r * 4 + c] = b3d_view.m[r][c];
    }
}

void b3d_get_proj_matrix(float out[16])
{
    if (!out)
        return;

    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c)
            out[r * 4 + c] = b3d_proj.m[r][c];
    }
}

bool b3d_is_initialized(void)
{
    return b3d_pixels != NULL && b3d_depth != NULL && b3d_width > 0 &&
           b3d_height > 0;
}

int b3d_get_width(void)
{
    return b3d_width;
}

int b3d_get_height(void)
{
    return b3d_height;
}

void b3d_set_light_direction(float x, float y, float z)
{
    /* Reject NaN/INF inputs */
    if (!isfinite(x) || !isfinite(y) || !isfinite(z))
        return;
    float len_sq = x * x + y * y + z * z;
    if (len_sq < B3D_EPSILON)
        return; /* Reject zero-length, keep previous direction */
    float inv_len = 1.0f / b3d_sqrtf(len_sq);
    b3d_light_dir = (b3d_vec_t) {x * inv_len, y * inv_len, z * inv_len, 0.0f};
}

void b3d_get_light_direction(float *x, float *y, float *z)
{
    if (x)
        *x = b3d_light_dir.x;
    if (y)
        *y = b3d_light_dir.y;
    if (z)
        *z = b3d_light_dir.z;
}

float b3d_get_ambient(void)
{
    return b3d_ambient;
}

void b3d_set_ambient(float ambient)
{
    /* Reject NaN/INF inputs */
    if (!isfinite(ambient))
        return;
    if (ambient < 0.0f)
        ambient = 0.0f;
    if (ambient > 1.0f)
        ambient = 1.0f;
    b3d_ambient = ambient;
}

static inline uint32_t b3d_shade_color(uint32_t c, float intensity)
{
    if (intensity < 0.0f)
        intensity = 0.0f;
    if (intensity > 1.0f)
        intensity = 1.0f;
    int r = (int) (((c >> 16) & 0xFF) * intensity);
    int g = (int) (((c >> 8) & 0xFF) * intensity);
    int b = (int) ((c & 0xFF) * intensity);
    return (uint32_t) ((r << 16) | (g << 8) | b);
}

bool b3d_triangle_lit(const b3d_tri_t *tri,
                      float nx,
                      float ny,
                      float nz,
                      uint32_t base_color)
{
    if (!tri || !b3d_is_finite_tri(tri))
        return false;
    if (!isfinite(nx) || !isfinite(ny) || !isfinite(nz))
        return false;

    /* Normalize the surface normal */
    b3d_vec_t n = {nx, ny, nz, 0.0f};
    n = b3d_vec_norm(n);

    /* Diffuse lighting with two-sided shading (abs of dot product) */
    float dot = b3d_vec_dot(n, b3d_light_dir);
    if (dot < 0.0f)
        dot = -dot;

    /* Preserve control flags (e.g., B3D_DRAW_BACKFACE) */
    uint32_t flags = base_color & 0xFF000000U;

    float intensity = b3d_ambient + (1.0f - b3d_ambient) * dot;
    uint32_t shaded = b3d_shade_color(base_color, intensity);

    return b3d_triangle(tri, shaded | flags);
}
