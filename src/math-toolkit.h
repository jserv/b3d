/*
 * B3D is freely redistributable under the MIT License. See the file
 * "LICENSE" for information on usage and redistribution of this file.
 */

/* Internal math types and functions. */

#ifndef B3D_MATH_H
#define B3D_MATH_H

#include <math.h>
#include <stdint.h>

/* Constants */
#define B3D_NEAR_DISTANCE 0.1f
#define B3D_FAR_DISTANCE 100.0f
#define B3D_EPSILON 1e-8f /* Near-zero threshold for division guards */
#define B3D_DEGEN_THRESHOLD                                                \
    0.0001f                      /* Degenerate triangle/scanline threshold \
                                  */
#define B3D_CULL_THRESHOLD 0.01f /* Back-face culling dot product threshold */
#define B3D_DEPTH_FAR 1e30f      /* Depth buffer clear value (far plane) */
#define B3D_PI 3.1415926536f     /* Pi constant for angle conversions */
#define B3D_CLIP_BUFFER_SIZE 32  /* Maximum triangles in clipping buffer */

/* Depth type (shared with public header) */
#ifndef B3D_DEPTH_T_DEFINED
#ifdef B3D_DEPTH_16BIT
typedef uint16_t b3d_depth_t;
#else
typedef float b3d_depth_t;
#endif
#define B3D_DEPTH_T_DEFINED
#endif

#ifdef B3D_DEPTH_16BIT
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
#define B3D_DEPTH_CLEAR B3D_DEPTH_FAR
static inline b3d_depth_t b3d_depth_from_float(float d)
{
    return d;
}
static inline float b3d_depth_to_float(b3d_depth_t d)
{
    return d;
}
#endif

/* Internal types */
typedef struct {
    float x, y, z, w;
} b3d_vec_t;

typedef struct {
    float m[4][4];
} b3d_mat_t;

typedef struct {
    b3d_vec_t p[3];
} b3d_triangle_t;

/* Vector operations */
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

/* Matrix operations */
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

static inline b3d_mat_t b3d_mat_qinv(b3d_mat_t m)
{
    b3d_mat_t o = (b3d_mat_t) {{
        [0][0] = m.m[0][0],
        [0][1] = m.m[1][0],
        [0][2] = m.m[2][0],
        [0][3] = 0,
        [1][0] = m.m[0][1],
        [1][1] = m.m[1][1],
        [1][2] = m.m[2][1],
        [1][3] = 0,
        [2][0] = m.m[0][2],
        [2][1] = m.m[1][2],
        [2][2] = m.m[2][2],
        [2][3] = 0,
    }};
    o.m[3][0] = -(m.m[3][0] * o.m[0][0] + m.m[3][1] * o.m[1][0] +
                  m.m[3][2] * o.m[2][0]);
    o.m[3][1] = -(m.m[3][0] * o.m[0][1] + m.m[3][1] * o.m[1][1] +
                  m.m[3][2] * o.m[2][1]);
    o.m[3][2] = -(m.m[3][0] * o.m[0][2] + m.m[3][1] * o.m[1][2] +
                  m.m[3][2] * o.m[2][2]);
    o.m[3][3] = 1;
    return o;
}

static inline b3d_mat_t b3d_mat_point_at(b3d_vec_t pos,
                                         b3d_vec_t target,
                                         b3d_vec_t up)
{
    b3d_vec_t forward = b3d_vec_sub(target, pos);
    forward = b3d_vec_norm(forward);
    b3d_vec_t a = b3d_vec_mul(forward, b3d_vec_dot(up, forward));
    up = b3d_vec_norm(b3d_vec_sub(up, a));
    b3d_vec_t right = b3d_vec_cross(up, forward);
    return (b3d_mat_t) {{
        [0][0] = right.x,
        [0][1] = right.y,
        [0][2] = right.z,
        [0][3] = 0,
        [1][0] = up.x,
        [1][1] = up.y,
        [1][2] = up.z,
        [1][3] = 0,
        [2][0] = forward.x,
        [2][1] = forward.y,
        [2][2] = forward.z,
        [2][3] = 0,
        [3][0] = pos.x,
        [3][1] = pos.y,
        [3][2] = pos.z,
        [3][3] = 1,
    }};
}

/* Geometry operations */
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
