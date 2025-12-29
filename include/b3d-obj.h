/*
 * B3D is freely redistributable under the MIT License. See the file
 * "LICENSE" for information on usage and redistribution of this file.
 */

/* Wavefront .obj file loader
 *
 * A simple .obj file loader that extracts triangle data.
 * Supports triangulated meshes and automatic fan triangulation for n-gons.
 * Based on parsing patterns from tinyobjloader-c.
 */

#ifndef B3D_OBJ_H
#define B3D_OBJ_H

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    float *triangles;   /* Triangle vertices: 9 floats per tri (ax,ay,az,...) */
    int triangle_count; /* Number of triangles */
    int vertex_count;   /* Total vertex components (triangle_count * 9) */
} b3d_mesh_t;

/* Internal: fast character classification */
#define B3D_OBJ_IS_SPACE(x) (((x) == ' ') || ((x) == '\t'))
#define B3D_OBJ_IS_DIGIT(x) ((unsigned) ((x) - '0') < 10u)
#define B3D_OBJ_IS_NEWLINE(x) (((x) == '\r') || ((x) == '\n') || ((x) == '\0'))

/* Internal: skip whitespace */
static inline void b3d_obj_skip_space(const char **p)
{
    while (B3D_OBJ_IS_SPACE(**p))
        (*p)++;
}

/* Internal: fast integer parser (handles sign) */
static inline int b3d_obj_parse_int(const char **p)
{
    int sign = 1, val = 0;

    b3d_obj_skip_space(p);
    if (**p == '-') {
        sign = -1;
        (*p)++;
    } else if (**p == '+') {
        (*p)++;
    }
    while (B3D_OBJ_IS_DIGIT(**p)) {
        val = val * 10 + (**p - '0');
        (*p)++;
    }
    return val * sign;
}

/* Internal: fast float parser */
static inline float b3d_obj_parse_float(const char **p)
{
    float sign = 1.0f, val = 0.0f, frac = 0.1f;
    int exp_sign = 1, exp = 0;

    b3d_obj_skip_space(p);
    if (**p == '-') {
        sign = -1.0f;
        (*p)++;
    } else if (**p == '+') {
        (*p)++;
    }

    /* Integer part */
    while (B3D_OBJ_IS_DIGIT(**p)) {
        val = val * 10.0f + (float) (**p - '0');
        (*p)++;
    }

    /* Fractional part */
    if (**p == '.') {
        (*p)++;
        while (B3D_OBJ_IS_DIGIT(**p)) {
            val += (float) (**p - '0') * frac;
            frac *= 0.1f;
            (*p)++;
        }
    }

    /* Exponent part */
    if (**p == 'e' || **p == 'E') {
        (*p)++;
        if (**p == '-') {
            exp_sign = -1;
            (*p)++;
        } else if (**p == '+') {
            (*p)++;
        }
        while (B3D_OBJ_IS_DIGIT(**p)) {
            exp = exp * 10 + (**p - '0');
            (*p)++;
        }
        /* Apply exponent */
        while (exp-- > 0)
            val = (exp_sign > 0) ? val * 10.0f : val * 0.1f;
    }

    return sign * val;
}

/* Internal: parse face vertex index (handles v, v/vt, v/vt/vn, v//vn) */
static inline int b3d_obj_parse_face_idx(const char **p, int num_vertices)
{
    b3d_obj_skip_space(p);
    int idx = b3d_obj_parse_int(p);

    /* Skip texture and normal indices if present */
    while (**p == '/') {
        (*p)++;
        /* Skip the number after slash */
        if (**p == '-' || **p == '+')
            (*p)++;
        while (B3D_OBJ_IS_DIGIT(**p))
            (*p)++;
    }

    /* Convert to 0-based, handle negative (relative) indices */
    if (idx > 0)
        return idx - 1;
    if (idx < 0)
        return num_vertices + idx;
    return 0; /* idx == 0 is invalid, treat as first vertex */
}

/* Internal: dynamic array growth factor */
#define B3D_OBJ_GROW_CAPACITY(cap) ((cap) < 8 ? 8 : (cap) * 2)

/* Load a mesh from an OBJ file.
 * @path:    path to the .obj file
 * @mesh:    pointer to mesh structure to fill
 *
 * Returns 0 on success, non-zero on error (1 = file not found,
 *         2 = memory allocation failed, 3 = invalid vertex index).
 * Supports automatic fan triangulation of n-gons. Free with b3d_free_mesh().
 */
static inline int b3d_load_obj(const char *path, b3d_mesh_t *mesh)
{
    char line[4096];
    int ret = 0;
    FILE *f = NULL;

    /* Vertex storage with capacity tracking */
    float *verts = NULL;
    int vert_count = 0, vert_cap = 0;

    /* Triangle storage with capacity tracking */
    float *tris = NULL;
    int tri_count = 0, tri_cap = 0;

    /* Face indices buffer for n-gon triangulation (max 64 vertices per face) */
    int face_idx[64];
    int face_count, face_truncated = 0;

    if (!path || !mesh)
        return 1;

    mesh->triangles = NULL;
    mesh->triangle_count = 0;
    mesh->vertex_count = 0;

    f = fopen(path, "r");
    if (!f)
        return 1;

    while (fgets(line, sizeof(line), f)) {
        const char *p = line;
        b3d_obj_skip_space(&p);

        /* Skip empty lines and comments */
        if (*p == '\0' || *p == '#' || *p == '\n' || *p == '\r')
            continue;

        /* Vertex: v x y z */
        if (p[0] == 'v' && B3D_OBJ_IS_SPACE(p[1])) {
            float x, y, z;
            p += 2;

            x = b3d_obj_parse_float(&p);
            y = b3d_obj_parse_float(&p);
            z = b3d_obj_parse_float(&p);

            /* Grow vertex buffer if needed */
            if (vert_count + 3 > vert_cap) {
                int new_cap = B3D_OBJ_GROW_CAPACITY(vert_cap);
                if (new_cap < vert_count + 3)
                    new_cap = vert_count + 3;
                if (new_cap > INT_MAX / (int) sizeof(float)) {
                    ret = 2;
                    goto cleanup;
                }
                float *tmp = realloc(verts, (size_t) new_cap * sizeof(float));
                if (!tmp) {
                    ret = 2;
                    goto cleanup;
                }
                verts = tmp;
                vert_cap = new_cap;
            }

            verts[vert_count++] = x;
            verts[vert_count++] = y;
            verts[vert_count++] = z;
            continue;
        }

        /* Face: f v1 v2 v3 ... (supports n-gons with fan triangulation) */
        if (p[0] == 'f' && B3D_OBJ_IS_SPACE(p[1])) {
            int i, num_verts = vert_count / 3;
            p += 2;

            /* Parse all face indices */
            face_count = 0;
            while (!B3D_OBJ_IS_NEWLINE(*p)) {
                b3d_obj_skip_space(&p);
                if (B3D_OBJ_IS_NEWLINE(*p))
                    break;

                if (face_count >= 64) {
                    face_truncated = 1;
                    while (!B3D_OBJ_IS_NEWLINE(*p))
                        p++;
                    break;
                }

                face_idx[face_count] = b3d_obj_parse_face_idx(&p, num_verts);

                if (face_idx[face_count] < 0 ||
                    face_idx[face_count] >= num_verts) {
                    ret = 3;
                    goto cleanup;
                }
                face_count++;
            }

            if (face_count < 3)
                continue;

            /* Fan triangulation: (0,1,2), (0,2,3), (0,3,4), ... */
            for (i = 2; i < face_count; i++) {
                int v0 = face_idx[0];
                int v1 = face_idx[i - 1];
                int v2 = face_idx[i];

                if (tri_count + 9 > tri_cap) {
                    int new_cap = B3D_OBJ_GROW_CAPACITY(tri_cap);
                    if (new_cap < tri_count + 9)
                        new_cap = tri_count + 9;
                    if (new_cap > INT_MAX / (int) sizeof(float)) {
                        ret = 2;
                        goto cleanup;
                    }
                    float *tmp =
                        realloc(tris, (size_t) new_cap * sizeof(float));
                    if (!tmp) {
                        ret = 2;
                        goto cleanup;
                    }
                    tris = tmp;
                    tri_cap = new_cap;
                }

                tris[tri_count++] = verts[v0 * 3 + 0];
                tris[tri_count++] = verts[v0 * 3 + 1];
                tris[tri_count++] = verts[v0 * 3 + 2];
                tris[tri_count++] = verts[v1 * 3 + 0];
                tris[tri_count++] = verts[v1 * 3 + 1];
                tris[tri_count++] = verts[v1 * 3 + 2];
                tris[tri_count++] = verts[v2 * 3 + 0];
                tris[tri_count++] = verts[v2 * 3 + 1];
                tris[tri_count++] = verts[v2 * 3 + 2];
            }
            continue;
        }
    }

    /* Shrink triangle buffer to exact size */
    if (tri_count > 0 && tri_count < tri_cap) {
        float *tmp = realloc(tris, (size_t) tri_count * sizeof(float));
        if (tmp)
            tris = tmp;
    }

    mesh->triangles = tris;
    mesh->vertex_count = tri_count;
    mesh->triangle_count = tri_count / 9;
    tris = NULL; /* Ownership transferred to mesh */

cleanup:
    (void) face_truncated;
    fclose(f);
    free(verts);
    free(tris);
    return ret;
}

/* Free a mesh loaded with b3d_load_obj() */
static inline void b3d_free_mesh(b3d_mesh_t *mesh)
{
    if (!mesh)
        return;

    if (mesh->triangles) {
        free(mesh->triangles);
        mesh->triangles = NULL;
    }
    mesh->triangle_count = 0;
    mesh->vertex_count = 0;
}

/* Calculate mesh bounds (useful for centering and scaling).
 * @mesh:       The mesh to analyze
 * @min_y:      Output: minimum Y coordinate
 * @max_y:      Output: maximum Y coordinate
 * @max_xz:     Output: maximum absolute X or Z coordinate
 */
static inline void b3d_mesh_bounds(const b3d_mesh_t *mesh,
                                   float *min_y,
                                   float *max_y,
                                   float *max_xz)
{
    if (!mesh || !mesh->triangles || mesh->vertex_count < 3) {
        if (min_y)
            *min_y = 0.0f;
        if (max_y)
            *max_y = 0.0f;
        if (max_xz)
            *max_xz = 0.0f;
        return;
    }

    float miny = mesh->triangles[1], maxy = mesh->triangles[1];
    float maxxz = 0;

    for (int i = 0; i < mesh->vertex_count; i += 3) {
        float x = mesh->triangles[i + 0];
        float y = mesh->triangles[i + 1];
        float z = mesh->triangles[i + 2];
        float ax = (x < 0) ? -x : x;
        float az = (z < 0) ? -z : z;

        if (y < miny)
            miny = y;
        if (y > maxy)
            maxy = y;
        if (ax > maxxz)
            maxxz = ax;
        if (az > maxxz)
            maxxz = az;
    }

    if (min_y)
        *min_y = miny;
    if (max_y)
        *max_y = maxy;
    if (max_xz)
        *max_xz = maxxz;
}

#endif /* B3D_OBJ_H */
