/*
 * B3D is freely redistributable under the MIT License. See the file
 * "LICENSE" for information on usage and redistribution of this file.
 */

/*
 * Wavefront .obj file loader
 *
 * A simple .obj file loader that extracts triangle data.
 * Only supports triangulated meshes (no quads or n-gons).
 */

#ifndef B3D_OBJ_H
#define B3D_OBJ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    float *triangles;   /* Array of triangle vertices (9 floats per triangle:
                           ax,ay,az,bx,by,bz,cx,cy,cz) */
    int triangle_count; /* Number of triangles */
    int vertex_count; /* Total number of vertex components (triangle_count * 9)
                       */
} b3d_mesh_t;

/* Load a mesh from an OBJ file.
 * @path : Path to the .obj file
 * @mesh : Pointer to mesh structure to fill
 *
 * Returns:
 *   0 on success, non-zero on error:
 *   1 = file not found
 *   2 = memory allocation failed
 *   3 = invalid vertex index in file
 *
 * Note: Only triangulated faces are supported. The mesh must be freed
 * with b3d_free_mesh() when no longer needed.
 */
static inline int b3d_load_obj(const char *path, b3d_mesh_t *mesh)
{
    if (!path || !mesh)
        return 1;

    mesh->triangles = NULL;
    mesh->triangle_count = 0;
    mesh->vertex_count = 0;

    FILE *obj_file = fopen(path, "r");
    if (!obj_file)
        return 1;

    int vi = 0, ti = 0;
    float *vertices = NULL;
    float *triangles = NULL;
    char line[1024];
    float x, y, z;

    /* First pass: read vertices */
    while (fgets(line, sizeof(line), obj_file)) {
        if (line[0] == 'v' && line[1] == ' ' &&
            sscanf(line, " v %f %f %f ", &x, &y, &z) == 3) {
            float *temp = realloc(vertices, (vi + 3) * sizeof(vertices[0]));
            if (!temp) {
                fclose(obj_file);
                free(vertices);
                return 2;
            }
            vertices = temp;
            vertices[vi++] = x;
            vertices[vi++] = y;
            vertices[vi++] = z;
        }
    }

    rewind(obj_file);

    /* Second pass: read faces */
    int a, b, c;
    while (fgets(line, sizeof(line), obj_file)) {
        if (line[0] == 'f') {
            /* Erase texture and normal indices (e.g., "1/2/3" -> "1") */
            for (size_t i = 0; i < sizeof(line) && line[i] != '\0'; ++i) {
                if (line[i] == '/') {
                    while (i < sizeof(line) && line[i] != '\0' &&
                           line[i] != ' ' && line[i] != '\n') {
                        line[i++] = ' ';
                    }
                }
            }
            if (sscanf(line, " f %d %d %d ", &a, &b, &c) == 3) {
                /* Validate vertex indices (OBJ uses 1-based indexing) */
                if (a <= 0 || b <= 0 || c <= 0 || a > vi / 3 || b > vi / 3 ||
                    c > vi / 3) {
                    fclose(obj_file);
                    free(vertices);
                    free(triangles);
                    return 3;
                }
                a--;
                b--;
                c--;
                float *temp =
                    realloc(triangles, (ti + 9) * sizeof(triangles[0]));
                if (!temp) {
                    fclose(obj_file);
                    free(vertices);
                    free(triangles);
                    return 2;
                }
                triangles = temp;
                triangles[ti++] = vertices[(a * 3) + 0];
                triangles[ti++] = vertices[(a * 3) + 1];
                triangles[ti++] = vertices[(a * 3) + 2];
                triangles[ti++] = vertices[(b * 3) + 0];
                triangles[ti++] = vertices[(b * 3) + 1];
                triangles[ti++] = vertices[(b * 3) + 2];
                triangles[ti++] = vertices[(c * 3) + 0];
                triangles[ti++] = vertices[(c * 3) + 1];
                triangles[ti++] = vertices[(c * 3) + 2];
            }
        }
    }

    fclose(obj_file);
    free(vertices);

    mesh->triangles = triangles;
    mesh->vertex_count = ti;
    mesh->triangle_count = ti / 9;

    return 0;
}

/* Free a mesh loaded with b3d_load_obj() */
static inline void b3d_free_mesh(b3d_mesh_t *mesh)
{
    if (!mesh)
        return;

    free(mesh->triangles);
    mesh->triangles = NULL;
    mesh->triangle_count = 0;
    mesh->vertex_count = 0;
}

/* Calculate mesh bounds (useful for centering and scaling).
 * @mesh   : The mesh to analyze
 * @min_y  : Output: minimum Y coordinate
 * @max_y  : Output: maximum Y coordinate
 * @max_xz : Output: maximum absolute X or Z coordinate
 */
static inline void b3d_mesh_bounds(const b3d_mesh_t *mesh,
                                   float *min_y,
                                   float *max_y,
                                   float *max_xz)
{
    if (!mesh || !mesh->triangles || mesh->vertex_count == 0)
        return;

    float miny = mesh->triangles[1];
    float maxy = mesh->triangles[1];
    float maxxz = 0;

    for (int i = 0; i < mesh->vertex_count; i += 3) {
        float y = mesh->triangles[i + 1];
        if (y < miny)
            miny = y;
        if (y > maxy)
            maxy = y;

        float ax = mesh->triangles[i + 0];
        float az = mesh->triangles[i + 2];
        if (ax < 0)
            ax = -ax;
        if (az < 0)
            az = -az;
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
