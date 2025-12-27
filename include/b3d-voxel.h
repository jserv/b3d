/*
 * B3D is freely redistributable under the MIT License. See the file
 * "LICENSE" for information on usage and redistribution of this file.
 */

/*
 * Mesh voxelizer for B3D
 *
 * Converts triangle meshes to voxel representations (axis-aligned cubes).
 * Based on voxelizer by Karim Naaji (MIT License).
 *
 * Reference:
 *   https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox2.txt
 */

#ifndef B3D_VOXEL_H
#define B3D_VOXEL_H

#include <stddef.h>
#include <stdint.h>

/* Voxel grid position */
typedef struct {
    float x, y, z;
} b3d_voxel_pos_t;

/* Voxel data: position and color */
typedef struct {
    b3d_voxel_pos_t pos;
    uint32_t color;
} b3d_voxel_t;

/* Axis-aligned bounding box */
typedef struct {
    b3d_voxel_pos_t min;
    b3d_voxel_pos_t max;
} b3d_aabb_t;

/*
 * Voxelize a triangle mesh.
 *
 * @triangles:   input triangle data (9 floats per triangle: ax,ay,az,bx,...)
 * @tri_count:   number of triangles
 * @voxel_size:  size of each voxel cube (must be > 0)
 * @color:       color for all voxels (0xRRGGBB)
 * @out_voxels:  output buffer for voxel data (caller-allocated), or NULL
 * @max_voxels:  maximum voxels the buffer can hold
 *
 * Returns:
 *   - Number of unique voxels generated when out_voxels is provided
 *   - max_voxels + 1 if output buffer was full (truncation indicator)
 *   - Upper bound estimate if out_voxels is NULL (for buffer sizing)
 *   - 0 on error (invalid input or voxel_size <= 0)
 *
 * Note: Estimation mode (NULL) may overestimate due to counting overlapping
 * triangle contributions without deduplication. Allocate 1.5-2x the estimate.
 */
size_t b3d_voxelize(const float *triangles,
                    size_t tri_count,
                    float voxel_size,
                    uint32_t color,
                    b3d_voxel_t *out_voxels,
                    size_t max_voxels);

/*
 * Render voxels as cubes using B3D (lit).
 *
 * @voxels:      array of voxels to render
 * @count:       number of voxels
 * @voxel_size:  size of each voxel cube (must be > 0)
 *
 * Renders each voxel as a cube with 12 triangles (6 faces x 2 triangles).
 * Uses b3d_triangle_lit for lighting calculations.
 * Silently returns if voxels is NULL, count is 0, or voxel_size <= 0.
 */
void b3d_voxel_render(const b3d_voxel_t *voxels, size_t count, float voxel_size);

/*
 * Render voxels as cubes using B3D (flat/unlit).
 *
 * @voxels:      array of voxels to render
 * @count:       number of voxels
 * @voxel_size:  size of each voxel cube (must be > 0)
 *
 * Renders each voxel as a cube without lighting calculations.
 * Silently returns if voxels is NULL, count is 0, or voxel_size <= 0.
 */
void b3d_voxel_render_flat(const b3d_voxel_t *voxels,
                           size_t count,
                           float voxel_size);

/*
 * Calculate mesh bounding box.
 *
 * @triangles: input triangle data (9 floats per triangle)
 * @tri_count: number of triangles
 * @aabb:      output bounding box
 */
void b3d_voxel_mesh_aabb(const float *triangles,
                         size_t tri_count,
                         b3d_aabb_t *aabb);

#endif /* B3D_VOXEL_H */
