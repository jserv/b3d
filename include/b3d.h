/*
 * B3D is freely redistributable under the MIT License. See the file
 * "LICENSE" for information on usage and redistribution of this file.
 */

#ifndef B3D_H
#define B3D_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* Optional compile-time configuration:
 *   B3D_NO_CULLING   : Disable back-face culling
 *   B3D_DEPTH_32BIT  : Default, 32-bit fixed-point depth (16.16)
 *   B3D_DEPTH_16BIT  : Optional, 16-bit depth buffer
 *   B3D_FLOAT_POINT  : Use floating-point math for comparisons
 */

#ifndef B3D_DEPTH_16BIT
#ifndef B3D_DEPTH_32BIT
#define B3D_DEPTH_32BIT 1
#endif
#endif

/* Depth buffer element type (depends on compile-time option) */
#ifdef B3D_FLOAT_POINT
typedef float b3d_depth_t;
#elif defined(B3D_DEPTH_16BIT)
typedef uint16_t b3d_depth_t;
#else
typedef int32_t b3d_depth_t;
#endif
#define B3D_DEPTH_T_DEFINED

/* Angle Unit Convention:
 *   b3d_rotate_x/y/z(angle) : angle in RADIANS
 *   b3d_set_camera(yaw, pitch, roll): angles in RADIANS
 *   b3d_set_fov(fov) : fov in DEGREES
 *   b3d_init(fov) : fov in DEGREES
 */

/* Matrix stack size for push/pop operations */
#define B3D_MATRIX_STACK_SIZE 16

/* Public API types for efficient parameter passing */

/* 3D point/vertex */
typedef struct {
    float x, y, z;
} b3d_point_t;

/* Triangle defined by 3 vertices */
typedef struct {
    b3d_point_t v[3];
} b3d_tri_t;

/* Camera parameters */
typedef struct {
    float x, y, z;          /* position */
    float yaw, pitch, roll; /* orientation in radians */
} b3d_camera_t;

/* Initialization and clearing */

/* Initialize the renderer with pixel/depth buffers and field of view (degrees).
 * @pixel_buffer: output pixel buffer (w * h * sizeof(uint32_t) bytes)
 * @depth_buffer: depth buffer (w * h * sizeof(b3d_depth_t) bytes)
 * @w:            buffer width in pixels
 * @h:            buffer height in pixels
 * @fov:          field of view in degrees
 *
 * Returns true on success, false on failure (invalid parameters or overflow).
 */
bool b3d_init(uint32_t *pixel_buffer,
              b3d_depth_t *depth_buffer,
              int w,
              int h,
              float fov);

/* Clear pixel buffer to black and depth buffer to far plane */
void b3d_clear(void);

/* Calculate safe buffer size to avoid integer overflow.
 * @w:         buffer width in pixels
 * @h:         buffer height in pixels
 * @elem_size: size of each element in bytes
 *
 * Returns 0 if parameters would overflow.
 */
size_t b3d_buffer_size(int w, int h, size_t elem_size);

/* Model transformation functions (angles in RADIANS) */

/* Reset model matrix to identity */
void b3d_reset(void);

/* Apply translation to model matrix
 * @x, @y, @z: translation vector components
 */
void b3d_translate(float x, float y, float z);

/* Apply rotation around X axis (angle in radians)
 * @angle: rotation angle in radians
 */
void b3d_rotate_x(float angle);

/* Apply rotation around Y axis (angle in radians)
 * @angle: rotation angle in radians
 */
void b3d_rotate_y(float angle);

/* Apply rotation around Z axis (angle in radians)
 * @angle: rotation angle in radians
 */
void b3d_rotate_z(float angle);

/* Apply scale to model matrix
 * @x, @y, @z: scale factors for each axis
 */
void b3d_scale(float x, float y, float z);

/* Matrix stack operations for nested transformations */

/* Push current model matrix onto stack.
 * Returns true on success, false if stack full.
 */
bool b3d_push_matrix(void);

/* Pop model matrix from stack.
 * Returns true on success, false if stack empty.
 */
bool b3d_pop_matrix(void);

/* Direct matrix access (row-major 4x4 matrix, 16 floats) */

/* Copy current model matrix to out[16]
 * @out: output array of 16 floats for matrix data
 */
void b3d_get_model_matrix(float out[16]);

/* Set model matrix from m[16]
 * @m: input array of 16 floats for matrix data
 */
void b3d_set_model_matrix(const float m[16]);

/* Camera functions (angles in RADIANS except where noted) */

/* Set camera position and orientation (yaw, pitch, roll in radians)
 * @cam: pointer to camera parameters
 */
void b3d_set_camera(const b3d_camera_t *cam);

/* Point camera at target position (uses current camera position).
 * Note: Updates view matrix but invalidates yaw/pitch/roll in b3d_get_camera.
 * @x, @y, @z: target position to look at
 */
void b3d_look_at(float x, float y, float z);

/* Set field of view (in DEGREES)
 * @fov_in_degrees: field of view angle in degrees
 */
void b3d_set_fov(float fov_in_degrees);

/* Get current camera parameters
 * @out: output camera parameters (position and orientation)
 */
void b3d_get_camera(b3d_camera_t *out);

/* Get current field of view in degrees */
float b3d_get_fov(void);

/* Get current view matrix (row-major 4x4)
 * @out: output array of 16 floats
 */
void b3d_get_view_matrix(float out[16]);

/* Get current projection matrix (row-major 4x4)
 * @out: output array of 16 floats
 */
void b3d_get_proj_matrix(float out[16]);

/* Lighting
 *
 * Coordinate space: Light direction and triangle normals are both expected
 * in MODEL SPACE. Lighting is computed before model-view transformation,
 * so shading rotates with the object (classic glxgears style).
 */

/* Set light direction (auto-normalized).
 * @x, @y, @z: light direction vector components (model space)
 *
 * Zero-length and NaN/INF vectors are rejected (previous direction kept).
 * Default: (0, 0, 1) pointing toward +Z.
 */
void b3d_set_light_direction(float x, float y, float z);

/* Get current light direction (normalized).
 * @x, @y, @z: output pointers (NULL-safe)
 */
void b3d_get_light_direction(float *x, float *y, float *z);

/* Set ambient light level.
 * @ambient: ambient intensity clamped to [0, 1]
 *
 * NaN/INF values are rejected (previous level kept).
 * Default: 0.2 (20% ambient).
 */
void b3d_set_ambient(float ambient);

/* Get current ambient light level.
 * Returns ambient intensity in [0, 1].
 */
float b3d_get_ambient(void);

/* Rendering */

/* Render a triangle.
 * @tri: pointer to triangle with 3 vertices
 * @c:   triangle color in 0xRRGGBB format
 *
 * Returns true if rendered, false if culled/clipped away.
 */
bool b3d_triangle(const b3d_tri_t *tri, uint32_t c);

/* Render a lit triangle with surface normal.
 * @tri:        pointer to triangle with 3 vertices
 * @nx, @ny, @nz: surface normal in model space (will be normalized)
 * @base_color: base color in 0xRRGGBB format
 *
 * Uses two-sided diffuse lighting with ambient. Lighting is computed in
 * model space before transformation (shading rotates with the object).
 * Returns true if rendered, false if culled/clipped away.
 */
bool b3d_triangle_lit(const b3d_tri_t *tri,
                      float nx,
                      float ny,
                      float nz,
                      uint32_t base_color);

/* Utility functions */

/* Project world coordinate to screen coordinate.
 * @x, @y, @z:    world coordinates to project
 * @sx, @sy:      output screen coordinates
 *
 * Returns true if in front of camera, false otherwise.
 */
bool b3d_to_screen(float x, float y, float z, int *sx, int *sy);

/* Number of triangles dropped during clipping due to buffer limits (reset by
 * b3d_clear).
 * Returns number of triangles dropped since last clear.
 */
size_t b3d_get_clip_drop_count(void);

/* State query functions */

/* Check if renderer is initialized and ready */
bool b3d_is_initialized(void);

/* Get current framebuffer width in pixels */
int b3d_get_width(void);

/* Get current framebuffer height in pixels */
int b3d_get_height(void);

#endif /* B3D_H */
