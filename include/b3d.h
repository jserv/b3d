/*
 * B3D is freely redistributable under the MIT License. See the file
 * "LICENSE" for information on usage and redistribution of this file.
 */

#ifndef B3D_H
#define B3D_H

#include <stdint.h>
#include <stddef.h>

/*
 * Optional compile-time configuration:
 *   B3D_NO_CULLING - Define to disable back-face culling
 *   B3D_DEPTH_16BIT      - Define to use 16-bit depth buffer (saves memory, less precision)
 */

/* Depth buffer element type (depends on compile-time option) */
#ifdef B3D_DEPTH_16BIT
typedef uint16_t b3d_depth_t;
#else
typedef float b3d_depth_t;
#endif
#define B3D_DEPTH_T_DEFINED

/*
 * Angle Unit Convention:
 *   - b3d_rotate_x/y/z(angle)        : angle in RADIANS
 *   - b3d_set_camera(yaw, pitch, roll): angles in RADIANS
 *   - b3d_set_fov(fov)               : fov in DEGREES
 *   - b3d_init(fov)                  : fov in DEGREES
 */

/* Matrix stack size for push/pop operations */
#define B3D_MATRIX_STACK_SIZE 16

/*
 * Initialization and clearing
 */

/* Initialize the renderer with pixel/depth buffers and field of view (degrees) */
void b3d_init(uint32_t *pixel_buffer, b3d_depth_t *depth_buffer, int w, int h, float fov);

/* Clear pixel buffer to black and depth buffer to far plane */
void b3d_clear(void);

/*
 * Calculate safe buffer size to avoid integer overflow.
 * Returns 0 if parameters would overflow.
 */
size_t b3d_buffer_size(int w, int h, size_t elem_size);

/*
 * Model transformation functions (angles in RADIANS)
 */

/* Reset model matrix to identity */
void b3d_reset(void);

/* Apply translation to model matrix */
void b3d_translate(float x, float y, float z);

/* Apply rotation around X axis (angle in radians) */
void b3d_rotate_x(float angle);

/* Apply rotation around Y axis (angle in radians) */
void b3d_rotate_y(float angle);

/* Apply rotation around Z axis (angle in radians) */
void b3d_rotate_z(float angle);

/* Apply scale to model matrix */
void b3d_scale(float x, float y, float z);

/*
 * Matrix stack operations for nested transformations
 */

/* Push current model matrix onto stack. Returns 1 on success, 0 if stack full. */
int b3d_push_matrix(void);

/* Pop model matrix from stack. Returns 1 on success, 0 if stack empty. */
int b3d_pop_matrix(void);

/*
 * Direct matrix access (row-major 4x4 matrix, 16 floats)
 */

/* Copy current model matrix to out[16] */
void b3d_get_model_matrix(float out[16]);

/* Set model matrix from m[16] */
void b3d_set_model_matrix(const float m[16]);

/*
 * Camera functions (angles in RADIANS except where noted)
 */

/* Set camera position and orientation (yaw, pitch, roll in radians) */
void b3d_set_camera(float x, float y, float z, float yaw, float pitch, float roll);

/* Point camera at target position (uses current camera position) */
void b3d_look_at(float x, float y, float z);

/* Set field of view (in DEGREES) */
void b3d_set_fov(float fov_in_degrees);

/*
 * Rendering
 */

/* Render a triangle. Returns 1 if rendered, 0 if culled/clipped away. */
int b3d_triangle(float ax, float ay, float az,
                 float bx, float by, float bz,
                 float cx, float cy, float cz,
                 uint32_t c);

/*
 * Utility functions
 */

/* Project world coordinate to screen coordinate. Returns 1 if in front of camera. */
int b3d_to_screen(float x, float y, float z, int *sx, int *sy);

/* Number of triangles dropped during clipping due to buffer limits (reset by b3d_clear). */
size_t b3d_get_clip_drop_count(void);

/*
 * Public globals (read-only recommended)
 */
extern int b3d_width, b3d_height;
extern uint32_t *b3d_pixels;
extern b3d_depth_t *b3d_depth;

#endif /* B3D_H */
