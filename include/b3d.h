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
 * @x, @y, @z:          camera position in world coordinates
 * @yaw, @pitch, @roll: camera orientation angles in radians
 */
void b3d_set_camera(float x,
                    float y,
                    float z,
                    float yaw,
                    float pitch,
                    float roll);

/* Point camera at target position (uses current camera position)
 * @x, @y, @z: target position to look at
 */
void b3d_look_at(float x, float y, float z);

/* Set field of view (in DEGREES)
 * @fov_in_degrees: field of view angle in degrees
 */
void b3d_set_fov(float fov_in_degrees);

/* Rendering */

/* Render a triangle.
 * @ax, @ay, @az:    first vertex coordinates
 * @bx, @by, @bz:    second vertex coordinates
 * @cx, @cy, @cz:    third vertex coordinates
 * @c:               triangle color in 0xRRGGBB format
 *
 * Returns true if rendered, false if culled/clipped away.
 */
bool b3d_triangle(float ax,
                  float ay,
                  float az,
                  float bx,
                  float by,
                  float bz,
                  float cx,
                  float cy,
                  float cz,
                  uint32_t c);

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

#endif /* B3D_H */
