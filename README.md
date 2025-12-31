# B3D

`B3D` is a minimal software 3D renderer in C99,
derived from [bootleg3D](https://github.com/malespiaut/bootleg3d).

![fps demo](examples/fps.png)

## Features
- Triangle rasterization with depth buffering and clipping
- Perspective projection with configurable FOV
- Model transforms (translate, rotate, scale) and camera control
- Zero heap allocations; depends only on `<stdbool.h>`, `<stddef.h>`, `<stdint.h>`, and `<string.h>`
- Optional 16-bit depth buffer (`B3D_DEPTH_16BIT`)

## Quick Start

```c
#include "b3d.h"

#define W 320
#define H 240

int main(void) {
    static uint32_t pixels[W * H];
    static b3d_depth_t depth[W * H];

    b3d_init(pixels, depth, W, H, 60.0f);
    b3d_set_camera(&(b3d_camera_t){0, 0, -3, 0, 0, 0});

    b3d_clear();
    b3d_tri_t tri = {{{-1, -1, 0}, {1, -1, 0}, {0, 1, 0}}};
    b3d_triangle(&tri, 0xFF0000);  /* red triangle */

    /* pixels[] now contains rendered image */
    return 0;
}
```

## Usage

```c
#include "b3d.h"
```

Compile-time options:
- `B3D_DEPTH_32BIT` - Default, 32-bit fixed-point depth (16.16)
- `B3D_DEPTH_16BIT` - Use 16-bit depth buffer
- `B3D_FLOAT_POINT` - Use floating-point math for comparisons

## Examples

| Example | Description |
|---------|-------------|
| ball | Bouncing ball demo |
| cubes | Benchmark: adds cubes until ~60 fps |
| donut | Torus with directional lighting |
| fps | First-person treasure hunt |
| gears | Interlocking animated gears |
| globe | Rotating Earth with Natural Earth country data |
| lena3d | Image rendered as 3D heightfield |
| lighting | Directional lighting showcase |
| microcraft | Minecraft-style voxel terrain with procedural generation |
| obj | OBJ viewer with auto-centering |
| terrain | Animated sine/cosine heightmap |
| voxel | Mesh voxelizer with cube rendering |

Build with `make check` to validate B3D implementation.
Build with `make all` (requires SDL2). Run headlessly with `--snapshot=PATH`.

### Globe Data Generation

The globe example uses [Natural Earth](https://www.naturalearthdata.com/) country boundary data (public domain).

```bash
pip3 install pyshp requests   # one-time setup
make update-globe-data        # download and process data
make globe                    # build example
```

## API

```c
// Types
typedef struct { float x, y, z; } b3d_point_t;
typedef struct { b3d_point_t v[3]; } b3d_tri_t;
typedef struct { float x, y, z, yaw, pitch, roll; } b3d_camera_t;

// Setup
bool b3d_init(uint32_t *pixels, b3d_depth_t *depth, int w, int h, float fov);
void b3d_clear(void);
size_t b3d_buffer_size(int w, int h, size_t elem_size);  // 0 on overflow

// Transforms (angles in radians)
void b3d_reset(void);
void b3d_translate(float x, float y, float z);
void b3d_rotate_x(float angle);
void b3d_rotate_y(float angle);
void b3d_rotate_z(float angle);
void b3d_scale(float x, float y, float z);

// Matrix stack
bool b3d_push_matrix(void);  // false if stack full
bool b3d_pop_matrix(void);   // false if stack empty

// Direct matrix access (row-major 4x4)
void b3d_get_model_matrix(float out[16]);
void b3d_set_model_matrix(const float m[16]);

// Camera (yaw/pitch/roll in radians, FOV in degrees)
void b3d_set_camera(const b3d_camera_t *cam);
void b3d_get_camera(b3d_camera_t *out);
void b3d_look_at(float x, float y, float z);
void b3d_set_fov(float degrees);
float b3d_get_fov(void);
void b3d_get_view_matrix(float out[16]);
void b3d_get_proj_matrix(float out[16]);

// Rendering (returns false if culled/clipped)
bool b3d_triangle(const b3d_tri_t *tri, uint32_t color);
bool b3d_to_screen(float x, float y, float z, int *sx, int *sy);

// State queries
bool b3d_is_initialized(void);
int b3d_get_width(void);
int b3d_get_height(void);
size_t b3d_get_clip_drop_count(void);
```

## Matrix Conventions

B3D uses DirectX-style matrix conventions, which differ from OpenGL:

| Aspect | B3D | OpenGL |
|--------|-----|--------|
| Coordinate system | Right-handed, Y-up, +Z forward | Right-handed, Y-up, -Z forward |
| Memory layout | Row-major `m[row][col]` | Column-major `m[col][row]` |
| Vector type | Row vectors | Column vectors |
| Multiplication | `v * M` (vector on left) | `M * v` (vector on right) |
| API multiply order | Post-multiply: `M = M * T` | Post-multiply: `M = M * T` |
| Effect order | Transforms apply in call order | Transforms apply in reverse call order |
| Translation storage | Row 3: `m[3][0..2]` | Column 3: `m[12..14]` in flat array |

Transform application order:
```c
// B3D: transforms applied in code order (scale -> rotate -> translate)
b3d_scale(2, 2, 2);
b3d_rotate_y(angle);
b3d_translate(x, y, z);

// OpenGL: same API order, but effects apply in reverse (translate -> rotate -> scale)
// To match B3D's effect order, reverse the calls:
glTranslatef(x, y, z);
glRotatef(angle, 0, 1, 0);
glScalef(2, 2, 2);
```

Matrix layout comparison:
```
B3D (row-major):          OpenGL (column-major):
[ Xx Xy Xz 0 ]            [ Xx Yx Zx Tx ]
[ Yx Yy Yz 0 ]            [ Xy Yy Zy Ty ]
[ Zx Zy Zz 0 ]            [ Xz Yz Zz Tz ]
[ Tx Ty Tz 1 ]            [ 0  0  0  1  ]
```

When porting from OpenGL: transpose matrices and reverse transform order.

## Color Format

Colors are 32-bit RGB (alpha channel ignored):
- Format: `0xRRGGBB` or `0xAARRGGBB`
- Examples: `0xFF0000` (red), `0x00FF00` (green), `0x0000FF` (blue)

Per-face control flags use upper bits:
- `B3D_DRAW_BACKFACE` (bit 31): Disable backface culling for this triangle
- Usage: `b3d_triangle(&tri, color | B3D_DRAW_BACKFACE)`
- Use case: Double-sided geometry (leaves, flags), skyboxes with outward normals

## Coordinate System

| Property | Value |
|----------|-------|
| Handedness | Right-handed |
| Up axis | +Y |
| Forward axis | +Z (into screen) |
| Screen origin | Top-left (0, 0) |
| Winding order | Counter-clockwise = front face |
| Near plane | 0.1 units |
| Far plane | 100.0 units |

## Return Values

| Function | Returns `false` when |
|----------|---------------------|
| `b3d_init` | NULL buffers, invalid size, overflow, or FOV ≤ 0 |
| `b3d_push_matrix` | Stack full (16 levels max) |
| `b3d_pop_matrix` | Stack empty |
| `b3d_triangle` | Triangle culled or fully clipped |
| `b3d_to_screen` | Point behind camera (w < epsilon) |

Note: `b3d_look_at` updates the view matrix but leaves `b3d_get_camera` orientation stale.
Use `b3d_get_view_matrix` for accurate view state after `b3d_look_at`.

## Thread Safety

B3D uses global state and is NOT thread-safe. Do not call B3D functions
from multiple threads simultaneously. For multi-threaded rendering,
serialize all B3D calls or use separate processes.

## Performance

- Depth buffer: 16-bit mode (`B3D_DEPTH_16BIT`) halves memory bandwidth
- Fixed-point: Default Q15.16 is faster on systems without FPU
- Culling: Back-face culling is enabled by default; use `B3D_DRAW_BACKFACE`
  per-face for skyboxes/particles, or `B3D_NO_CULLING` compile-time for all
- Batching: Minimize `b3d_push_matrix`/`b3d_pop_matrix` pairs
- Clipping: Use `b3d_get_clip_drop_count()` to detect buffer overflow

## License
`B3D` is available under a permissive MIT-style license.
Use of this source code is governed by a MIT license that can be found in the [LICENSE](LICENSE) file.
