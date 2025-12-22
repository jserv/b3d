# B3D

`B3D` is a minimal software 3D renderer in C99,
derived from [bootleg3D](https://github.com/malespiaut/bootleg3d).

![fps demo](examples/fps.png)

## Features
- Triangle rasterization with depth buffering and clipping
- Perspective projection with configurable FOV
- Model transforms (translate, rotate, scale) and camera control
- Zero heap allocations; depends only on `<stdint.h>` and `<string.h>`
- Optional 16-bit depth buffer (`B3D_DEPTH_16BIT`)

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
| cubes | Benchmark: adds cubes until ~60 fps |
| obj | OBJ viewer with auto-centering |
| fps | First-person treasure hunt |
| terrain | Animated sine/cosine heightmap |
| donut | Torus with directional lighting |
| lena3d | Image rendered as 3D heightfield |

Build with `make check` to validate B3D implementation.
Build with `make all` (requires SDL2). Run headlessly with `--snapshot=PATH`.

## API

```c
// Setup
void b3d_init(uint32_t *pixels, b3d_depth_t *depth, int w, int h, float fov);
void b3d_clear(void);

// Transforms (model matrix)
void b3d_reset(void);
void b3d_translate(float x, float y, float z);
void b3d_rotate_x(float angle);  // radians
void b3d_rotate_y(float angle);
void b3d_rotate_z(float angle);
void b3d_scale(float x, float y, float z);

// Camera
void b3d_set_camera(float x, float y, float z, float yaw, float pitch, float roll);
void b3d_look_at(float x, float y, float z);
void b3d_set_fov(float degrees);

// Rendering
void b3d_triangle(float ax, float ay, float az,
                  float bx, float by, float bz,
                  float cx, float cy, float cz, uint32_t color);
int b3d_to_screen(float x, float y, float z, int *sx, int *sy);

// Diagnostics
size_t b3d_get_clip_drop_count(void);
```

## License
`B3D` is available under a permissive MIT-style license.
Use of this source code is governed by a MIT license that can be found in the [LICENSE](LICENSE) file.
