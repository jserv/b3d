/*
 * Orthographic vs Perspective Projection Demo
 *
 * Demonstrates the visual difference between orthographic and perspective
 * projection modes. Press 'O' to toggle between modes.
 *
 * Key differences:
 *   Perspective:   Objects shrink with distance (foreshortening)
 *                  Parallel lines converge to vanishing points
 *                  Natural human vision simulation
 *
 *   Orthographic:  Objects retain size regardless of distance
 *                  Parallel lines stay parallel (no vanishing points)
 *                  Essential for CAD, 2D games, UI, isometric views
 *
 * Controls:
 *   O     - Toggle between orthographic and perspective
 *   ESC   - Quit
 *
 * Supports headless snapshots with '--snapshot=PATH' or B3D_SNAPSHOT.
 */

#include <SDL.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "b3d-math.h"
#include "b3d.h"
#include "utils.h"

/* Scene configuration */
#define CUBE_COUNT 5
#define CUBE_SIZE 1.0f
#define CUBE_SPACING_Z 3.0f
#define CUBE_SPACING_X 2.0f
#define GRID_EXTENT 6.0f
#define GRID_LINES 9
#define GRID_Y -1.5f
#define LINE_WIDTH 0.03f

static void print_mode_hint(bool use_ortho)
{
    if (use_ortho) {
        printf(
            "Mode: Orthographic (parallel lines stay parallel; distant cubes "
            "keep size)\n");
    } else {
        printf(
            "Mode: Perspective (parallel lines converge; distant cubes "
            "shrink)\n");
    }
}

static void update_window_title(SDL_Window *window, bool use_ortho)
{
    const char *mode = use_ortho ? "Orthographic" : "Perspective";
    const char *hint =
        use_ortho ? "parallel lines stay parallel" : "parallel lines converge";
    char title[96];
    snprintf(title, sizeof(title), "%s - %s (press O to toggle)", mode, hint);
    SDL_SetWindowTitle(window, title);
}

static void cleanup_sdl(SDL_Texture *texture,
                        SDL_Renderer *renderer,
                        SDL_Window *window)
{
    if (texture)
        SDL_DestroyTexture(texture);
    if (renderer)
        SDL_DestroyRenderer(renderer);
    if (window)
        SDL_DestroyWindow(window);
    SDL_Quit();
}

/* Cube face definitions (12 triangles = 6 faces x 2 triangles) */
static const b3d_tri_t cube_faces[12] = {
    /* Front face (-Z) */
    {
        {{-0.5f, -0.5f, -0.5f}, {-0.5f, 0.5f, -0.5f}, {0.5f, 0.5f, -0.5f}},
    },
    {
        {{-0.5f, -0.5f, -0.5f}, {0.5f, 0.5f, -0.5f}, {0.5f, -0.5f, -0.5f}},
    },
    /* Right face (+X) */
    {
        {{0.5f, -0.5f, -0.5f}, {0.5f, 0.5f, -0.5f}, {0.5f, 0.5f, 0.5f}},
    },
    {
        {{0.5f, -0.5f, -0.5f}, {0.5f, 0.5f, 0.5f}, {0.5f, -0.5f, 0.5f}},
    },
    /* Back face (+Z) */
    {
        {{0.5f, -0.5f, 0.5f}, {0.5f, 0.5f, 0.5f}, {-0.5f, 0.5f, 0.5f}},
    },
    {
        {{0.5f, -0.5f, 0.5f}, {-0.5f, 0.5f, 0.5f}, {-0.5f, -0.5f, 0.5f}},
    },
    /* Left face (-X) */
    {
        {{-0.5f, -0.5f, 0.5f}, {-0.5f, 0.5f, 0.5f}, {-0.5f, 0.5f, -0.5f}},
    },
    {
        {{-0.5f, -0.5f, 0.5f}, {-0.5f, 0.5f, -0.5f}, {-0.5f, -0.5f, -0.5f}},
    },
    /* Top face (+Y) */
    {
        {{-0.5f, 0.5f, -0.5f}, {-0.5f, 0.5f, 0.5f}, {0.5f, 0.5f, 0.5f}},
    },
    {
        {{-0.5f, 0.5f, -0.5f}, {0.5f, 0.5f, 0.5f}, {0.5f, 0.5f, -0.5f}},
    },
    /* Bottom face (-Y) */
    {{{0.5f, -0.5f, 0.5f}, {-0.5f, -0.5f, 0.5f}, {-0.5f, -0.5f, -0.5f}}},
    {
        {{0.5f, -0.5f, 0.5f}, {-0.5f, -0.5f, -0.5f}, {0.5f, -0.5f, -0.5f}},
    },
};

/* Face normals for lighting (one per face, used for both triangles) */
static const float cube_normals[6][3] = {
    {0, 0, -1}, /* Front */
    {1, 0, 0},  /* Right */
    {0, 0, 1},  /* Back */
    {-1, 0, 0}, /* Left */
    {0, 1, 0},  /* Top */
    {0, -1, 0}, /* Bottom */
};

/* Face base colors */
static const uint32_t cube_colors[6] = {
    0x5465ff, /* Front: blue */
    0x53917e, /* Right: teal */
    0xf19c79, /* Back: orange */
    0x6d1a36, /* Left: maroon */
    0xfcd0a1, /* Top: peach */
    0xb1b695, /* Bottom: olive */
};

/* Draw a single lit cube at position with rotation
 * @x, @y, @z: cube center position
 * @rot:       rotation angle in radians
 * @scale:     cube scale factor
 */
static void draw_cube(float x, float y, float z, float rot, float scale)
{
    b3d_reset();
    b3d_translate(x, y, z);
    b3d_rotate_y(rot);
    b3d_rotate_x(rot * 0.7f);
    b3d_scale(scale, scale, scale);

    for (int f = 0; f < 6; ++f) {
        const int ti = f * 2; /* Triangle index */
        const float *n = cube_normals[f];
        const uint32_t color = cube_colors[f];
        b3d_triangle_lit(&cube_faces[ti], n[0], n[1], n[2], color);
        b3d_triangle_lit(&cube_faces[ti + 1], n[0], n[1], n[2], color);
    }
}

/* Draw a thin quad as a line segment on the floor grid
 * @x0, @z0: start position
 * @x1, @z1: end position
 * @color:   line color
 */
static void draw_grid_line(float x0,
                           float z0,
                           float x1,
                           float z1,
                           uint32_t color)
{
    /* Calculate perpendicular direction for line width */
    float dx = x1 - x0;
    float dz = z1 - z0;
    float len = b3d_sqrtf(dx * dx + dz * dz);
    if (len < 0.001f)
        return;

    /* Perpendicular offset for line thickness */
    float px = -dz / len * LINE_WIDTH;
    float pz = dx / len * LINE_WIDTH;

    b3d_tri_t tri1 = {
        {{x0 - px, GRID_Y, z0 - pz},
         {x0 + px, GRID_Y, z0 + pz},
         {x1 + px, GRID_Y, z1 + pz}},
    };
    b3d_tri_t tri2 = {
        {{x0 - px, GRID_Y, z0 - pz},
         {x1 + px, GRID_Y, z1 + pz},
         {x1 - px, GRID_Y, z1 - pz}},
    };

    b3d_triangle(&tri1, color);
    b3d_triangle(&tri2, color);
}

/* Draw floor grid to demonstrate parallel line behavior
 * In perspective: lines converge to vanishing points
 * In orthographic: lines remain perfectly parallel
 */
static void draw_grid(void)
{
    b3d_reset();
    const uint32_t grid_color = 0x606060;
    float step = (2.0f * GRID_EXTENT) / (float) (GRID_LINES - 1);
    float z_far = (float) (CUBE_COUNT - 1) * CUBE_SPACING_Z + 2.0f;

    /* Lines along Z axis (demonstrate depth convergence) */
    for (int i = 0; i < GRID_LINES; ++i) {
        float x = -GRID_EXTENT + (float) i * step;
        draw_grid_line(x, -2.0f, x, z_far, grid_color);
    }

    /* Lines along X axis (cross lines for grid reference) */
    int z_lines = (int) ((z_far + 2.0f) / 2.0f) + 1;
    for (int i = 0; i < z_lines; ++i) {
        float z = -2.0f + (float) i * 2.0f;
        draw_grid_line(-GRID_EXTENT, z, GRID_EXTENT, z, grid_color);
    }
}

/* Render scene with cubes at varying distances
 * @pixels:    output pixel buffer
 * @depth:     depth buffer
 * @width:     buffer width
 * @height:    buffer height
 * @use_ortho: true for orthographic, false for perspective
 * @t:         time parameter for animation
 */
static void render_scene(uint32_t *pixels,
                         b3d_depth_t *depth,
                         int width,
                         int height,
                         bool use_ortho,
                         float t)
{
    const float aspect = (float) width / (float) height;

    b3d_init(pixels, depth, width, height, 60);
    b3d_clear();

    /* Set up lighting for better depth perception */
    b3d_set_light_direction(-0.5f, 0.8f, -0.3f);
    b3d_set_ambient(0.3f);

    /* Set up projection mode */
    if (use_ortho) {
        /* Orthographic projection:
         * - Fixed view volume, no perspective distortion
         * - Objects at all distances appear same size
         * - Parallel lines stay parallel (no vanishing point)
         */
        float ortho_size = 4.0f;
        b3d_ortho(-ortho_size * aspect, ortho_size * aspect, -ortho_size,
                  ortho_size, 0.1f, 100.0f);
    } else {
        /* Perspective projection:
         * - 60 degree field of view
         * - Natural depth perception (foreshortening)
         * - Parallel lines converge to vanishing points
         */
        b3d_set_fov(60.0f);
    }

    /* Camera position: looking toward +Z from negative Z */
    b3d_set_camera(&(b3d_camera_t) {0, 1.0f, -6, 0, 0, 0});

    /* Draw floor grid first (behind cubes) */
    draw_grid();

    /* Draw cubes at increasing Z distances
     * Key observation:
     *   Perspective:   distant cubes appear smaller
     *   Orthographic:  all cubes appear identical size
     */
    for (int i = 0; i < CUBE_COUNT; ++i) {
        float z = (float) i * CUBE_SPACING_Z;
        float x = ((float) i - (CUBE_COUNT - 1) / 2.0f) * CUBE_SPACING_X;
        draw_cube(x, 0, z, t + (float) i * 0.5f, CUBE_SIZE);
    }
}

int main(int argc, char **argv)
{
    int width = 800, height = 600;
    const char *snapshot = get_snapshot_path(argc, argv);
    bool use_ortho = false;

    const size_t pixel_count = (size_t) width * height;
    uint32_t *pixel_buffer = malloc(pixel_count * sizeof(uint32_t));
    b3d_depth_t *depth_buffer = malloc(pixel_count * sizeof(b3d_depth_t));

    if (!pixel_buffer || !depth_buffer) {
        free(pixel_buffer);
        free(depth_buffer);
        return 1;
    }

    /* Headless snapshot mode: render orthographic view */
    if (snapshot) {
        render_scene(pixel_buffer, depth_buffer, width, height, true, 0.5f);
        write_png(snapshot, pixel_buffer, width, height);
        free(pixel_buffer);
        free(depth_buffer);
        return 0;
    }

    /* Interactive SDL2 mode */
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        free(pixel_buffer);
        free(depth_buffer);
        return 1;
    }

    SDL_Window *window =
        SDL_CreateWindow("Perspective", SDL_WINDOWPOS_CENTERED,
                         SDL_WINDOWPOS_CENTERED, width, height, 0);
    if (!window) {
        fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
        free(pixel_buffer);
        free(depth_buffer);
        SDL_Quit();
        return 1;
    }

    SDL_Renderer *renderer = SDL_CreateRenderer(
        window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!renderer) {
        fprintf(stderr, "SDL_CreateRenderer failed: %s\n", SDL_GetError());
        free(pixel_buffer);
        free(depth_buffer);
        cleanup_sdl(NULL, NULL, window);
        return 1;
    }

    SDL_Texture *texture =
        SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888,
                          SDL_TEXTUREACCESS_STREAMING, width, height);
    if (!texture) {
        fprintf(stderr, "SDL_CreateTexture failed: %s\n", SDL_GetError());
        free(pixel_buffer);
        free(depth_buffer);
        cleanup_sdl(NULL, renderer, window);
        return 1;
    }

    printf("Orthographic vs Perspective Demo\n");
    printf("  O   - Toggle projection mode\n");
    printf("  ESC - Quit\n");
    printf(
        "Look for: grid lines converging in perspective and constant cube size "
        "in orthographic.\n\n");
    print_mode_hint(use_ortho);
    update_window_title(window, use_ortho);

    int quit = 0;
    while (!quit) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT ||
                (event.type == SDL_KEYDOWN &&
                 event.key.keysym.scancode == SDL_SCANCODE_ESCAPE)) {
                quit = 1;
                break;
            }
            /* Toggle orthographic mode with 'O' key */
            if (event.type == SDL_KEYDOWN &&
                event.key.keysym.scancode == SDL_SCANCODE_O) {
                use_ortho = !use_ortho;
                print_mode_hint(use_ortho);
                update_window_title(window, use_ortho);
            }
        }
        if (quit)
            break;

        float t = (float) SDL_GetTicks() * 0.001f;
        render_scene(pixel_buffer, depth_buffer, width, height, use_ortho, t);

        SDL_RenderClear(renderer);
        SDL_UpdateTexture(texture, NULL, pixel_buffer,
                          width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }

    free(pixel_buffer);
    free(depth_buffer);
    cleanup_sdl(texture, renderer, window);
    return 0;
}
