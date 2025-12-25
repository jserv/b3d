/*
 * Lighting demo - demonstrates B3D's basic lighting support.
 *
 * Controls:
 *   Arrow keys: rotate light direction
 *   +/-: adjust ambient level
 *   Space: toggle auto-rotate light
 *   1/2: switch between cube and pyramid
 *   ESC: quit
 *
 * Supports headless snapshots with '--snapshot=PATH' or B3D_SNAPSHOT.
 */

#include <SDL.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "b3d-math.h"
#include "b3d.h"
#include "utils.h"

#define WIDTH 640
#define HEIGHT 480

/* Cube geometry with face normals */
static const struct {
    b3d_tri_t tri;
    float nx, ny, nz;
} cube_faces[] = {
    /* Front face (+Z) */
    {{{
         {-0.5f, -0.5f, 0.5f},
         {0.5f, -0.5f, 0.5f},
         {0.5f, 0.5f, 0.5f},
     }},
     0.0f,
     0.0f,
     1.0f},
    {{{
         {-0.5f, -0.5f, 0.5f},
         {0.5f, 0.5f, 0.5f},
         {-0.5f, 0.5f, 0.5f},
     }},
     0.0f,
     0.0f,
     1.0f},
    /* Back face (-Z) */
    {{{
         {0.5f, -0.5f, -0.5f},
         {-0.5f, -0.5f, -0.5f},
         {-0.5f, 0.5f, -0.5f},
     }},
     0.0f,
     0.0f,
     -1.0f},
    {{{
         {0.5f, -0.5f, -0.5f},
         {-0.5f, 0.5f, -0.5f},
         {0.5f, 0.5f, -0.5f},
     }},
     0.0f,
     0.0f,
     -1.0f},
    /* Right face (+X) */
    {{{
         {0.5f, -0.5f, 0.5f},
         {0.5f, -0.5f, -0.5f},
         {0.5f, 0.5f, -0.5f},
     }},
     1.0f,
     0.0f,
     0.0f},
    {{{
         {0.5f, -0.5f, 0.5f},
         {0.5f, 0.5f, -0.5f},
         {0.5f, 0.5f, 0.5f},
     }},
     1.0f,
     0.0f,
     0.0f},
    /* Left face (-X) */
    {{{
         {-0.5f, -0.5f, -0.5f},
         {-0.5f, -0.5f, 0.5f},
         {-0.5f, 0.5f, 0.5f},
     }},
     -1.0f,
     0.0f,
     0.0f},
    {{{
         {-0.5f, -0.5f, -0.5f},
         {-0.5f, 0.5f, 0.5f},
         {-0.5f, 0.5f, -0.5f},
     }},
     -1.0f,
     0.0f,
     0.0f},
    /* Top face (+Y) */
    {{{
         {-0.5f, 0.5f, 0.5f},
         {0.5f, 0.5f, 0.5f},
         {0.5f, 0.5f, -0.5f},
     }},
     0.0f,
     1.0f,
     0.0f},
    {{{
         {-0.5f, 0.5f, 0.5f},
         {0.5f, 0.5f, -0.5f},
         {-0.5f, 0.5f, -0.5f},
     }},
     0.0f,
     1.0f,
     0.0f},
    /* Bottom face (-Y) */
    {{{
         {-0.5f, -0.5f, -0.5f},
         {0.5f, -0.5f, -0.5f},
         {0.5f, -0.5f, 0.5f},
     }},
     0.0f,
     -1.0f,
     0.0f},
    {{{
         {-0.5f, -0.5f, -0.5f},
         {0.5f, -0.5f, 0.5f},
         {-0.5f, -0.5f, 0.5f},
     }},
     0.0f,
     -1.0f,
     0.0f},
};
#define CUBE_FACE_COUNT (sizeof(cube_faces) / sizeof(cube_faces[0]))

/* Pyramid geometry with face normals */
static const struct {
    b3d_tri_t tri;
    float nx, ny, nz;
} pyramid_faces[] = {
    /* Base (two triangles) */
    {{{
         {-0.5f, -0.5f, -0.5f},
         {0.5f, -0.5f, -0.5f},
         {0.5f, -0.5f, 0.5f},
     }},
     0.0f,
     -1.0f,
     0.0f},
    {{{
         {-0.5f, -0.5f, -0.5f},
         {0.5f, -0.5f, 0.5f},
         {-0.5f, -0.5f, 0.5f},
     }},
     0.0f,
     -1.0f,
     0.0f},
    /* Front face */
    {{{
         {-0.5f, -0.5f, 0.5f},
         {0.5f, -0.5f, 0.5f},
         {0.0f, 0.5f, 0.0f},
     }},
     0.0f,
     0.4472f,
     0.8944f},
    /* Right face */
    {{{
         {0.5f, -0.5f, 0.5f},
         {0.5f, -0.5f, -0.5f},
         {0.0f, 0.5f, 0.0f},
     }},
     0.8944f,
     0.4472f,
     0.0f},
    /* Back face */
    {{{
         {0.5f, -0.5f, -0.5f},
         {-0.5f, -0.5f, -0.5f},
         {0.0f, 0.5f, 0.0f},
     }},
     0.0f,
     0.4472f,
     -0.8944f},
    /* Left face */
    {{{
         {-0.5f, -0.5f, -0.5f},
         {-0.5f, -0.5f, 0.5f},
         {0.0f, 0.5f, 0.0f},
     }},
     -0.8944f,
     0.4472f,
     0.0f},
};
#define PYRAMID_FACE_COUNT (sizeof(pyramid_faces) / sizeof(pyramid_faces[0]))

/* Light direction angles (spherical coordinates) */
static float light_theta = 0.5f; /* Azimuth angle */
static float light_phi = 0.5f;   /* Elevation angle */
static bool auto_rotate = true;
static int shape = 0; /* 0 = cube, 1 = pyramid */

/* Update light direction from spherical coordinates */
static void update_light(void)
{
    float sin_phi, cos_phi, sin_theta, cos_theta;
    b3d_sincosf(light_phi, &sin_phi, &cos_phi);
    b3d_sincosf(light_theta, &sin_theta, &cos_theta);
    float x = cos_phi * sin_theta;
    float y = sin_phi;
    float z = cos_phi * cos_theta;
    b3d_set_light_direction(x, y, z);
}

/* Render the scene */
static void render(uint32_t *pixels, b3d_depth_t *depth, float t)
{
    b3d_init(pixels, depth, WIDTH, HEIGHT, 60);
    b3d_set_camera(&(b3d_camera_t) {0, 0, -3, 0, 0, 0});
    b3d_clear();

    /* Rotate the object */
    b3d_reset();
    b3d_rotate_y(t * 0.5f);
    b3d_rotate_x(t * 0.3f);

    /* Draw lit geometry */
    uint32_t color = 0x4488FF; /* Blue base color */

    if (shape == 0) {
        for (size_t i = 0; i < CUBE_FACE_COUNT; i++)
            b3d_triangle_lit(&cube_faces[i].tri, cube_faces[i].nx,
                             cube_faces[i].ny, cube_faces[i].nz, color);
    } else {
        color = 0xFF8844; /* Orange for pyramid */
        for (size_t i = 0; i < PYRAMID_FACE_COUNT; i++)
            b3d_triangle_lit(&pyramid_faces[i].tri, pyramid_faces[i].nx,
                             pyramid_faces[i].ny, pyramid_faces[i].nz, color);
    }
}

int main(int argc, char **argv)
{
    const char *snapshot = get_snapshot_path(argc, argv);

    /* Allocate buffers */
    uint32_t *pixels = malloc((size_t) WIDTH * HEIGHT * sizeof(uint32_t));
    b3d_depth_t *depth = malloc((size_t) WIDTH * HEIGHT * sizeof(b3d_depth_t));
    if (!pixels || !depth) {
        fprintf(stderr, "Failed to allocate buffers\n");
        return 1;
    }

    /* Headless snapshot mode */
    if (snapshot) {
        update_light();
        render(pixels, depth, 0.8f);
        write_png(snapshot, pixels, WIDTH, HEIGHT);
        printf("Snapshot saved to %s\n", snapshot);
        free(pixels);
        free(depth);
        return 0;
    }

    /* SDL initialization */
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL init failed: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window *window =
        SDL_CreateWindow("B3D Lighting Demo", SDL_WINDOWPOS_CENTERED,
                         SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, 0);
    if (!window) {
        fprintf(stderr, "Window creation failed: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, 0);
    SDL_Texture *texture =
        SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888,
                          SDL_TEXTUREACCESS_STREAMING, WIDTH, HEIGHT);

    bool running = true;
    Uint32 start_time = SDL_GetTicks();
    float ambient = 0.2f;

    printf("Lighting Demo Controls:\n");
    printf("  Arrows: rotate light   +/-: ambient   Space: auto-rotate\n");
    printf("  1/2: cube/pyramid      ESC: quit\n");

    while (running) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT)
                running = false;
            else if (event.type == SDL_KEYDOWN) {
                switch (event.key.keysym.sym) {
                case SDLK_ESCAPE:
                    running = false;
                    break;
                case SDLK_LEFT:
                    light_theta -= 0.1f;
                    auto_rotate = false;
                    break;
                case SDLK_RIGHT:
                    light_theta += 0.1f;
                    auto_rotate = false;
                    break;
                case SDLK_UP:
                    light_phi += 0.1f;
                    if (light_phi > 1.5f)
                        light_phi = 1.5f;
                    auto_rotate = false;
                    break;
                case SDLK_DOWN:
                    light_phi -= 0.1f;
                    if (light_phi < -1.5f)
                        light_phi = -1.5f;
                    auto_rotate = false;
                    break;
                case SDLK_EQUALS:
                case SDLK_PLUS:
                    ambient += 0.05f;
                    if (ambient > 1.0f)
                        ambient = 1.0f;
                    b3d_set_ambient(ambient);
                    printf("Ambient: %.0f%%\n", ambient * 100);
                    break;
                case SDLK_MINUS:
                    ambient -= 0.05f;
                    if (ambient < 0.0f)
                        ambient = 0.0f;
                    b3d_set_ambient(ambient);
                    printf("Ambient: %.0f%%\n", ambient * 100);
                    break;
                case SDLK_SPACE:
                    auto_rotate = !auto_rotate;
                    printf("Auto-rotate: %s\n", auto_rotate ? "ON" : "OFF");
                    break;
                case SDLK_1:
                    shape = 0;
                    printf("Shape: Cube\n");
                    break;
                case SDLK_2:
                    shape = 1;
                    printf("Shape: Pyramid\n");
                    break;
                }
            }
        }

        float t = (SDL_GetTicks() - start_time) / 1000.0f;

        /* Auto-rotate light */
        if (auto_rotate)
            light_theta = t * 0.8f;

        update_light();
        render(pixels, depth, t);

        SDL_UpdateTexture(texture, NULL, pixels, WIDTH * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);

        SDL_Delay(16);
    }

    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    free(pixels);
    free(depth);

    return 0;
}
