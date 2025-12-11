/*
    Animated heightmap inspired by the terrain and benchmark scenes in the
    external software renderers. Builds a simple sine/cosine height field and
    lets bootleg3d churn through thousands of small triangles.
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "SDL.h"
#include "b3d.h"
#include "pngwrite.h"

#define GRID_SIZE 64
#define CELL_SIZE 0.5f

static const char *get_snapshot_path(int argc, char **argv)
{
    const char *env = getenv("B3D_SNAPSHOT");
    if (env && env[0]) return env;
    for (int i = 1; i < argc; ++i) {
        const char *flag = "--snapshot=";
        size_t len = strlen(flag);
        if (strncmp(argv[i], flag, len) == 0) {
            return argv[i] + (int)len;
        }
    }
    return NULL;
}

static void write_png(const char *path, const uint32_t *rgba, int width, int height)
{
    FILE *file = fopen(path, "wb");
    if (!file) return;
    uint8_t *out = malloc((size_t)width * (size_t)height * 4);
    if (!out) {
        fclose(file);
        return;
    }
    size_t idx = 0;
    for (int i = 0; i < width * height; ++i) {
        uint32_t p = rgba[i];
        out[idx++] = (uint8_t)((p >> 16) & 0xff);
        out[idx++] = (uint8_t)((p >> 8) & 0xff);
        out[idx++] = (uint8_t)(p & 0xff);
        out[idx++] = 0xff;
    }
    png_write(file, (unsigned)width, (unsigned)height, out, true);
    free(out);
    fclose(file);
}

static float height_at(int x, int z, float t)
{
    float fx = (float)x * 0.3f;
    float fz = (float)z * 0.25f;
    return sinf(fx * 0.6f + t * 0.7f) * 0.6f + cosf(fz * 0.5f + t * 1.1f) * 0.4f;
}

static uint32_t height_color(float h)
{
    float shade = (h + 1.2f) * 0.4f;
    if (shade < 0.0f) shade = 0.0f;
    if (shade > 1.0f) shade = 1.0f;
    uint32_t r = (uint32_t)(80 + shade * 100.0f);
    uint32_t g = (uint32_t)(140 + shade * 110.0f);
    uint32_t b = (uint32_t)(90 + shade * 80.0f);
    return (r << 16) | (g << 8) | b;
}

static void render_heightmap(uint32_t *pixels, b3d_depth_t *depth, int width, int height, float t)
{
    b3d_init(pixels, depth, width, height, 70);
    b3d_set_camera(0.0f, 1.5f, -8.0f, 0.0f, 0.0f, 0.0f);

    float half_grid = (GRID_SIZE - 1) * CELL_SIZE * 0.5f;
    b3d_clear();

    /* Tilt the whole patch and slowly orbit it. */
    b3d_reset();
    b3d_rotate_y(t * 0.15f);
    b3d_rotate_x(-0.55f);
    b3d_translate(0.0f, -1.4f, 12.0f);

    for (int z = 0; z < GRID_SIZE - 1; ++z) {
        for (int x = 0; x < GRID_SIZE - 1; ++x) {
            float h00 = height_at(x,     z,     t);
            float h10 = height_at(x + 1, z,     t);
            float h01 = height_at(x,     z + 1, t);
            float h11 = height_at(x + 1, z + 1, t);

            float fx = (x * CELL_SIZE) - half_grid;
            float fz = (z * CELL_SIZE);
            float fx_next = fx + CELL_SIZE;
            float fz_next = fz + CELL_SIZE;

            uint32_t c0 = height_color((h00 + h10 + h11) / 3.0f);
            uint32_t c1 = height_color((h00 + h11 + h01) / 3.0f);

            /* Winding flipped so culling keeps the patch visible when tilted. */
            b3d_triangle(fx,      h00, fz,      fx_next, h11, fz_next, fx_next, h10, fz,      c0);
            b3d_triangle(fx,      h00, fz,      fx,      h01, fz_next, fx_next, h11, fz_next, c1);
        }
    }
}

int main(int argc, char **argv)
{
    int width = 960;
    int height = 720;
    const char *snapshot = get_snapshot_path(argc, argv);

    uint32_t *pixels = malloc(width * height * sizeof(pixels[0]));
    b3d_depth_t *depth = malloc(width * height * sizeof(depth[0]));

    /* Headless snapshot mode for CI/docs. */
    if (snapshot) {
        render_heightmap(pixels, depth, width, height, 1.0f);
        write_png(snapshot, pixels, width, height);
        free(pixels);
        free(depth);
        return 0;
    }

    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window *window = SDL_CreateWindow("Heightmap", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, 0);
    SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC);
    SDL_Texture *texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, width, height);

    SDL_SetWindowTitle(window, "Heightmap (sine/cosine)");

    int quit = 0;
    while (!quit) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT ||
               (event.type == SDL_KEYDOWN && event.key.keysym.scancode == SDL_SCANCODE_ESCAPE)) {
                quit = 1;
                break;
            }
        }
        if (quit) break;

        float t = SDL_GetTicks() * 0.001f;
        render_heightmap(pixels, depth, width, height, t);

        SDL_RenderClear(renderer);
        SDL_UpdateTexture(texture, NULL, pixels, width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }

    free(pixels);
    free(depth);
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
