/*
    Torus (“donut”) demo inspired by externals/donut.c. Generates a parametric
    torus, applies simple directional lighting, and renders via the public
    b3d API. Supports headless snapshots with --snapshot=PATH or
    B3D_SNAPSHOT.
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "SDL.h"
#include "b3d.h"
#include "pngwrite.h"

/* Warm-to-cool gradient for nicer lighting */
static const uint32_t palette[] = {
    0x0f1028, 0x14163b, 0x1a1d4e, 0x1f245f, 0x232b70, 0x29327f, 0x2f3990,
    0x3541a1, 0x3c49b3, 0x4352c4, 0x4b5bd4, 0x5365e3, 0x5d6eec, 0x6677f3,
    0x7281f8, 0x7f8bfb, 0x8c94fa, 0x9b9ff5, 0xaaa8ec, 0xbab2e0, 0xcbbbd0,
    0xdcc5bc, 0xeecfa5, 0xf7d88e, 0xfde07a, 0xfdd567, 0xfbc556, 0xf7b445,
    0xf3a235, 0xee9028, 0xe77d1b, 0xdd6911, 0xd05509,
};

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

static uint32_t shade_color(float dot)
{
    if (dot < 0.0f) dot = 0.0f;
    if (dot > 1.0f) dot = 1.0f;
    int idx = (int)(dot * (sizeof(palette) / sizeof(palette[0]) - 1));
    return palette[idx];
}

static void render_frame(uint32_t *pixels, b3d_depth_t *depth, int width, int height, float t)
{
    b3d_init(pixels, depth, width, height, 70.0f);
    b3d_set_camera(0.0f, 0.0f, -6.0f, 0.0f, 0.0f, 0.0f);
    b3d_clear();

    b3d_reset();
    b3d_rotate_y(t * 0.6f);
    b3d_rotate_x(t * 0.35f);

    const float R = 2.0f;   /* major radius */
    const float r = 0.7f;   /* minor radius */
    const int segU = 96;    /* higher tessellation for smoother surface */
    const int segV = 64;
    const float du = (2.0f * 3.14159265f) / segU;
    const float dv = (2.0f * 3.14159265f) / segV;
    float lx = 0.3f, ly = 0.8f, lz = -0.6f;
    float llen = sqrtf(lx * lx + ly * ly + lz * lz);
    if (llen > 0.0f) {
        lx /= llen; ly /= llen; lz /= llen;
    }

    for (int iu = 0; iu < segU; ++iu) {
        float u0 = iu * du;
        float u1 = (iu + 1) * du;
        float cu0 = cosf(u0), su0 = sinf(u0);
        float cu1 = cosf(u1), su1 = sinf(u1);
        for (int iv = 0; iv < segV; ++iv) {
            float v0 = iv * dv;
            float v1 = (iv + 1) * dv;
            float cv0 = cosf(v0), sv0 = sinf(v0);
            float cv1 = cosf(v1), sv1 = sinf(v1);

            /* Four corners of the quad */
            float x00 = (R + r * cv0) * cu0;
            float y00 = (R + r * cv0) * su0;
            float z00 = r * sv0;

            float x10 = (R + r * cv0) * cu1;
            float y10 = (R + r * cv0) * su1;
            float z10 = r * sv0;

            float x01 = (R + r * cv1) * cu0;
            float y01 = (R + r * cv1) * su0;
            float z01 = r * sv1;

            float x11 = (R + r * cv1) * cu1;
            float y11 = (R + r * cv1) * su1;
            float z11 = r * sv1;

            /* Per-vertex normals of parametric torus */
            float nx00 = cu0 * cv0;
            float ny00 = su0 * cv0;
            float nz00 = sv0;

            float nx10 = cu1 * cv0;
            float ny10 = su1 * cv0;
            float nz10 = sv0;

            float nx01 = cu0 * cv1;
            float ny01 = su0 * cv1;
            float nz01 = sv1;

            float nx11 = cu1 * cv1;
            float ny11 = su1 * cv1;
            float nz11 = sv1;

            /* Two triangles per quad; average vertex normals for shading */
            float dot0 = (nx00 * lx + ny00 * ly + nz00 * lz +
                          nx10 * lx + ny10 * ly + nz10 * lz +
                          nx11 * lx + ny11 * ly + nz11 * lz) / 3.0f;
            float dot1 = (nx00 * lx + ny00 * ly + nz00 * lz +
                          nx11 * lx + ny11 * ly + nz11 * lz +
                          nx01 * lx + ny01 * ly + nz01 * lz) / 3.0f;
            /* Gentle gamma to reduce banding across palette steps */
            dot0 = powf(fmaxf(dot0, 0.0f), 0.8f);
            dot1 = powf(fmaxf(dot1, 0.0f), 0.8f);

            uint32_t c0 = shade_color(dot0);
            uint32_t c1 = shade_color(dot1);

            b3d_triangle(x00, y00, z00, x10, y10, z10, x11, y11, z11, c0);
            b3d_triangle(x00, y00, z00, x11, y11, z11, x01, y01, z01, c1);
        }
    }
}

int main(int argc, char **argv)
{
    int width = 960;
    int height = 720;
    const char *snapshot = get_snapshot_path(argc, argv);

    uint32_t *pixels = malloc((size_t)width * (size_t)height * sizeof(pixels[0]));
    b3d_depth_t *depth = malloc((size_t)width * (size_t)height * sizeof(depth[0]));

    if (snapshot) {
        render_frame(pixels, depth, width, height, 1.4f);
        write_png(snapshot, pixels, width, height);
        free(pixels);
        free(depth);
        return 0;
    }

    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window *window = SDL_CreateWindow("Donut", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, 0);
    SDL_Renderer *renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC);
    SDL_Texture *texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, width, height);
    SDL_SetWindowTitle(window, "Donut (torus) demo");

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
        render_frame(pixels, depth, width, height, t);

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
