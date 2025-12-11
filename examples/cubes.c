/*
 * This program draws spinning cubes using SDL2.
 *
 * This can be used as a crude benchmark, as it will add more cubes to the scene
 * until the framerate is around 60fps. It supports headless snapshots with
 * '--snapshot=PATH' or B3D_SNAPSHOT.
 */

#include <SDL.h>
#include <stdlib.h>
#include <string.h>

#include "b3d.h"
#include "pngwrite.h"

static const char *get_snapshot_path(int argc, char **argv)
{
    const char *env = getenv("B3D_SNAPSHOT");
    if (env && env[0])
        return env;

    for (int i = 1; i < argc; ++i) {
        const char *flag = "--snapshot=";
        size_t len = strlen(flag);
        if (!strncmp(argv[i], flag, len))
            return argv[i] + (int) len;
    }
    return NULL;
}

static void write_png(const char *path,
                      const uint32_t *rgba,
                      int width,
                      int height)
{
    FILE *file = fopen(path, "wb");
    if (!file)
        return;
    uint8_t *out = malloc((size_t) width * (size_t) height * 4);
    if (!out) {
        fclose(file);
        return;
    }
    size_t idx = 0;
    for (int i = 0; i < width * height; ++i) {
        uint32_t p = rgba[i];
        out[idx++] = (uint8_t) ((p >> 16) & 0xff);
        out[idx++] = (uint8_t) ((p >> 8) & 0xff);
        out[idx++] = (uint8_t) (p & 0xff);
        out[idx++] = 0xff;
    }
    png_write(file, (unsigned) width, (unsigned) height, out, true);
    free(out);
    fclose(file);
}

static void render_cubes(uint32_t *pixels,
                         b3d_depth_t *depth,
                         int width,
                         int height,
                         int cube_count,
                         float t)
{
    b3d_init(pixels, depth, width, height, 60);
    b3d_set_camera(0, 0, -2, 0, 0, 0);
    b3d_clear();

    for (int i = 0; i < cube_count; ++i) {
        b3d_reset();
        b3d_rotate_z(t);
        b3d_rotate_y(t);
        b3d_rotate_x(t);
        b3d_rotate_y(i * 0.1);
        b3d_translate(1, 1, fmodf(i * 0.1, 100));
        b3d_rotate_z(i + t);

        b3d_triangle(-0.5, -0.5, -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5,
                     0xfcd0a1);
        b3d_triangle(-0.5, -0.5, -0.5, 0.5, 0.5, -0.5, 0.5, -0.5, -0.5,
                     0xb1b695);
        b3d_triangle(0.5, -0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5, 0x53917e);
        b3d_triangle(0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0.5, 0x63535b);
        b3d_triangle(0.5, -0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0x6d1a36);
        b3d_triangle(0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0xd4e09b);
        b3d_triangle(-0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, -0.5,
                     0xf6f4d2);
        b3d_triangle(-0.5, -0.5, 0.5, -0.5, 0.5, -0.5, -0.5, -0.5, -0.5,
                     0xcbdfbd);
        b3d_triangle(-0.5, 0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0xf19c79);
        b3d_triangle(-0.5, 0.5, -0.5, 0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0xa44a3f);
        b3d_triangle(0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, -0.5,
                     0x5465ff);
        b3d_triangle(0.5, -0.5, 0.5, -0.5, -0.5, -0.5, 0.5, -0.5, -0.5,
                     0x788bff);
    }
}

int main(int argc, char **argv)
{
    int width = 800, height = 600;
    const char *snapshot = get_snapshot_path(argc, argv);

    uint32_t *pixel_buffer = malloc(width * height * sizeof(uint32_t));
    b3d_depth_t *depth_buffer = malloc(width * height * sizeof(b3d_depth_t));

    if (snapshot) {
        render_cubes(pixel_buffer, depth_buffer, width, height, 100, 1.2f);
        write_png(snapshot, pixel_buffer, width, height);
        free(pixel_buffer);
        free(depth_buffer);
        return 0;
    }

    // Set up SDL2 stuff, including whats needed to display a buffer of pixels.
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window *window = SDL_CreateWindow(
        "", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, 0);
    SDL_Renderer *renderer =
        SDL_CreateRenderer(window, -1, 0 /*SDL_RENDERER_PRESENTVSYNC*/);
    SDL_Texture *texture =
        SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888,
                          SDL_TEXTUREACCESS_STREAMING, width, height);

    b3d_init(pixel_buffer, depth_buffer, width, height, 60);

    // For framerate counting.
    double freq = SDL_GetPerformanceFrequency();
    uint32_t next_update = 0;
#define FPS_SAMPLES 100
    float average_fps[FPS_SAMPLES];
    int average_index = 0;
    int have_enough_samples = 0;

    int cube_count = 100;

    b3d_set_camera(0, 0, -2, 0, 0, 0);

    int quit = 0;
    while (!quit) {
        uint64_t time_stamp = SDL_GetPerformanceCounter();

        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT ||
                (event.type == SDL_KEYDOWN &&
                 event.key.keysym.scancode == SDL_SCANCODE_ESCAPE)) {
                quit = 1;
                break;
            }
        }
        if (quit)
            break;

        // Current time in milliseconds.
        float t = SDL_GetTicks() * 0.001f;

        b3d_clear();

        for (int i = 0; i < cube_count; ++i) {
            // Reset transformations back to the origin.
            b3d_reset();

            b3d_rotate_z(t);
            b3d_rotate_y(t);
            b3d_rotate_x(t);

            b3d_rotate_y(i * 0.1);
            b3d_translate(1, 1, fmodf(i * 0.1, 100));
            b3d_rotate_z(i + t);

            b3d_triangle(-0.5, -0.5, -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5,
                         0xfcd0a1);
            b3d_triangle(-0.5, -0.5, -0.5, 0.5, 0.5, -0.5, 0.5, -0.5, -0.5,
                         0xb1b695);
            b3d_triangle(0.5, -0.5, -0.5, 0.5, 0.5, -0.5, 0.5, 0.5, 0.5,
                         0x53917e);
            b3d_triangle(0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0.5,
                         0x63535b);
            b3d_triangle(0.5, -0.5, 0.5, 0.5, 0.5, 0.5, -0.5, 0.5, 0.5,
                         0x6d1a36);
            b3d_triangle(0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5,
                         0xd4e09b);
            b3d_triangle(-0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5, -0.5,
                         0xf6f4d2);
            b3d_triangle(-0.5, -0.5, 0.5, -0.5, 0.5, -0.5, -0.5, -0.5, -0.5,
                         0xcbdfbd);
            b3d_triangle(-0.5, 0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                         0xf19c79);
            b3d_triangle(-0.5, 0.5, -0.5, 0.5, 0.5, 0.5, 0.5, 0.5, -0.5,
                         0xa44a3f);
            b3d_triangle(0.5, -0.5, 0.5, -0.5, -0.5, 0.5, -0.5, -0.5, -0.5,
                         0x5465ff);
            b3d_triangle(0.5, -0.5, 0.5, -0.5, -0.5, -0.5, 0.5, -0.5, -0.5,
                         0x788bff);
        }

        // Display the pixel buffer on screen (using a streaming texture).
        SDL_RenderClear(renderer);
        SDL_UpdateTexture(texture, NULL, pixel_buffer,
                          width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);

        // Display the average framerate in the window title.
        if (SDL_GetTicks() > next_update && have_enough_samples) {
            char title[32];
            float fps = 0;
            for (int i = 0; i < FPS_SAMPLES; ++i)
                fps += average_fps[i];
            fps /= FPS_SAMPLES;
            cube_count += (fps > 60) ? 50 : -50;
            snprintf(title, 32, "%d tris, %.1f fps", cube_count * 12, fps);
            SDL_SetWindowTitle(window, title);
            next_update = SDL_GetTicks() + 250;
        }
        average_fps[average_index++] =
            1.0 / ((SDL_GetPerformanceCounter() - time_stamp) / freq);
        if (average_index >= FPS_SAMPLES) {
            average_index = 0;
            have_enough_samples = 1;
        }
    }

    free(pixel_buffer);
    free(depth_buffer);
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
