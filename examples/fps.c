/*
 * A weird little demo game where you collect golden heads. It demonstrates how
 * immediate-mode rendering of triangles can be used to throw together things
 * with simple code.
 *
 * The head model is loaded from an OBJ file, and the rest of the terrain is
 * generated randomly. This uses a trick where the random number generator seed
 * is reset before use to recall the same set of numbers without having to store
 * them. It supports headless snapshots with --snapshot=PATH or B3D_SNAPSHOT.
 *
 * Use WASD or the arrow keys to move, and mouse to look.
 */

#include <SDL.h>
#include <string.h>
#include <time.h>

#include "b3d.h"
#include "b3d_obj.h"
#include "pngwrite.h"

#define RND (rand() / (float) RAND_MAX)

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

int main(int argument_count, char **arguments)
{
    int width = 640;
    int height = 480;
    const char *snapshot = get_snapshot_path(argument_count, arguments);

    setbuf(stdout, NULL);

    uint32_t *pixels = malloc(width * height * sizeof(pixels[0]));
    b3d_depth_t *depth = malloc(width * height * sizeof(depth[0]));
    b3d_init(pixels, depth, width, height, 90);

    /* Load OBJ file */
    const char *file_name = "assets/moai.obj";
    for (int i = 1; i < argument_count; ++i) {
        if (strncmp(arguments[i], "--snapshot=", 11) != 0) {
            file_name = arguments[i];
            break;
        }
    }

    b3d_mesh_t mesh;
    int err = b3d_load_obj(file_name, &mesh);
    if (err == 1) {
        printf("Failed to load file '%s'.\n", file_name);
        exit(1);
    } else if (err == 2) {
        printf("Memory allocation failed.\n");
        exit(1);
    } else if (err == 3) {
        printf("Invalid vertex index in OBJ file.\n");
        exit(1);
    }

    float world_size = 20;
    float boundary = 2.5f;
    int seed = 12345; /* Fixed seed for reproducible snapshot */

    if (snapshot) {
        srand(seed);
        float heads[16];
        for (int i = 0; i < 16; ++i)
            heads[i] = (int) ((RND * (world_size - boundary) * 2) -
                              (world_size - boundary)) |
                       1;

        b3d_clear();
        b3d_set_camera(-1.0f, 1.0f, -3.0f, 1.0f, 0.0f, 0);

        /* Draw checkerboard floor */
        b3d_reset();
        for (int z = -world_size; z < world_size; ++z) {
            for (int x = -world_size; x < world_size; ++x) {
                uint32_t c = (x + z) & 1 ? 0x424C88 : 0xF7C396;
                b3d_triangle(x + .5, 0.0, z + .5, x - .5, 0.0, z - .5, x - .5,
                             0.0, z + .5, c);
                b3d_triangle(x + .5, 0.0, z + .5, x + .5, 0.0, z - .5, x - .5,
                             0.0, z - .5, c);
            }
        }

        /* Draw golden heads */
        float t = 1.5f;
        for (int h = 0; h < 16; h += 2) {
            float hx = heads[h], hz = heads[h + 1];
            b3d_reset();
            b3d_rotate_y(h + t * 3);
            b3d_scale(0.4, 0.4, 0.4);
            b3d_translate(hx, 0.4 + sinf(h + t * 3) * 0.1, hz);
            srand(h);
            for (int i = 0; i < mesh.vertex_count; i += 9) {
                uint32_t r = 200 + (int) (RND * 50);
                uint32_t g = 150 + (int) (RND * 50);
                uint32_t b = 50 + (int) (RND * 50);
                b3d_triangle(mesh.triangles[i + 0], mesh.triangles[i + 1],
                             mesh.triangles[i + 2], mesh.triangles[i + 3],
                             mesh.triangles[i + 4], mesh.triangles[i + 5],
                             mesh.triangles[i + 6], mesh.triangles[i + 7],
                             mesh.triangles[i + 8], (r << 16 | g << 8 | b));
            }
        }

        /* Draw border cubes */
        srand(seed);
        for (int i = -world_size; i < world_size; i += 2) {
            for (int j = 0; j < 4; ++j) {
                float bx = i, bz = i;
                switch (j) {
                case 0:
                    bx = -world_size;
                    break;
                case 1:
                    bx = world_size;
                    break;
                case 2:
                    bz = -world_size;
                    break;
                case 3:
                    bz = world_size;
                    break;
                }
                b3d_reset();
                b3d_rotate_y(RND * 3.14159f);
                b3d_rotate_x(RND * 3.14159f);
                b3d_rotate_z(RND * 3.14159f);
                b3d_scale(1.0f + RND * 2, 1.0f + RND * 8, 1.0f + RND * 2);
                b3d_translate(bx, 0.5, bz);
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
        }

        /* Draw pyramids */
        srand(seed);
        for (int i = 0; i < 20; ++i) {
            b3d_reset();
            b3d_scale(1, 1 + RND * 3, 1);
            b3d_rotate_y(RND * 3.14159f);
            b3d_translate((RND * world_size * 2) - world_size, 0,
                          (RND * world_size * 2) - world_size);
            b3d_triangle(0.0, 2.0, 0.0, -1.0, 0.0, 1.0, 1.0, 0.0, 1.0,
                         0x004749);
            b3d_triangle(0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, -1.0,
                         0x00535a);
            b3d_triangle(0.0, 2.0, 0.0, 1.0, 0.0, -1.0, -1.0, 0.0, -1.0,
                         0x00746b);
            b3d_triangle(0.0, 2.0, 0.0, -1.0, 0.0, -1.0, -1.0, 0.0, 1.0,
                         0x00945c);
        }

        write_png(snapshot, pixels, width, height);
        free(pixels);
        free(depth);
        b3d_free_mesh(&mesh);
        return 0;
    }

    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window *window = SDL_CreateWindow(
        "", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, 0);
    SDL_Renderer *renderer =
        SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC);
    SDL_Texture *texture =
        SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888,
                          SDL_TEXTUREACCESS_STREAMING, width, height);

    float player_x = -1.0f;
    float player_z = -3.0f;
    float player_height = 1.0f;
    float player_yaw = 1.0f;
    float player_pitch = 0.0f;
    float player_forward_speed = 0.0f;
    float player_strafe_speed = 0.0f;
    float mouse_sensitivity = 0.001f;
    int up = 0, down = 0, left = 0, right = 0, crouch = 0;
    SDL_SetRelativeMouseMode(1);

    seed = time(NULL);
    srand(seed);
    float heads[16];
    for (int i = 0; i < 16; ++i) {
        heads[i] = (int) ((RND * (world_size - boundary) * 2) -
                          (world_size - boundary)) |
                   1;
    }
    int heads_found = 0;
    float head_radius = 0.5f;
    float head_confetti_y = 2.0f;

    SDL_SetWindowTitle(window, "Find The Golden Heads");

    int quit = 0;
    while (!quit) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
                break;
            }

            if (event.type == SDL_KEYDOWN || event.type == SDL_KEYUP) {
                int sc = event.key.keysym.scancode;
                if (sc == SDL_SCANCODE_UP || sc == SDL_SCANCODE_W) {
                    up = event.key.state;
                } else if (sc == SDL_SCANCODE_DOWN || sc == SDL_SCANCODE_S) {
                    down = event.key.state;
                } else if (sc == SDL_SCANCODE_LEFT || sc == SDL_SCANCODE_A) {
                    left = event.key.state;
                } else if (sc == SDL_SCANCODE_RIGHT || sc == SDL_SCANCODE_D) {
                    right = event.key.state;
                } else if (sc == SDL_SCANCODE_LSHIFT ||
                           sc == SDL_SCANCODE_RSHIFT || sc == SDL_SCANCODE_C) {
                    crouch = event.key.state;
                } else if (sc == SDL_SCANCODE_ESCAPE) {
                    quit = 1;
                    break;
                }
            } else if (event.type == SDL_MOUSEMOTION) {
                player_yaw -= event.motion.xrel * mouse_sensitivity;
                player_pitch += event.motion.yrel * mouse_sensitivity;
                if (player_pitch < -1.570796f)
                    player_pitch = -1.57f;
                if (player_pitch > 1.570796f)
                    player_pitch = 1.57f;
            }
        }
        if (quit)
            break;

        b3d_clear();
        float t = SDL_GetTicks() * 0.001f;

        /* A basic version of standard (mouse-look) FPS controls */
        if (up)
            player_forward_speed = 0.1f;
        if (down)
            player_forward_speed = -0.1f;
        if (left)
            player_strafe_speed = 0.1f;
        if (right)
            player_strafe_speed = -0.1f;
        player_height += ((crouch ? 0.5f : 1.0f) - player_height) * 0.1f;
        player_x -= cosf(player_yaw - 1.570796f) * player_forward_speed;
        player_z -= sinf(player_yaw - 1.570796f) * player_forward_speed;
        player_x -= cosf(player_yaw) * player_strafe_speed;
        player_z -= sinf(player_yaw) * player_strafe_speed;
        player_forward_speed *= 0.9f;
        player_strafe_speed *= 0.9f;

        /* Collision with world edge */
        if (player_x < -(world_size - boundary))
            player_x = -world_size + boundary;
        if (player_x > (world_size - boundary))
            player_x = world_size - boundary;
        if (player_z < -(world_size - boundary))
            player_z = -world_size + boundary;
        if (player_z > (world_size - boundary))
            player_z = world_size - boundary;

        b3d_set_camera(player_x, player_height, player_z, player_yaw,
                       player_pitch, 0);

        /* Draw a checkerboard floor */
        b3d_reset();
        for (int z = -world_size; z < world_size; ++z) {
            for (int x = -world_size; x < world_size; ++x) {
                uint32_t c = (x + z) & 1 ? 0x424C88 : 0xF7C396;
                b3d_triangle(x + .5, 0.0, z + .5, x - .5, 0.0, z - .5, x - .5,
                             0.0, z + .5, c);
                b3d_triangle(x + .5, 0.0, z + .5, x + .5, 0.0, z - .5, x - .5,
                             0.0, z - .5, c);
            }
        }

        /* Draw the golden heads */
        for (int h = 0; h < 16; h += 2) {
            float x = heads[h];
            float z = heads[h + 1];
            /* If the player is close to a head, mark it as found using NAN */
            if (fabsf(player_x - x) < head_radius &&
                fabsf(player_z - z) < head_radius) {
                heads[h] = NAN;
                heads[h + 1] = NAN;
                ++heads_found;
                char title[64];
                snprintf(title, 64, "%d / 8 heads found", heads_found);
                SDL_SetWindowTitle(window, title);
            }
            /* Draw the remaining heads */
            if (!isnan(heads[h])) {
                b3d_reset();
                b3d_rotate_y(h + t * 3);
                b3d_scale(0.4, 0.4, 0.4);
                b3d_translate(x, 0.4 + sinf(h + t * 3) * .1, z);
                srand(h);
                for (int i = 0; i < mesh.vertex_count; i += 9) {
                    uint32_t r = 200 + (RND * 50);
                    uint32_t g = 150 + (RND * 50);
                    uint32_t b = 50 + (RND * 50);
                    b3d_triangle(mesh.triangles[i + 0], mesh.triangles[i + 1],
                                 mesh.triangles[i + 2], mesh.triangles[i + 3],
                                 mesh.triangles[i + 4], mesh.triangles[i + 5],
                                 mesh.triangles[i + 6], mesh.triangles[i + 7],
                                 mesh.triangles[i + 8], (r << 16 | g << 8 | b));
                }
            }
        }

        /* Make a jagged border around the world using stretched cubes */
        srand(seed);
        for (int i = -world_size; i < world_size; i += 2) {
            for (int j = 0; j < 4; ++j) {
                float x = i, z = i;
                switch (j) {
                case 0:
                    x = -world_size;
                    break;
                case 1:
                    x = world_size;
                    break;
                case 2:
                    z = -world_size;
                    break;
                case 3:
                    z = world_size;
                    break;
                }
                b3d_reset();
                b3d_rotate_y(RND * 3.14159f);
                b3d_rotate_x(RND * 3.14159f);
                b3d_rotate_z(RND * 3.14159f);
                b3d_scale(1.0f + RND * 2, 1.0f + RND * 8, 1.0f + RND * 2);
                b3d_translate(x, 0.5, z);
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
        }

        /* Scatter some pyramids around in the world */
        srand(seed);
        for (int i = 0; i < 20; ++i) {
            b3d_reset();
            b3d_scale(1, 1 + RND * 3, 1);
            b3d_rotate_y(RND * 3.14159f);
            b3d_translate((RND * world_size * 2) - world_size, 0,
                          (RND * world_size * 2) - world_size);
            b3d_triangle(0.0, 2.0, 0.0, -1.0, 0.0, 1.0, 1.0, 0.0, 1.0,
                         0x004749);
            b3d_triangle(0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, -1.0,
                         0x00535a);
            b3d_triangle(0.0, 2.0, 0.0, 1.0, 0.0, -1.0, -1.0, 0.0, -1.0,
                         0x00746b);
            b3d_triangle(0.0, 2.0, 0.0, -1.0, 0.0, -1.0, -1.0, 0.0, 1.0,
                         0x00945c);
        }

        /* Draw confetti-like random triangles if all heads have been found */
        if (heads_found == 8) {
            srand(seed);
            for (int i = 0; i < 1000; ++i) {
                b3d_reset();
                b3d_scale(RND * 5, RND * 5, RND * 5);
                b3d_translate((RND * world_size * 2) - world_size,
                              (RND * world_size * 2) - world_size,
                              (RND * world_size * 2) - world_size);
                float y = RND * head_confetti_y;
                y = fmodf(y, 50);
                b3d_triangle(RND - .5, RND - .5 + y, RND - .5, RND - .5,
                             RND - .5 + y, RND - .5, RND - .5, RND - .5 + y,
                             RND - .5, RND * 0xffffff);
            }
            head_confetti_y += 0.1f;
        }

        /* Reset the depth buffer to draw the UI on top of the world */
        memset(b3d_depth, 0x7f, b3d_width * b3d_height * sizeof(b3d_depth[0]));

        /* Draw some UI to show how many heads have been collected */
        for (int h = 0; h < 8; ++h) {
            b3d_reset();
            b3d_set_camera(0, 0, 0, 0, 0, 0);
            b3d_scale(.05, .05, .05);
            b3d_rotate_y(t);
            b3d_translate(h * .1 - .35, -.4 + sinf(h + t * 5) * 0.01, .5);
            srand(h);
            for (int i = 0; i < mesh.vertex_count; i += 9) {
                uint32_t r = 200 + (RND * 50);
                uint32_t g = 150 + (RND * 50);
                uint32_t b = 50 + (RND * 50);
                b3d_triangle(
                    mesh.triangles[i + 0], mesh.triangles[i + 1],
                    mesh.triangles[i + 2], mesh.triangles[i + 3],
                    mesh.triangles[i + 4], mesh.triangles[i + 5],
                    mesh.triangles[i + 6], mesh.triangles[i + 7],
                    mesh.triangles[i + 8],
                    h < heads_found ? (r << 16 | g << 8 | b) : 0x444444);
            }
        }

        /* Draw crosshair */
        b3d_pixels[(b3d_width / 2 - 5) + (b3d_height / 2) * b3d_width] =
            0xffffff;
        b3d_pixels[(b3d_width / 2) + (b3d_height / 2) * b3d_width] = 0xffffff;
        b3d_pixels[(b3d_width / 2 + 5) + (b3d_height / 2) * b3d_width] =
            0xffffff;
        b3d_pixels[(b3d_width / 2) + (b3d_height / 2 - 5) * b3d_width] =
            0xffffff;
        b3d_pixels[(b3d_width / 2) + (b3d_height / 2 + 5) * b3d_width] =
            0xffffff;

        /* Display the pixel buffer on the screen */
        SDL_Delay(1);
        SDL_RenderClear(renderer);
        SDL_UpdateTexture(texture, NULL, pixels, width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }

    free(pixels);
    free(depth);
    b3d_free_mesh(&mesh);
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
