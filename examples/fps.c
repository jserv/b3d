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

#include "b3d-math.h"
#include "b3d-obj.h"
#include "b3d.h"
#include "utils.h"

#define RND (rand() / (float) RAND_MAX)

int main(int argument_count, char **arguments)
{
    int width = 640, height = 480;
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
        b3d_set_camera(&(b3d_camera_t) {-1.0f, 1.0f, -3.0f, 1.0f, 0.0f, 0});

        /* Draw checkerboard floor */
        b3d_reset();
        for (int z = -world_size; z < world_size; ++z) {
            for (int x = -world_size; x < world_size; ++x) {
                uint32_t c = (x + z) & 1 ? 0x424C88 : 0xF7C396;
                b3d_triangle(&(b3d_tri_t) {{{x + .5f, 0, z + .5f},
                                            {x - .5f, 0, z - .5f},
                                            {x - .5f, 0, z + .5f}}},
                             c);
                b3d_triangle(&(b3d_tri_t) {{{x + .5f, 0, z + .5f},
                                            {x + .5f, 0, z - .5f},
                                            {x - .5f, 0, z - .5f}}},
                             c);
            }
        }

        /* Draw golden heads */
        float t = 1.5f;
        for (int h = 0; h < 16; h += 2) {
            float hx = heads[h], hz = heads[h + 1];
            b3d_reset();
            b3d_rotate_y(h + t * 3);
            b3d_scale(0.4, 0.4, 0.4);
            b3d_translate(hx, 0.4 + b3d_sinf(h + t * 3) * 0.1, hz);
            srand(h);
            for (int i = 0; i < mesh.vertex_count; i += 9) {
                uint32_t r = 200 + (int) (RND * 50);
                uint32_t g = 150 + (int) (RND * 50);
                uint32_t bv = 50 + (int) (RND * 50);
                b3d_triangle(
                    &(b3d_tri_t) {
                        {{mesh.triangles[i + 0], mesh.triangles[i + 1],
                          mesh.triangles[i + 2]},
                         {mesh.triangles[i + 3], mesh.triangles[i + 4],
                          mesh.triangles[i + 5]},
                         {mesh.triangles[i + 6], mesh.triangles[i + 7],
                          mesh.triangles[i + 8]}}},
                    (r << 16 | g << 8 | bv));
            }
        }

        /* Cube and pyramid geometry */
        static const b3d_tri_t cube_faces[12] = {
            {
                {
                    {-0.5, -0.5, -0.5},
                    {-0.5, 0.5, -0.5},
                    {0.5, 0.5, -0.5},
                },
            },
            {
                {
                    {-0.5, -0.5, -0.5},
                    {0.5, 0.5, -0.5},
                    {0.5, -0.5, -0.5},
                },
            },
            {
                {
                    {0.5, -0.5, -0.5},
                    {0.5, 0.5, -0.5},
                    {0.5, 0.5, 0.5},
                },
            },
            {
                {
                    {0.5, -0.5, -0.5},
                    {0.5, 0.5, 0.5},
                    {0.5, -0.5, 0.5},
                },
            },
            {
                {
                    {0.5, -0.5, 0.5},
                    {0.5, 0.5, 0.5},
                    {-0.5, 0.5, 0.5},
                },
            },
            {
                {
                    {0.5, -0.5, 0.5},
                    {-0.5, 0.5, 0.5},
                    {-0.5, -0.5, 0.5},
                },
            },
            {
                {
                    {-0.5, -0.5, 0.5},
                    {-0.5, 0.5, 0.5},
                    {-0.5, 0.5, -0.5},
                },
            },
            {
                {
                    {-0.5, -0.5, 0.5},
                    {-0.5, 0.5, -0.5},
                    {-0.5, -0.5, -0.5},
                },
            },
            {
                {
                    {-0.5, 0.5, -0.5},
                    {-0.5, 0.5, 0.5},
                    {0.5, 0.5, 0.5},
                },
            },
            {
                {
                    {-0.5, 0.5, -0.5},
                    {0.5, 0.5, 0.5},
                    {0.5, 0.5, -0.5},
                },
            },
            {
                {
                    {0.5, -0.5, 0.5},
                    {-0.5, -0.5, 0.5},
                    {-0.5, -0.5, -0.5},
                },
            },
            {
                {
                    {0.5, -0.5, 0.5},
                    {-0.5, -0.5, -0.5},
                    {0.5, -0.5, -0.5},
                },
            },
        };
        static const uint32_t cube_colors[12] = {
            0xfcd0a1, 0xb1b695, 0x53917e, 0x63535b, 0x6d1a36, 0xd4e09b,
            0xf6f4d2, 0xcbdfbd, 0xf19c79, 0xa44a3f, 0x5465ff, 0x788bff,
        };
        static const b3d_tri_t pyramid_faces[4] = {
            {
                {
                    {0.0, 2.0, 0.0},
                    {-1.0, 0.0, 1.0},
                    {1.0, 0.0, 1.0},
                },
            },
            {
                {
                    {0.0, 2.0, 0.0},
                    {1.0, 0.0, 1.0},
                    {1.0, 0.0, -1.0},
                },
            },
            {
                {
                    {0.0, 2.0, 0.0},
                    {1.0, 0.0, -1.0},
                    {-1.0, 0.0, -1.0},
                },
            },
            {
                {
                    {0.0, 2.0, 0.0},
                    {-1.0, 0.0, -1.0},
                    {-1.0, 0.0, 1.0},
                },
            },
        };
        static const uint32_t pyramid_colors[4] = {0x004749, 0x00535a, 0x00746b,
                                                   0x00945c};

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
                for (int f = 0; f < 12; ++f)
                    b3d_triangle(&cube_faces[f], cube_colors[f]);
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
            for (int f = 0; f < 4; ++f)
                b3d_triangle(&pyramid_faces[f], pyramid_colors[f]);
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
        float sin_fwd, cos_fwd, sin_str, cos_str;
        b3d_sincosf(player_yaw - 1.570796f, &sin_fwd, &cos_fwd);
        b3d_sincosf(player_yaw, &sin_str, &cos_str);
        player_x -= cos_fwd * player_forward_speed;
        player_z -= sin_fwd * player_forward_speed;
        player_x -= cos_str * player_strafe_speed;
        player_z -= sin_str * player_strafe_speed;
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

        b3d_set_camera(&(b3d_camera_t) {player_x, player_height, player_z,
                                        player_yaw, player_pitch, 0});

        /* Draw a checkerboard floor */
        b3d_reset();
        for (int z = -world_size; z < world_size; ++z) {
            for (int x = -world_size; x < world_size; ++x) {
                uint32_t c = (x + z) & 1 ? 0x424C88 : 0xF7C396;
                b3d_triangle(&(b3d_tri_t) {{{x + .5f, 0, z + .5f},
                                            {x - .5f, 0, z - .5f},
                                            {x - .5f, 0, z + .5f}}},
                             c);
                b3d_triangle(&(b3d_tri_t) {{{x + .5f, 0, z + .5f},
                                            {x + .5f, 0, z - .5f},
                                            {x - .5f, 0, z - .5f}}},
                             c);
            }
        }

        /* Draw the golden heads */
        for (int h = 0; h < 16; h += 2) {
            float x = heads[h];
            float z = heads[h + 1];
            /* If the player is close to a head, mark it as found using NAN */
            if (b3d_fabsf(player_x - x) < head_radius &&
                b3d_fabsf(player_z - z) < head_radius) {
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
                b3d_translate(x, 0.4 + b3d_sinf(h + t * 3) * .1, z);
                srand(h);
                for (int i = 0; i < mesh.vertex_count; i += 9) {
                    uint32_t r = 200 + (RND * 50);
                    uint32_t g = 150 + (RND * 50);
                    uint32_t bv = 50 + (RND * 50);
                    b3d_triangle(
                        &(b3d_tri_t) {
                            {{mesh.triangles[i + 0], mesh.triangles[i + 1],
                              mesh.triangles[i + 2]},
                             {mesh.triangles[i + 3], mesh.triangles[i + 4],
                              mesh.triangles[i + 5]},
                             {mesh.triangles[i + 6], mesh.triangles[i + 7],
                              mesh.triangles[i + 8]}}},
                        (r << 16 | g << 8 | bv));
                }
            }
        }

        /* Cube and pyramid geometry */
        static const b3d_tri_t cube_faces[12] = {
            {
                {
                    {-0.5, -0.5, -0.5},
                    {-0.5, 0.5, -0.5},
                    {0.5, 0.5, -0.5},
                },
            },
            {
                {
                    {-0.5, -0.5, -0.5},
                    {0.5, 0.5, -0.5},
                    {0.5, -0.5, -0.5},
                },
            },
            {
                {
                    {0.5, -0.5, -0.5},
                    {0.5, 0.5, -0.5},
                    {0.5, 0.5, 0.5},
                },
            },
            {
                {
                    {0.5, -0.5, -0.5},
                    {0.5, 0.5, 0.5},
                    {0.5, -0.5, 0.5},
                },
            },
            {
                {
                    {0.5, -0.5, 0.5},
                    {0.5, 0.5, 0.5},
                    {-0.5, 0.5, 0.5},
                },
            },
            {
                {
                    {0.5, -0.5, 0.5},
                    {-0.5, 0.5, 0.5},
                    {-0.5, -0.5, 0.5},
                },
            },
            {
                {
                    {-0.5, -0.5, 0.5},
                    {-0.5, 0.5, 0.5},
                    {-0.5, 0.5, -0.5},
                },
            },
            {
                {
                    {-0.5, -0.5, 0.5},
                    {-0.5, 0.5, -0.5},
                    {-0.5, -0.5, -0.5},
                },
            },
            {
                {
                    {-0.5, 0.5, -0.5},
                    {-0.5, 0.5, 0.5},
                    {0.5, 0.5, 0.5},
                },
            },
            {
                {
                    {-0.5, 0.5, -0.5},
                    {0.5, 0.5, 0.5},
                    {0.5, 0.5, -0.5},
                },
            },
            {
                {
                    {0.5, -0.5, 0.5},
                    {-0.5, -0.5, 0.5},
                    {-0.5, -0.5, -0.5},
                },
            },
            {
                {
                    {0.5, -0.5, 0.5},
                    {-0.5, -0.5, -0.5},
                    {0.5, -0.5, -0.5},
                },
            },
        };
        static const uint32_t cube_colors[12] = {
            0xfcd0a1, 0xb1b695, 0x53917e, 0x63535b, 0x6d1a36, 0xd4e09b,
            0xf6f4d2, 0xcbdfbd, 0xf19c79, 0xa44a3f, 0x5465ff, 0x788bff,
        };
        static const b3d_tri_t pyramid_faces[4] = {
            {
                {
                    {0.0, 2.0, 0.0},
                    {-1.0, 0.0, 1.0},
                    {1.0, 0.0, 1.0},
                },
            },
            {
                {
                    {0.0, 2.0, 0.0},
                    {1.0, 0.0, 1.0},
                    {1.0, 0.0, -1.0},
                },
            },
            {
                {
                    {0.0, 2.0, 0.0},
                    {1.0, 0.0, -1.0},
                    {-1.0, 0.0, -1.0},
                },
            },
            {
                {
                    {0.0, 2.0, 0.0},
                    {-1.0, 0.0, -1.0},
                    {-1.0, 0.0, 1.0},
                },
            },
        };
        static const uint32_t pyramid_colors[4] = {0x004749, 0x00535a, 0x00746b,
                                                   0x00945c};

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
                for (int f = 0; f < 12; ++f)
                    b3d_triangle(&cube_faces[f], cube_colors[f]);
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
            for (int f = 0; f < 4; ++f)
                b3d_triangle(&pyramid_faces[f], pyramid_colors[f]);
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
                b3d_triangle(
                    &(b3d_tri_t) {{{RND - .5f, RND - .5f + y, RND - .5f},
                                   {RND - .5f, RND - .5f + y, RND - .5f},
                                   {RND - .5f, RND - .5f + y, RND - .5f}}},
                    (uint32_t) (RND * 0xffffff));
            }
            head_confetti_y += 0.1f;
        }

        /* Reset the depth buffer to draw the UI on top of the world */
        memset(depth, 0x7f, width * height * sizeof(depth[0]));

        /* Draw some UI to show how many heads have been collected */
        for (int h = 0; h < 8; ++h) {
            b3d_reset();
            b3d_set_camera(&(b3d_camera_t) {0, 0, 0, 0, 0, 0});
            b3d_scale(.05, .05, .05);
            b3d_rotate_y(t);
            b3d_translate(h * .1 - .35, -.4 + b3d_sinf(h + t * 5) * 0.01, .5);
            srand(h);
            for (int i = 0; i < mesh.vertex_count; i += 9) {
                uint32_t r = 200 + (uint32_t) (RND * 50);
                uint32_t g = 150 + (uint32_t) (RND * 50);
                uint32_t bv = 50 + (uint32_t) (RND * 50);
                b3d_triangle(
                    &(b3d_tri_t) {
                        {{mesh.triangles[i + 0], mesh.triangles[i + 1],
                          mesh.triangles[i + 2]},
                         {mesh.triangles[i + 3], mesh.triangles[i + 4],
                          mesh.triangles[i + 5]},
                         {mesh.triangles[i + 6], mesh.triangles[i + 7],
                          mesh.triangles[i + 8]}}},
                    h < heads_found ? (r << 16 | g << 8 | bv) : 0x444444);
            }
        }

        /* Draw crosshair */
        int cx = width / 2, cy = height / 2;
        pixels[(cx - 5) + cy * width] = 0xffffff;
        pixels[cx + cy * width] = 0xffffff;
        pixels[(cx + 5) + cy * width] = 0xffffff;
        pixels[cx + (cy - 5) * width] = 0xffffff;
        pixels[cx + (cy + 5) * width] = 0xffffff;

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
