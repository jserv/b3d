/*
 * OBJ model viewer with file browser
 *
 * Loads and renders Wavefront .obj models from the assets/ directory.
 * Supports automatic fan triangulation for quads and n-gons.
 *
 * Controls:
 *   Left/Right arrows - Switch between OBJ files in assets/
 *   ESC or Q          - Quit
 *
 * Command line:
 *   ./obj [file.obj]           - Load specific file
 *   ./obj --snapshot=PATH      - Headless PNG output
 */

#include <SDL.h>
#include <ctype.h>
#include <dirent.h>
#include <stdlib.h>
#include <string.h>

#include "b3d-obj.h"
#include "b3d.h"
#include "utils.h"

#define MAX_OBJ_FILES 64
#define ASSETS_DIR "assets/"

/* Global state for mesh and camera */
static b3d_mesh_t mesh;
static float min_y, max_y, max_xz;
static float y_offset, z_offset;

/* File list */
static char *obj_files[MAX_OBJ_FILES];
static int obj_file_count = 0;
static int current_file = 0;

/* Case-insensitive extension check */
static int has_obj_extension(const char *name, size_t len)
{
    if (len <= 4)
        return 0;
    const char *ext = name + len - 4;
    return (ext[0] == '.') && (tolower(ext[1]) == 'o') &&
           (tolower(ext[2]) == 'b') && (tolower(ext[3]) == 'j');
}

/* Scan assets directory for .obj files */
static void scan_obj_files(void)
{
    DIR *dir = opendir(ASSETS_DIR);
    if (!dir)
        return;

    struct dirent *entry;
    while ((entry = readdir(dir)) && obj_file_count < MAX_OBJ_FILES) {
        const char *name = entry->d_name;
        size_t len = strlen(name);

        /* Check for .obj/.OBJ extension (case-insensitive) */
        if (has_obj_extension(name, len)) {
            /* Build full path: assets/filename.obj */
            size_t path_len = strlen(ASSETS_DIR) + len + 1;
            char *path = malloc(path_len);
            if (!path)
                continue; /* Skip on allocation failure */
            snprintf(path, path_len, "%s%s", ASSETS_DIR, name);
            obj_files[obj_file_count] = path;
            obj_file_count++;
        }
    }
    closedir(dir);

    /* Sort alphabetically for consistent ordering */
    for (int i = 0; i < obj_file_count - 1; i++) {
        for (int j = i + 1; j < obj_file_count; j++) {
            if (strcmp(obj_files[i], obj_files[j]) > 0) {
                char *tmp = obj_files[i];
                obj_files[i] = obj_files[j];
                obj_files[j] = tmp;
            }
        }
    }
}

/* Free file list */
static void free_obj_files(void)
{
    for (int i = 0; i < obj_file_count; i++)
        free(obj_files[i]);
    obj_file_count = 0;
}

/* Load mesh and setup camera. Returns 0 on success. */
static int load_mesh(const char *path)
{
    b3d_mesh_t new_mesh;
    int err = b3d_load_obj(path, &new_mesh);

    if (err != 0) {
        printf("Failed to load '%s' (error %d)\n", path, err);
        return err;
    }

    /* Free previous mesh if any */
    b3d_free_mesh(&mesh);
    mesh = new_mesh;

    /* Calculate bounds for camera positioning */
    b3d_mesh_bounds(&mesh, &min_y, &max_y, &max_xz);

    /* Center model and position camera */
    y_offset = (min_y + max_y) / 2;
    z_offset = -((max_y - min_y) + max_xz);
    if (z_offset > -1.0f)
        z_offset = -1.0f; /* Minimum distance */

    b3d_set_camera(&(b3d_camera_t) {0, y_offset, z_offset, 0, 0, 0});

    printf("[%d/%d] %s: %d triangles\n", current_file + 1, obj_file_count, path,
           mesh.triangle_count);

    return 0;
}

/* Update window title with current file info */
static void update_title(SDL_Window *window)
{
    char title[256];
    const char *filename = strrchr(obj_files[current_file], '/');
    filename = filename ? filename + 1 : obj_files[current_file];

    snprintf(title, sizeof(title),
             "OBJ Viewer [%d/%d] %s (%d tris) - Left/Right to switch",
             current_file + 1, obj_file_count, filename, mesh.triangle_count);
    SDL_SetWindowTitle(window, title);
}

/* Render the current mesh */
static void render_mesh(float time)
{
    b3d_clear();
    b3d_reset();
    b3d_rotate_y(time * 0.3f);

    float height_range = max_y - min_y;
    if (height_range < 0.001f)
        height_range = 1.0f;

    for (int i = 0; i < mesh.vertex_count; i += 9) {
        float avg_y = (mesh.triangles[i + 1] + mesh.triangles[i + 4] +
                       mesh.triangles[i + 7]) /
                      3;
        float brightness = (avg_y - min_y) / height_range;
        uint32_t c = (50 + (int) (brightness * 200)) & 0xff;

        b3d_triangle(
            &(b3d_tri_t) {{{mesh.triangles[i + 0], mesh.triangles[i + 1],
                            mesh.triangles[i + 2]},
                           {mesh.triangles[i + 3], mesh.triangles[i + 4],
                            mesh.triangles[i + 5]},
                           {mesh.triangles[i + 6], mesh.triangles[i + 7],
                            mesh.triangles[i + 8]}}},
            (c << 16 | c << 8 | c));
    }
}

int main(int argument_count, char **arguments)
{
    int width = 800, height = 600;
    int ret = 1;
    const char *snapshot = get_snapshot_path(argument_count, arguments);
    const char *explicit_file = NULL;

    uint32_t *pixels = NULL;
    b3d_depth_t *depth = NULL;
    SDL_Window *window = NULL;
    SDL_Renderer *renderer = NULL;
    SDL_Texture *texture = NULL;

    /* Check for explicit file argument */
    for (int i = 1; i < argument_count; ++i) {
        if (strncmp(arguments[i], "--snapshot=", 11) != 0) {
            explicit_file = arguments[i];
            break;
        }
    }

    /* Scan for available OBJ files */
    scan_obj_files();

    if (obj_file_count == 0 && !explicit_file) {
        printf("No .obj files found in %s\n", ASSETS_DIR);
        goto cleanup_files;
    }

    /* If explicit file given, find it in the list or use directly */
    if (explicit_file) {
        int found = 0;
        for (int i = 0; i < obj_file_count; i++) {
            if (strcmp(obj_files[i], explicit_file) == 0 ||
                strcmp(strrchr(obj_files[i], '/') + 1, explicit_file) == 0) {
                current_file = i;
                found = 1;
                break;
            }
        }
        if (!found && obj_file_count < MAX_OBJ_FILES) {
            obj_files[obj_file_count] = strdup(explicit_file);
            current_file = obj_file_count;
            obj_file_count++;
        } else if (!found) {
            printf("Warning: file list full, cannot add '%s'\n", explicit_file);
        }
    }

    pixels = malloc(width * height * sizeof(pixels[0]));
    depth = malloc(width * height * sizeof(depth[0]));
    if (!pixels || !depth) {
        fprintf(stderr, "Failed to allocate framebuffer\n");
        goto cleanup_buffers;
    }

    b3d_init(pixels, depth, width, height, 70);

    if (load_mesh(obj_files[current_file]) != 0)
        goto cleanup_buffers;

    /* Snapshot mode: render single frame and exit */
    if (snapshot) {
        render_mesh(0.8f);
        write_png(snapshot, pixels, width, height);
        ret = 0;
        goto cleanup_mesh;
    }

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        goto cleanup_mesh;
    }

    window = SDL_CreateWindow("", SDL_WINDOWPOS_CENTERED,
                              SDL_WINDOWPOS_CENTERED, width, height, 0);
    if (!window) {
        fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
        goto cleanup_sdl;
    }

    renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC);
    if (!renderer) {
        fprintf(stderr, "SDL_CreateRenderer failed: %s\n", SDL_GetError());
        goto cleanup_window;
    }

    texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888,
                                SDL_TEXTUREACCESS_STREAMING, width, height);
    if (!texture) {
        fprintf(stderr, "SDL_CreateTexture failed: %s\n", SDL_GetError());
        goto cleanup_renderer;
    }

    update_title(window);

    int quit = 0;
    while (!quit) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
                break;
            }
            if (event.type == SDL_KEYDOWN) {
                switch (event.key.keysym.sym) {
                case SDLK_ESCAPE:
                case SDLK_q:
                    quit = 1;
                    break;

                case SDLK_LEFT:
                case SDLK_UP:
                    if (obj_file_count > 1) {
                        int prev = (current_file - 1 + obj_file_count) %
                                   obj_file_count;
                        if (load_mesh(obj_files[prev]) == 0) {
                            current_file = prev;
                            update_title(window);
                        }
                    }
                    break;

                case SDLK_RIGHT:
                case SDLK_DOWN:
                    if (obj_file_count > 1) {
                        int next = (current_file + 1) % obj_file_count;
                        if (load_mesh(obj_files[next]) == 0) {
                            current_file = next;
                            update_title(window);
                        }
                    }
                    break;

                case SDLK_HOME:
                    if (current_file != 0 && load_mesh(obj_files[0]) == 0) {
                        current_file = 0;
                        update_title(window);
                    }
                    break;

                case SDLK_END:
                    if (current_file != obj_file_count - 1) {
                        int last = obj_file_count - 1;
                        if (load_mesh(obj_files[last]) == 0) {
                            current_file = last;
                            update_title(window);
                        }
                    }
                    break;
                }
            }
        }
        if (quit)
            break;

        float t = SDL_GetTicks() * 0.001f;
        render_mesh(t);

        SDL_RenderClear(renderer);
        SDL_UpdateTexture(texture, NULL, pixels, width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }

    ret = 0;

    SDL_DestroyTexture(texture);
cleanup_renderer:
    SDL_DestroyRenderer(renderer);
cleanup_window:
    SDL_DestroyWindow(window);
cleanup_sdl:
    SDL_Quit();
cleanup_mesh:
    b3d_free_mesh(&mesh);
cleanup_buffers:
    free(pixels);
    free(depth);
cleanup_files:
    free_obj_files();
    return ret;
}
