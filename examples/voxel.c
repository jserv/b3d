/*
 * Mesh voxelizer demo
 *
 * Loads an OBJ model and converts it to voxels, then renders the voxelized
 * version. Demonstrates the b3d-voxel.h voxelization functionality.
 *
 * Usage: ./voxel [model.obj] [--voxel-size=N]
 *
 * Supports headless snapshots with --snapshot=PATH or B3D_SNAPSHOT.
 */

#include <SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "b3d-obj.h"
#include "b3d-voxel.h"
#include "b3d.h"
#include "utils.h"

/* Default voxel size (smaller = more voxels, higher resolution) */
#define DEFAULT_VOXEL_SIZE 0.05f

/* Maximum voxels to allocate */
#define MAX_VOXELS 100000

/* Parse voxel size from command line */
static float get_voxel_size(int argc, char **argv)
{
    for (int i = 1; i < argc; ++i) {
        const char *flag = "--voxel-size=";
        size_t len = strlen(flag);
        if (!strncmp(argv[i], flag, len)) {
            const char *str = argv[i] + len;
            char *endptr;
            float size = strtof(str, &endptr);
            /* Validate: consumed input, no error, within valid range */
            if (endptr != str && *endptr == '\0' && size > 0.001f &&
                size < 1.0f)
                return size;
            fprintf(stderr, "Warning: invalid voxel size '%s', using default\n",
                    str);
        }
    }
    return DEFAULT_VOXEL_SIZE;
}

/* Get OBJ file path from arguments */
static const char *get_obj_path(int argc, char **argv)
{
    for (int i = 1; i < argc; ++i) {
        if (strncmp(argv[i], "--", 2) != 0)
            return argv[i];
    }
    return "assets/moai.obj";
}

int main(int argc, char **argv)
{
    int width = 800, height = 600;
    const char *snapshot = get_snapshot_path(argc, argv);
    const char *obj_path = get_obj_path(argc, argv);
    float voxel_size = get_voxel_size(argc, argv);

    uint32_t *pixels = malloc(width * height * sizeof(pixels[0]));
    b3d_depth_t *depth = malloc(width * height * sizeof(depth[0]));

    if (!pixels || !depth) {
        fprintf(stderr, "Failed to allocate framebuffer\n");
        return 1;
    }

    /* Load OBJ file */
    b3d_mesh_t mesh;
    int err = b3d_load_obj(obj_path, &mesh);
    if (err) {
        fprintf(stderr, "Failed to load '%s' (error %d)\n", obj_path, err);
        free(pixels);
        free(depth);
        return 1;
    }

    printf("Loaded %d triangles from '%s'\n", mesh.triangle_count, obj_path);
    printf("Voxel size: %.3f\n", voxel_size);

    /* Allocate voxel buffer */
    b3d_voxel_t *voxels = malloc(MAX_VOXELS * sizeof(b3d_voxel_t));
    if (!voxels) {
        fprintf(stderr, "Failed to allocate voxel buffer\n");
        b3d_free_mesh(&mesh);
        free(pixels);
        free(depth);
        return 1;
    }

    /* Voxelize the mesh */
    printf("Voxelizing...\n");
    size_t voxel_count =
        b3d_voxelize(mesh.triangles, (size_t) mesh.triangle_count, voxel_size,
                     0x8080FF, /* Light blue */
                     voxels, MAX_VOXELS);

    /* Check for truncation (returns max + 1 when buffer full) */
    if (voxel_count > MAX_VOXELS) {
        printf("Warning: voxel buffer full, output truncated to %d voxels\n",
               MAX_VOXELS);
        voxel_count = MAX_VOXELS;
    } else {
        printf("Generated %zu voxels\n", voxel_count);
    }

    /* Calculate bounds for camera positioning */
    float min_y, max_y, max_xz;
    b3d_mesh_bounds(&mesh, &min_y, &max_y, &max_xz);

    float y_offset = (min_y + max_y) / 2;
    float z_offset = -((max_y - min_y) + max_xz) * 1.2f;

    b3d_init(pixels, depth, width, height, 70);
    b3d_set_camera(&(b3d_camera_t) {0, y_offset, z_offset, 0, 0, 0});
    b3d_set_light_direction(0.5f, 1.0f, 0.7f);
    b3d_set_ambient(0.3f);

    if (snapshot) {
        float t = 0.8f;
        b3d_clear();
        b3d_reset();
        b3d_rotate_y(t * 0.5f);

        b3d_voxel_render(voxels, voxel_count, voxel_size);

        write_png(snapshot, pixels, width, height);
        printf("Snapshot saved to '%s'\n", snapshot);

        free(voxels);
        free(pixels);
        free(depth);
        b3d_free_mesh(&mesh);
        return 0;
    }

    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        fprintf(stderr, "SDL_Init failed: %s\n", SDL_GetError());
        free(voxels);
        free(pixels);
        free(depth);
        b3d_free_mesh(&mesh);
        return 1;
    }

    SDL_Window *window =
        SDL_CreateWindow("B3D Voxelizer", SDL_WINDOWPOS_CENTERED,
                         SDL_WINDOWPOS_CENTERED, width, height, 0);
    if (!window) {
        fprintf(stderr, "SDL_CreateWindow failed: %s\n", SDL_GetError());
        SDL_Quit();
        free(voxels);
        free(pixels);
        free(depth);
        b3d_free_mesh(&mesh);
        return 1;
    }

    SDL_Renderer *renderer =
        SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC);
    if (!renderer) {
        fprintf(stderr, "SDL_CreateRenderer failed: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        free(voxels);
        free(pixels);
        free(depth);
        b3d_free_mesh(&mesh);
        return 1;
    }

    SDL_Texture *texture =
        SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888,
                          SDL_TEXTUREACCESS_STREAMING, width, height);
    if (!texture) {
        fprintf(stderr, "SDL_CreateTexture failed: %s\n", SDL_GetError());
        SDL_DestroyRenderer(renderer);
        SDL_DestroyWindow(window);
        SDL_Quit();
        free(voxels);
        free(pixels);
        free(depth);
        b3d_free_mesh(&mesh);
        return 1;
    }

    printf("Controls: ESC to quit, Space to toggle lighting\n");

    int quit = 0;
    int use_lighting = 1;

    while (!quit) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT)
                quit = 1;
            else if (event.type == SDL_KEYDOWN) {
                if (event.key.keysym.sym == SDLK_ESCAPE)
                    quit = 1;
                else if (event.key.keysym.sym == SDLK_SPACE)
                    use_lighting = !use_lighting;
            }
        }

        float t = SDL_GetTicks() * 0.001f;
        b3d_clear();
        b3d_reset();
        b3d_rotate_y(t * 0.5f);

        if (use_lighting)
            b3d_voxel_render(voxels, voxel_count, voxel_size);
        else
            b3d_voxel_render_flat(voxels, voxel_count, voxel_size);

        SDL_RenderClear(renderer);
        SDL_UpdateTexture(texture, NULL, pixels, width * sizeof(uint32_t));
        SDL_RenderCopy(renderer, texture, NULL, NULL);
        SDL_RenderPresent(renderer);
    }

    free(voxels);
    free(pixels);
    free(depth);
    b3d_free_mesh(&mesh);
    SDL_DestroyTexture(texture);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}
