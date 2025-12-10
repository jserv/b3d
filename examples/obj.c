/*
    Loads a model from a '.obj' file and renders it using bootleg3d.
    You can pass the path to a file as the first argument, otherwise
    it will load a file (moai.obj) provided in this examples folder.

    Right now the '.obj' loading code will fail on models containing
    non-triangular faces.
*/

#include "SDL.h"
#include "b3d.h"
#include "b3d_obj.h"

int main(int argument_count, char ** arguments) {
    int width = 800;
    int height = 600;

    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window * window = SDL_CreateWindow("", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, width, height, 0);
    SDL_Renderer * renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC);
    SDL_Texture * texture = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888, SDL_TEXTUREACCESS_STREAMING, width, height);

    uint32_t * pixels = malloc(width * height * sizeof(pixels[0]));
    b3d_depth_t * depth = malloc(width * height * sizeof(depth[0]));
    b3d_init(pixels, depth, width, height, 70);

    /* Load OBJ file */
    const char * file_name = "assets/moai.obj";
    if (argument_count == 2) file_name = arguments[1];

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

    printf("Loaded %d triangles from file '%s'.\n", mesh.triangle_count, file_name);

    /* Calculate mesh bounds for camera positioning */
    float min_y, max_y, max_xz;
    b3d_mesh_bounds(&mesh, &min_y, &max_y, &max_xz);

    /* Centre the model along the y axis */
    float y_offset = (min_y + max_y) / 2;

    /* Move backwards based on how tall and wide the model is */
    float z_offset = -((max_y - min_y) + max_xz);

    b3d_set_camera(0, y_offset, z_offset, 0, 0, 0);

    int quit = 0;
    while (!quit) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
                break;
            }
        }
        if (quit) break;

        float t = SDL_GetTicks() * 0.001f;
        b3d_clear();
        b3d_reset();
        b3d_rotate_y(t * 0.3);

        for (int i = 0; i < mesh.vertex_count; i += 9) {
            float avg_y = (mesh.triangles[i + 1] + mesh.triangles[i + 4] + mesh.triangles[i + 7]) / 3;
            float brightness = (avg_y - min_y) / max_y;
            uint32_t c = (50 + (int)(brightness * 200)) & 0xff;
            b3d_triangle(
                mesh.triangles[i + 0], mesh.triangles[i + 1], mesh.triangles[i + 2],
                mesh.triangles[i + 3], mesh.triangles[i + 4], mesh.triangles[i + 5],
                mesh.triangles[i + 6], mesh.triangles[i + 7], mesh.triangles[i + 8],
                (c << 16 | c << 8 | c)
            );
        }

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
