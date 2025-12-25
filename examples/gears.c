/*
 * Gears demo
 *
 * Based on the classic glxgears by Brian Paul (1999).
 * Supports headless snapshots with --snapshot=PATH or B3D_SNAPSHOT.
 */

#include <SDL.h>
#include <stdio.h>
#include <stdlib.h>

#include "b3d-math.h"
#include "b3d.h"
#include "utils.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

/* Convert degrees to radians */
#define DEG2RAD(x) ((x) * M_PI / 180.0f)

/* Draw a gear
 * @inner_radius: radius of hole at center
 * @outer_radius: radius at center of teeth
 * @width: width of gear
 * @teeth: number of teeth
 * @tooth_depth: depth of tooth
 * @color: base color in 0xRRGGBB format
 *
 * Generates gear geometry matching the original OpenGL gears demo.
 */
static void gear(float inner_radius,
                 float outer_radius,
                 float width,
                 int teeth,
                 float tooth_depth,
                 uint32_t color)
{
    float r0 = inner_radius;
    float r1 = outer_radius - tooth_depth / 2.0f;
    float r2 = outer_radius + tooth_depth / 2.0f;

    float da = 2.0f * M_PI / (float) teeth / 4.0f;

    /* Draw front face (ring + tooth fronts) - normal points +Z */
    for (int i = 0; i < teeth; i++) {
        float angle = (float) i * 2.0f * M_PI / (float) teeth;

        float s0, c0, s1, c1, s2, c2, s3, c3, s4, c4;
        b3d_sincosf(angle, &s0, &c0);
        b3d_sincosf(angle + 1.0f * da, &s1, &c1);
        b3d_sincosf(angle + 2.0f * da, &s2, &c2);
        b3d_sincosf(angle + 3.0f * da, &s3, &c3);
        b3d_sincosf(angle + 4.0f * da, &s4, &c4);

        /* Ring sections (four per tooth, mirroring glxgears quad strip) */
        float cs[] = {c0, c1, c2, c3};
        float ss[] = {s0, s1, s2, s3};
        float ce[] = {c1, c2, c3, c4};
        float se[] = {s1, s2, s3, s4};
        for (int j = 0; j < 4; ++j) {
            b3d_triangle_lit(&(b3d_tri_t) {{
                                 {r0 * cs[j], r0 * ss[j], width * 0.5f},
                                 {r1 * cs[j], r1 * ss[j], width * 0.5f},
                                 {r1 * ce[j], r1 * se[j], width * 0.5f},
                             }},
                             0.0f, 0.0f, 1.0f, color);
            b3d_triangle_lit(&(b3d_tri_t) {{
                                 {r0 * cs[j], r0 * ss[j], width * 0.5f},
                                 {r1 * ce[j], r1 * se[j], width * 0.5f},
                                 {r0 * ce[j], r0 * se[j], width * 0.5f},
                             }},
                             0.0f, 0.0f, 1.0f, color);
        }

        /* Front sides of teeth */
        b3d_triangle_lit(&(b3d_tri_t) {{
                             {r1 * c0, r1 * s0, width * 0.5f},
                             {r2 * c1, r2 * s1, width * 0.5f},
                             {r2 * c2, r2 * s2, width * 0.5f},
                         }},
                         0.0f, 0.0f, 1.0f, color);
        b3d_triangle_lit(&(b3d_tri_t) {{
                             {r1 * c0, r1 * s0, width * 0.5f},
                             {r2 * c2, r2 * s2, width * 0.5f},
                             {r1 * c3, r1 * s3, width * 0.5f},
                         }},
                         0.0f, 0.0f, 1.0f, color);
    }

    /* Draw back face - normal points -Z */
    for (int i = 0; i < teeth; i++) {
        float angle = (float) i * 2.0f * M_PI / (float) teeth;

        float s0, c0, s1, c1, s2, c2, s3, c3, s4, c4;
        b3d_sincosf(angle, &s0, &c0);
        b3d_sincosf(angle + 1.0f * da, &s1, &c1);
        b3d_sincosf(angle + 2.0f * da, &s2, &c2);
        b3d_sincosf(angle + 3.0f * da, &s3, &c3);
        b3d_sincosf(angle + 4.0f * da, &s4, &c4);

        /* Back face triangles - reverse winding */
        float cs[] = {c1, c2, c3, c4};
        float ss[] = {s1, s2, s3, s4};
        float ce[] = {c0, c1, c2, c3};
        float se[] = {s0, s1, s2, s3};
        for (int j = 0; j < 4; ++j) {
            b3d_triangle_lit(&(b3d_tri_t) {{
                                 {r1 * cs[j], r1 * ss[j], -width * 0.5f},
                                 {r1 * ce[j], r1 * se[j], -width * 0.5f},
                                 {r0 * ce[j], r0 * se[j], -width * 0.5f},
                             }},
                             0.0f, 0.0f, -1.0f, color);
            b3d_triangle_lit(&(b3d_tri_t) {{
                                 {r0 * cs[j], r0 * ss[j], -width * 0.5f},
                                 {r1 * cs[j], r1 * ss[j], -width * 0.5f},
                                 {r0 * ce[j], r0 * se[j], -width * 0.5f},
                             }},
                             0.0f, 0.0f, -1.0f, color);
        }

        /* Back sides of teeth - reverse winding */
        b3d_triangle_lit(&(b3d_tri_t) {{
                             {r2 * c1, r2 * s1, -width * 0.5f},
                             {r1 * c0, r1 * s0, -width * 0.5f},
                             {r1 * c3, r1 * s3, -width * 0.5f},
                         }},
                         0.0f, 0.0f, -1.0f, color);
        b3d_triangle_lit(&(b3d_tri_t) {{
                             {r2 * c2, r2 * s2, -width * 0.5f},
                             {r2 * c1, r2 * s1, -width * 0.5f},
                             {r1 * c3, r1 * s3, -width * 0.5f},
                         }},
                         0.0f, 0.0f, -1.0f, color);
    }

    /* Draw outward faces of teeth (lead edge, top, trail edge, valley) */
    for (int i = 0; i < teeth; i++) {
        float angle = (float) i * 2.0f * M_PI / (float) teeth;

        float s0, c0, s1, c1, s2, c2, s3, c3, s4, c4;
        b3d_sincosf(angle, &s0, &c0);
        b3d_sincosf(angle + 1.0f * da, &s1, &c1);
        b3d_sincosf(angle + 2.0f * da, &s2, &c2);
        b3d_sincosf(angle + 3.0f * da, &s3, &c3);
        b3d_sincosf(angle + 4.0f * da, &s4, &c4);

        float hw = width * 0.5f;
        float u, v, len;
        float nx, ny;

        /* Leading edge of tooth */
        u = r2 * c1 - r1 * c0;
        v = r2 * s1 - r1 * s0;
        len = b3d_sqrtf(u * u + v * v);
        if (len > 0.0f) {
            nx = v / len;
            ny = -u / len;
            b3d_triangle_lit(&(b3d_tri_t) {{
                                 {r1 * c0, r1 * s0, hw},
                                 {r1 * c0, r1 * s0, -hw},
                                 {r2 * c1, r2 * s1, -hw},
                             }},
                             nx, ny, 0.0f, color);
            b3d_triangle_lit(&(b3d_tri_t) {{
                                 {r1 * c0, r1 * s0, hw},
                                 {r2 * c1, r2 * s1, -hw},
                                 {r2 * c1, r2 * s1, hw},
                             }},
                             nx, ny, 0.0f, color);
        }

        /* Tooth top face */
        float s_top, c_top;
        b3d_sincosf(angle + 1.5f * da, &s_top, &c_top);
        b3d_triangle_lit(&(b3d_tri_t) {{{r2 * c1, r2 * s1, hw},
                                        {r2 * c1, r2 * s1, -hw},
                                        {r2 * c2, r2 * s2, -hw}}},
                         c_top, s_top, 0.0f, color);
        b3d_triangle_lit(&(b3d_tri_t) {{{r2 * c1, r2 * s1, hw},
                                        {r2 * c2, r2 * s2, -hw},
                                        {r2 * c2, r2 * s2, hw}}},
                         c_top, s_top, 0.0f, color);

        /* Trailing edge of tooth */
        u = r1 * c3 - r2 * c2;
        v = r1 * s3 - r2 * s2;
        len = b3d_sqrtf(u * u + v * v);
        if (len > 0.0f) {
            nx = v / len;
            ny = -u / len;
            b3d_triangle_lit(&(b3d_tri_t) {{
                                 {r2 * c2, r2 * s2, hw},
                                 {r2 * c2, r2 * s2, -hw},
                                 {r1 * c3, r1 * s3, -hw},
                             }},
                             nx, ny, 0.0f, color);
            b3d_triangle_lit(&(b3d_tri_t) {{
                                 {r2 * c2, r2 * s2, hw},
                                 {r1 * c3, r1 * s3, -hw},
                                 {r1 * c3, r1 * s3, hw},
                             }},
                             nx, ny, 0.0f, color);
        }

        /* Valley between teeth (outer rim at r1) */
        float s_rim, c_rim;
        b3d_sincosf(angle + 3.5f * da, &s_rim, &c_rim);
        b3d_triangle_lit(&(b3d_tri_t) {{
                             {r1 * c3, r1 * s3, hw},
                             {r1 * c3, r1 * s3, -hw},
                             {r1 * c4, r1 * s4, -hw},
                         }},
                         c_rim, s_rim, 0.0f, color);
        b3d_triangle_lit(&(b3d_tri_t) {{
                             {r1 * c3, r1 * s3, hw},
                             {r1 * c4, r1 * s4, -hw},
                             {r1 * c4, r1 * s4, hw},
                         }},
                         c_rim, s_rim, 0.0f, color);
    }

    /* Draw inside radius cylinder */
    int segs = teeth * 4;
    for (int i = 0; i < segs; i++) {
        float a0 = (float) i * 2.0f * M_PI / (float) segs;
        float a1 = (float) (i + 1) * 2.0f * M_PI / (float) segs;

        float s0, c0, s1, c1;
        b3d_sincosf(a0, &s0, &c0);
        b3d_sincosf(a1, &s1, &c1);
        float hw = width * 0.5f;

        /* Normal points inward (average of two edge normals) */
        float nx = -(c0 + c1);
        float ny = -(s0 + s1);

        b3d_triangle_lit(&(b3d_tri_t) {{
                             {r0 * c0, r0 * s0, -hw},
                             {r0 * c0, r0 * s0, hw},
                             {r0 * c1, r0 * s1, hw},
                         }},
                         nx, ny, 0.0f, color);
        b3d_triangle_lit(&(b3d_tri_t) {{
                             {r0 * c0, r0 * s0, -hw},
                             {r0 * c1, r0 * s1, hw},
                             {r0 * c1, r0 * s1, -hw},
                         }},
                         nx, ny, 0.0f, color);
    }
}

/* Scene view angles (can be adjusted with keyboard in interactive mode) */
static float view_rotx = 20.0f;
static float view_roty = 30.0f;
static float view_rotz = 0.0f;

/* Gear animation angle (in degrees, matching original) */
static float angle_deg = 90.0f;

/* Render a frame
 * @pixels: output pixel buffer
 * @depth:  depth buffer
 * @width:  buffer width
 * @height: buffer height
 */
static void render_frame(uint32_t *pixels,
                         b3d_depth_t *depth,
                         int width,
                         int height)
{
    b3d_init(pixels, depth, width, height, 60.0f);
    b3d_set_camera(&(b3d_camera_t) {0.0f, 0.0f, -18.0f, 0.0f, 0.0f, 0.0f});
    b3d_clear();

    /* Set light direction: normalized (1, 1, 2) - matches original glxgears */
    b3d_set_light_direction(0.408248f, 0.408248f, 0.816497f);

    /* Apply view rotation */
    b3d_reset();
    b3d_rotate_x(DEG2RAD(view_rotx));
    b3d_rotate_y(DEG2RAD(view_roty));
    b3d_rotate_z(DEG2RAD(view_rotz));

    /* Gear 1 - large red gear with 20 teeth */
    b3d_push_matrix();
    b3d_rotate_z(DEG2RAD(angle_deg));
    b3d_translate(-3.0f, -2.0f, 0.0f);
    gear(1.0f, 4.0f, 1.0f, 20, 0.7f, 0xcc1900);
    b3d_pop_matrix();

    /* Gear 2 - medium green gear with 10 teeth
     * Positioned to mesh with gear 1 on the right
     * Rotates at 2:1 ratio (opposite direction) with phase offset */
    b3d_push_matrix();
    b3d_rotate_z(DEG2RAD(-2.0f * angle_deg - 9.0f));
    b3d_translate(3.1f, -2.0f, 0.0f);
    gear(0.5f, 2.0f, 2.0f, 10, 0.7f, 0x00cc33);
    b3d_pop_matrix();

    /* Gear 3 - small blue gear with 10 teeth
     * Positioned to mesh with gear 1 on top
     * Rotates at 2:1 ratio (opposite direction) with phase offset */
    b3d_push_matrix();
    b3d_rotate_z(DEG2RAD(-2.0f * angle_deg - 25.0f));
    b3d_translate(-3.1f, 4.2f, 0.0f);
    gear(1.3f, 2.0f, 0.5f, 10, 0.7f, 0x3333ff);
    b3d_pop_matrix();
}

int main(int argc, char **argv)
{
    int width = 800, height = 600;
    const char *snapshot = get_snapshot_path(argc, argv);

    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(pixels[0]));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(depth[0]));

    if (!pixels || !depth) {
        fprintf(stderr, "Failed to allocate buffers\n");
        free(pixels);
        free(depth);
        return 1;
    }

    if (snapshot) {
        angle_deg = 45.0f; /* Static snapshot angle */
        render_frame(pixels, depth, width, height);
        write_png(snapshot, pixels, width, height);
        free(pixels);
        free(depth);
        return 0;
    }

    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window *window =
        SDL_CreateWindow("B3D Gears", SDL_WINDOWPOS_CENTERED,
                         SDL_WINDOWPOS_CENTERED, width, height, 0);
    SDL_Renderer *renderer =
        SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC);
    SDL_Texture *texture =
        SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888,
                          SDL_TEXTUREACCESS_STREAMING, width, height);

    int quit = 0;
    Uint32 last_time = SDL_GetTicks();

    while (!quit) {
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = 1;
            } else if (event.type == SDL_KEYDOWN) {
                switch (event.key.keysym.scancode) {
                case SDL_SCANCODE_ESCAPE:
                    quit = 1;
                    break;
                case SDL_SCANCODE_UP:
                    view_rotx += 5.0f;
                    break;
                case SDL_SCANCODE_DOWN:
                    view_rotx -= 5.0f;
                    break;
                case SDL_SCANCODE_LEFT:
                    view_roty += 5.0f;
                    break;
                case SDL_SCANCODE_RIGHT:
                    view_roty -= 5.0f;
                    break;
                default:
                    break;
                }
            }
        }
        if (quit)
            break;

        /* Animate: ~70 degrees per second (matches glxgears) */
        Uint32 now = SDL_GetTicks();
        angle_deg += (now - last_time) * 0.07f;
        last_time = now;

        render_frame(pixels, depth, width, height);

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
