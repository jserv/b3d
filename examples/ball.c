/*
 * Amiga Bouncing Ball Demo
 *
 * Recreation of the iconic Amiga 1000 "Boing Ball" demo from 1984.
 * Features a red/white checkered geodesic sphere bouncing with gravity
 * simulation against a grid background.
 *
 * The original demo was created by Dale Luck and R.J. Mical in 1984
 * to showcase the Amiga's graphics capabilities at CES.
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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Grid configuration - original Amiga colors from Config.h */
#define GRID_COLS 15
#define GRID_ROWS 12
#define GRID_COLOR_BG 0xa0a0a0    /* RGB(160, 160, 160) gray */
#define GRID_LINE_COLOR 0x962896  /* RGB(150, 40, 150) purple/magenta */
#define GRID_SHADOW_BG 0x5a5a5a   /* RGB(90, 90, 90) darker gray */
#define GRID_SHADOW_LINE 0x500a50 /* RGB(80, 10, 80) darker purple */

/* Ball colors */
#define BALL_RED 0xff0000
#define BALL_WHITE 0xffffff

/* Physics constants */
#define BALL_RADIUS 0.8f
#define GRAVITY 9.8f
#define BOUNCE_DAMPING 0.95f
#define CEILING_DAMPING 0.8f
#define WALL_LEFT -1.8f
#define WALL_RIGHT 1.8f
#define FLOOR_Y -1.2f
#define CEILING_Y 1.5f

/* Ball spin constants */
#define MIN_BOUNCE_VELOCITY 2.2f
#define SPIN_ON_BOUNCE_SCALE 0.15f
#define SPIN_ON_BOUNCE_ADD 0.6f
#define SPIN_DECAY 0.995f
#define MIN_SPIN_SPEED 0.5f
#define SPIN_VX_FACTOR 1.2f

/* Lighting (shared with shadow projection) */
#define LIGHT_DIR_X 0.5f
#define LIGHT_DIR_Y 0.8f
#define LIGHT_DIR_Z -0.5f

/* Shadow scaling constants */
#define SHADOW_HEIGHT_FACTOR 0.35f
#define SHADOW_MAX_SCALE 1.35f
#define SHADOW_MIN_SCALE 0.6f

/* Ball rendering constants */
#define BALL_AXIAL_TILT 0.4f
#define AMBIENT_LIGHT 0.3f

/* Ball state */
typedef struct {
    float x, y, z;
    float vx, vy;
    float rotation;
    float spin_speed;
} ball_state_t;

/* UV sphere parameters for clean latitude/longitude-aligned tiles */
#define UV_SEGMENTS 30 /* Longitudinal divisions (div by BALL_STRIPES) */
#define UV_RINGS 18    /* Latitudinal divisions (divisible by BALL_BANDS) */
#define MAX_SPHERE_TRIS (UV_SEGMENTS * UV_RINGS * 2 + UV_SEGMENTS * 2)

/* Classic Amiga ball pattern: 6x6 tiles matching the original 1984 demo */
#define BALL_STRIPES 6 /* Longitudinal tile divisions */
#define BALL_BANDS 6   /* Latitudinal tile divisions */

static b3d_tri_t sphere_tris[MAX_SPHERE_TRIS];
static float sphere_normals[MAX_SPHERE_TRIS][3];
static int sphere_colors[MAX_SPHERE_TRIS];
static int sphere_tri_count = 0;

static float clampf(float v, float lo, float hi)
{
    if (v < lo)
        return lo;
    if (v > hi)
        return hi;
    return v;
}

/* Determine tile color based on longitude and latitude indices
 * @lon_idx: longitude segment index (0 to BALL_STRIPES-1)
 * @lat_idx: latitude band index (0 to BALL_BANDS-1)
 */
static uint32_t get_tile_color(int lon_idx, int lat_idx)
{
    /* Checkerboard pattern: alternate based on sum of indices */
    if ((lon_idx + lat_idx) % 2 == 0)
        return BALL_RED;
    return BALL_WHITE;
}

/* Add a triangle to the sphere mesh
 * @v0, @v1, @v2: triangle vertices
 * @color: face color
 */
static void add_sphere_tri(const float v0[3],
                           const float v1[3],
                           const float v2[3],
                           uint32_t color)
{
    if (sphere_tri_count >= MAX_SPHERE_TRIS)
        return;

    b3d_tri_t *tri = &sphere_tris[sphere_tri_count];
    tri->v[0].x = v0[0];
    tri->v[0].y = v0[1];
    tri->v[0].z = v0[2];
    tri->v[1].x = v1[0];
    tri->v[1].y = v1[1];
    tri->v[1].z = v1[2];
    tri->v[2].x = v2[0];
    tri->v[2].y = v2[1];
    tri->v[2].z = v2[2];

    /* Calculate face normal from centroid (for sphere, normal = position) */
    float nx = (v0[0] + v1[0] + v2[0]) / 3.0f;
    float ny = (v0[1] + v1[1] + v2[1]) / 3.0f;
    float nz = (v0[2] + v1[2] + v2[2]) / 3.0f;
    float len = b3d_sqrtf(nx * nx + ny * ny + nz * nz);
    if (len > 0.0001f) {
        nx /= len;
        ny /= len;
        nz /= len;
    }

    sphere_normals[sphere_tri_count][0] = nx;
    sphere_normals[sphere_tri_count][1] = ny;
    sphere_normals[sphere_tri_count][2] = nz;
    sphere_colors[sphere_tri_count] = color;
    sphere_tri_count++;
}

/* Generate UV sphere with latitude/longitude-aligned tiles
 * Creates a sphere where each tile boundary aligns with the color pattern,
 * resulting in clean checkerboard squares.
 */
static void generate_sphere(void)
{
    sphere_tri_count = 0;

    /* Calculate how many geometry segments per color tile */
    int segs_per_stripe = UV_SEGMENTS / BALL_STRIPES;
    int rings_per_band = UV_RINGS / BALL_BANDS;

    for (int lat = 0; lat < UV_RINGS; lat++) {
        /* Latitude angles (phi): from top pole (0) to bottom pole (π) */
        float phi0 = (float) M_PI * lat / UV_RINGS;
        float phi1 = (float) M_PI * (lat + 1) / UV_RINGS;

        float sin_phi0, cos_phi0, sin_phi1, cos_phi1;
        b3d_sincosf(phi0, &sin_phi0, &cos_phi0);
        b3d_sincosf(phi1, &sin_phi1, &cos_phi1);

        /* Determine which color band this latitude ring belongs to */
        int lat_band = lat / rings_per_band;
        if (lat_band >= BALL_BANDS)
            lat_band = BALL_BANDS - 1;

        for (int lon = 0; lon < UV_SEGMENTS; lon++) {
            /* Longitude angles (theta): full circle with rotation */
            float theta0 = 2.0f * (float) M_PI * lon / UV_SEGMENTS;
            float theta1 = 2.0f * (float) M_PI * (lon + 1) / UV_SEGMENTS;

            float sin_theta0, cos_theta0, sin_theta1, cos_theta1;
            b3d_sincosf(theta0, &sin_theta0, &cos_theta0);
            b3d_sincosf(theta1, &sin_theta1, &cos_theta1);

            /* Four corners of the quad */
            float v00[3] = {sin_phi0 * cos_theta0, cos_phi0,
                            sin_phi0 * sin_theta0};
            float v01[3] = {sin_phi0 * cos_theta1, cos_phi0,
                            sin_phi0 * sin_theta1};
            float v10[3] = {sin_phi1 * cos_theta0, cos_phi1,
                            sin_phi1 * sin_theta0};
            float v11[3] = {sin_phi1 * cos_theta1, cos_phi1,
                            sin_phi1 * sin_theta1};

            /* Determine which color stripe this longitude segment belongs to */
            int lon_stripe = lon / segs_per_stripe;
            if (lon_stripe >= BALL_STRIPES)
                lon_stripe = BALL_STRIPES - 1;

            uint32_t color = get_tile_color(lon_stripe, lat_band);

            /* Top cap (first ring connects to pole) */
            if (lat == 0) {
                float pole[3] = {0.0f, 1.0f, 0.0f};
                add_sphere_tri(pole, v01, v00, color);
            }
            /* Bottom cap (last ring connects to pole) */
            else if (lat == UV_RINGS - 1) {
                float pole[3] = {0.0f, -1.0f, 0.0f};
                add_sphere_tri(v00, v01, pole, color);
            }
            /* Regular quad (two triangles) */
            else {
                add_sphere_tri(v00, v01, v11, color);
                add_sphere_tri(v00, v11, v10, color);
            }
        }
    }
}

/* Draw classic Amiga-style flat 2D grid background
 * @pixels: pixel buffer
 * @width: buffer width
 * @height: buffer height
 *
 * The original Amiga demo used a simple flat rectangular grid
 * with purple lines on a gray background.
 */
static void draw_grid_background(uint32_t *pixels, int width, int height)
{
    /* Fill background with solid gray color */
    for (int i = 0; i < width * height; i++)
        pixels[i] = GRID_COLOR_BG;

    /* Calculate cell size to fill the screen */
    int cell_w = width / GRID_COLS;
    int cell_h = height / GRID_ROWS;

    /* Draw vertical lines */
    for (int col = 0; col <= GRID_COLS; col++) {
        int x = col * cell_w;
        if (x >= width)
            x = width - 1;
        for (int y = 0; y < height; y++)
            pixels[y * width + x] = GRID_LINE_COLOR;
    }

    /* Draw horizontal lines */
    for (int row = 0; row <= GRID_ROWS; row++) {
        int y = row * cell_h;
        if (y >= height)
            y = height - 1;
        for (int x = 0; x < width; x++)
            pixels[y * width + x] = GRID_LINE_COLOR;
    }
}

/* Update ball physics
 * @ball: ball state to update
 * @dt: time delta in seconds
 */
static void update_ball_physics(ball_state_t *ball, float dt)
{
    float ground = FLOOR_Y + BALL_RADIUS;

    /* Apply gravity */
    ball->vy -= GRAVITY * dt;

    /* Update position */
    ball->x += ball->vx * dt;
    ball->y += ball->vy * dt;

    /* Bounce off floor */
    if (ball->y < ground) {
        ball->y = ground;
        ball->vy = -ball->vy * BOUNCE_DAMPING;
        /* Keep the ball lively like the original demo */
        if (b3d_fabsf(ball->vy) < MIN_BOUNCE_VELOCITY)
            ball->vy = MIN_BOUNCE_VELOCITY;
        ball->spin_speed +=
            b3d_fabsf(ball->vx) * SPIN_ON_BOUNCE_SCALE + SPIN_ON_BOUNCE_ADD;
    }

    /* Bounce off ceiling (optional, for dramatic effect) */
    if (ball->y > CEILING_Y) {
        ball->y = CEILING_Y;
        ball->vy = -ball->vy * CEILING_DAMPING;
    }

    /* Bounce off walls */
    if (ball->x < WALL_LEFT) {
        ball->x = WALL_LEFT;
        ball->vx = -ball->vx;
    }
    if (ball->x > WALL_RIGHT) {
        ball->x = WALL_RIGHT;
        ball->vx = -ball->vx;
    }

    /* Spin decays slowly but stays responsive to movement */
    ball->spin_speed *= SPIN_DECAY;
    if (ball->spin_speed < MIN_SPIN_SPEED)
        ball->spin_speed = MIN_SPIN_SPEED;

    ball->rotation += (ball->vx * SPIN_VX_FACTOR + ball->spin_speed) * dt;
}

/* Draw ball shadow on grid
 * @pixels: pixel buffer
 * @width: buffer width
 * @height: buffer height
 * @ball_x: ball x position in world coords
 * @ball_y: ball y position in world coords
 * @ball_z: ball z position in world coords
 *
 * The shadow is drawn by darkening the grid cells it covers,
 * preserving the grid line pattern like in the original demo.
 * Shadow position is calculated by projecting from light direction.
 */
static void draw_shadow(uint32_t *pixels,
                        int width,
                        int height,
                        float ball_x,
                        float ball_y,
                        float ball_z)
{
    /* Calculate shadow position based on light direction.
     * Trace ray from ball center in opposite direction of light until floor.
     * Ray: ball_pos + t * (-light_dir), solve for y = FLOOR_Y
     */
    float height_above_floor = ball_y - FLOOR_Y;
    if (height_above_floor <= 0.0f || LIGHT_DIR_Y <= 0.0f)
        return;

    float shadow_x = ball_x - height_above_floor * (LIGHT_DIR_X / LIGHT_DIR_Y);
    float shadow_z = ball_z - height_above_floor * (LIGHT_DIR_Z / LIGHT_DIR_Y);

    int shadow_sx, shadow_sy;
    if (!b3d_to_screen(shadow_x, FLOOR_Y, shadow_z, &shadow_sx, &shadow_sy))
        return;

    /* Clamp shadow to stay on screen */
    if (shadow_sy > height - 2)
        shadow_sy = height - 2;
    if (shadow_sy < 0)
        shadow_sy = 0;

    /* Scale shadow size inversely with ball height above the floor */
    float contact_scale = clampf(
        SHADOW_MAX_SCALE / (SHADOW_HEIGHT_FACTOR * height_above_floor + 1.0f),
        SHADOW_MIN_SCALE, SHADOW_MAX_SCALE);

    /* Derive ellipse axes from projecting a circle on the floor plane.
     * Project points around the shadow center to get proper perspective.
     */
    int base_rx = 0, base_ry = 0;
    int proj_x_sx, proj_x_sy, proj_z_sx, proj_z_sy;

    if (b3d_to_screen(shadow_x + BALL_RADIUS, FLOOR_Y, shadow_z, &proj_x_sx,
                      &proj_x_sy)) {
        base_rx = abs(proj_x_sx - shadow_sx);
    }
    if (b3d_to_screen(shadow_x, FLOOR_Y, shadow_z + BALL_RADIUS, &proj_z_sx,
                      &proj_z_sy)) {
        /* Use the screen-space vertical offset from z extent to shape ellipse
         */
        base_ry = abs(proj_z_sy - shadow_sy);
    }

    /* Fallback sizes if projection degenerates */
    if (base_rx < 10)
        base_rx = 60;
    if (base_ry < 6)
        base_ry = 20;

    int radius_x = (int) (base_rx * contact_scale);
    int radius_y = (int) (base_ry * contact_scale);
    if (radius_x < 20)
        radius_x = 20;
    if (radius_y < 8)
        radius_y = 8;

    /* Cell dimensions for grid line detection */
    int cell_w = width / GRID_COLS;
    int cell_h = height / GRID_ROWS;

    /* Draw elliptical shadow, darkening grid cells while preserving lines */
    for (int dy = -radius_y; dy <= radius_y; dy++) {
        float y_frac = (float) dy / (float) radius_y;
        float dx_max = b3d_sqrtf(1.0f - y_frac * y_frac) * (float) radius_x;
        for (int dx = -(int) dx_max; dx <= (int) dx_max; dx++) {
            int px = shadow_sx + dx;
            int py = shadow_sy + dy;
            if (px >= 0 && px < width && py >= 0 && py < height) {
                /* Check if this pixel is on a grid line */
                int on_vline = (px % cell_w == 0);
                int on_hline = (py % cell_h == 0);
                if (on_vline || on_hline)
                    pixels[py * width + px] = GRID_SHADOW_LINE;
                else
                    pixels[py * width + px] = GRID_SHADOW_BG;
            }
        }
    }
}

/* Render scene
 * @pixels: pixel buffer
 * @depth: depth buffer
 * @width: buffer width
 * @height: buffer height
 * @ball: ball state
 */
static void render_scene(uint32_t *pixels,
                         b3d_depth_t *depth,
                         int width,
                         int height,
                         ball_state_t *ball)
{
    /* Initialize B3D */
    b3d_init(pixels, depth, width, height, 60);
    b3d_set_camera(&(b3d_camera_t) {0, 0, -4, 0, 0, 0});

    /* Clear depth buffer and draw grid background */
    b3d_clear();
    draw_grid_background(pixels, width, height);

    /* Draw shadow under the ball */
    draw_shadow(pixels, width, height, ball->x, ball->y, ball->z);

    /* Set up lighting from upper-left-front */
    b3d_set_light_direction(LIGHT_DIR_X, LIGHT_DIR_Y, LIGHT_DIR_Z);
    b3d_set_ambient(AMBIENT_LIGHT);

    /* Transform and render ball */
    b3d_reset();
    b3d_translate(ball->x, ball->y, ball->z);
    b3d_scale(BALL_RADIUS, BALL_RADIUS, BALL_RADIUS);
    b3d_rotate_y(ball->rotation);
    /* Tilt the ball like in the original demo */
    b3d_rotate_z(BALL_AXIAL_TILT);

    for (int i = 0; i < sphere_tri_count; i++) {
        b3d_triangle_lit(&sphere_tris[i], sphere_normals[i][0],
                         sphere_normals[i][1], sphere_normals[i][2],
                         sphere_colors[i]);
    }
}

int main(int argc, char **argv)
{
    int width = 800, height = 600;
    const char *snapshot = get_snapshot_path(argc, argv);

    uint32_t *pixels = malloc(width * height * sizeof(uint32_t));
    if (!pixels) {
        fprintf(stderr, "Failed to allocate pixel buffer\n");
        return 1;
    }
    b3d_depth_t *depth = malloc(width * height * sizeof(b3d_depth_t));
    if (!depth) {
        fprintf(stderr, "Failed to allocate depth buffer\n");
        free(pixels);
        return 1;
    }
    generate_sphere();

    /* Initialize ball state */
    ball_state_t ball = {
        .x = 0.0f,
        .y = 1.0f,
        .z = 0.0f,
        .vx = 1.2f,
        .vy = 0.0f,
        .rotation = 0.0f,
        .spin_speed = 2.4f,
    };

    if (snapshot) {
        /* Render single frame for snapshot */
        ball.y = 0.8f;
        ball.rotation = 0.5f;
        render_scene(pixels, depth, width, height, &ball);
        write_png(snapshot, pixels, width, height);
        free(pixels);
        free(depth);
        return 0;
    }

    /* SDL2 setup */
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window *window =
        SDL_CreateWindow("Amiga Boing Ball", SDL_WINDOWPOS_CENTERED,
                         SDL_WINDOWPOS_CENTERED, width, height, 0);
    SDL_Renderer *renderer =
        SDL_CreateRenderer(window, -1, SDL_RENDERER_PRESENTVSYNC);
    SDL_Texture *texture =
        SDL_CreateTexture(renderer, SDL_PIXELFORMAT_ARGB8888,
                          SDL_TEXTUREACCESS_STREAMING, width, height);

    uint32_t last_time = SDL_GetTicks();

    int quit = 0;
    while (!quit) {
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

        /* Calculate delta time */
        uint32_t current_time = SDL_GetTicks();
        float dt = (current_time - last_time) / 1000.0f;
        last_time = current_time;

        /* Cap delta time to avoid physics explosion */
        if (dt > 0.1f)
            dt = 0.1f;

        /* Update physics */
        update_ball_physics(&ball, dt);

        /* Render */
        render_scene(pixels, depth, width, height, &ball);

        /* Display */
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
