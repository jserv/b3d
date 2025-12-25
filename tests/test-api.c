/*
 * B3D API validation test using ASCII rendering
 * Validates core B3D APIs by rendering a rotating cube and comparing ASCII
 * output
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../include/b3d.h"

/* ANSI color codes for terminal output */
#define ANSI_GREEN "\033[32m"
#define ANSI_RED "\033[31m"
#define ANSI_BOLD "\033[1m"
#define ANSI_RESET "\033[0m"

static int tests_run = 0;
static int tests_passed = 0;
static int section_passed = 0;
static int section_total = 0;

#define TEST(name) static int test_##name(void)

#define RUN_TEST(name)                        \
    do {                                      \
        tests_run++, section_total++;         \
        if (test_##name())                    \
            tests_passed++, section_passed++; \
    } while (0)

#define SECTION_BEGIN(name)                      \
    do {                                         \
        section_passed = 0, section_total = 0;   \
        printf("  %-44s", name), fflush(stdout); \
    } while (0)

#define SECTION_END()                                            \
    do {                                                         \
        if (section_passed == section_total)                     \
            printf("[ " ANSI_GREEN "OK" ANSI_RESET " ]\n");      \
        else                                                     \
            printf("[ " ANSI_RED "FAIL %d/%d" ANSI_RESET " ]\n", \
                   section_passed, section_total);               \
    } while (0)

#define ASSERT(cond)  \
    do {              \
        if (!(cond))  \
            return 0; \
    } while (0)

static const char *palette = " .:-=+*#%@";

static char color_to_char(uint32_t c)
{
    size_t len = strlen(palette);
    if (len == 0)
        return ' ';

    unsigned r = (c >> 16) & 0xff;
    unsigned g = (c >> 8) & 0xff;
    unsigned b = (c >> 0) & 0xff;
    float lum = (0.299f * r + 0.587f * g + 0.114f * b) / 255.0f;
    int idx = (int) (lum * (float) (len - 1) + 0.5f);
    if (idx < 0)
        idx = 0;
    if ((size_t) idx >= len)
        idx = (int) len - 1;
    return palette[idx];
}

/* Test basic initialization and buffer management */
TEST(api_init)
{
    const int width = 32, height = 32;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    ASSERT(pixels);
    ASSERT(depth);

    /* Test b3d_init and b3d_clear */
    bool init_ok = b3d_init(pixels, depth, width, height, 65.0f);
    ASSERT(init_ok);
    b3d_clear();
    /* Verify buffer was cleared (all pixels should be 0 after clear) */
    const size_t pixel_count = (size_t) width * (size_t) height;
    size_t non_zero = 0;
    for (size_t i = 0; i < pixel_count; i++) {
        if (pixels[i] != 0)
            non_zero++;
    }
    /* After clear, pixels should be 0, depth should be at far plane */
    ASSERT(non_zero == 0);

    free(pixels);
    free(depth);
    return 1;
}

/* Test transformation APIs */
TEST(api_transform)
{
    const int width = 32, height = 32;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    /* Check allocation success */
    if (!pixels || !depth) {
        free(pixels);
        free(depth);
        return 0;
    }

    b3d_init(pixels, depth, width, height, 65.0f);
    b3d_clear();

    /* Test b3d_reset (the b3d_init already calls b3d_reset at the end) */
    float initial_matrix[16];
    b3d_get_model_matrix(initial_matrix);
    /* Check if it's identity matrix (diagonal should be 1, others 0) */
    ASSERT(fabsf(initial_matrix[0] - 1.0f) < 0.01f);  /* m[0][0] */
    ASSERT(fabsf(initial_matrix[5] - 1.0f) < 0.01f);  /* m[1][1] */
    ASSERT(fabsf(initial_matrix[10] - 1.0f) < 0.01f); /* m[2][2] */
    ASSERT(fabsf(initial_matrix[15] - 1.0f) < 0.01f); /* m[3][3] */
    /* Check that off-diagonal elements are approximately 0 */
    ASSERT(fabsf(initial_matrix[1]) < 0.01f);  /* m[0][1] */
    ASSERT(fabsf(initial_matrix[2]) < 0.01f);  /* m[0][2] */
    ASSERT(fabsf(initial_matrix[3]) < 0.01f);  /* m[0][3] */
    ASSERT(fabsf(initial_matrix[4]) < 0.01f);  /* m[1][0] */
    ASSERT(fabsf(initial_matrix[6]) < 0.01f);  /* m[1][2] */
    ASSERT(fabsf(initial_matrix[7]) < 0.01f);  /* m[1][3] */
    ASSERT(fabsf(initial_matrix[8]) < 0.01f);  /* m[2][0] */
    ASSERT(fabsf(initial_matrix[9]) < 0.01f);  /* m[2][1] */
    ASSERT(fabsf(initial_matrix[11]) < 0.01f); /* m[2][3] */
    ASSERT(fabsf(initial_matrix[12]) < 0.01f); /* m[3][0] */
    ASSERT(fabsf(initial_matrix[13]) < 0.01f); /* m[3][1] */
    ASSERT(fabsf(initial_matrix[14]) < 0.01f); /* m[3][2] */

    /* Test b3d_translate */
    b3d_translate(1.0f, 2.0f, 3.0f);
    float matrix[16];
    b3d_get_model_matrix(matrix);
    /* Translation should affect the last column */
    ASSERT(fabsf(matrix[12] - 1.0f) < 0.01f); /* tx */
    ASSERT(fabsf(matrix[13] - 2.0f) < 0.01f); /* ty */
    ASSERT(fabsf(matrix[14] - 3.0f) < 0.01f); /* tz */

    /* Test b3d_rotate_x - use smaller angle for more predictable results */
    b3d_reset();
    b3d_rotate_x(0.1f); /* Small angle */
    float rot_x_matrix[16];
    b3d_get_model_matrix(rot_x_matrix);
    /* Based on debug output, after b3d_rotate_x(0.1f): */
    /* [5]~0.995004 (cos), [6]~0.099833 (sin), [9]~-0.099833 (-sin), */
    /* [10]~0.995004 (cos) */
    ASSERT(fabsf(rot_x_matrix[5] - 0.995004f) <
           0.01f); /* m[1][1] should be ~cos(0.1) */
    ASSERT(fabsf(rot_x_matrix[6] - 0.099833f) <
           0.01f); /* m[1][2] should be ~sin(0.1) */
    ASSERT(fabsf(rot_x_matrix[9] + 0.099833f) <
           0.01f); /* m[2][1] should be ~-sin(0.1) */
    ASSERT(fabsf(rot_x_matrix[10] - 0.995004f) <
           0.01f); /* m[2][2] should be ~cos(0.1) */

    /* Test b3d_rotate_y - use smaller angle for more predictable results */
    b3d_reset();
    b3d_rotate_y(0.1f); /* Small angle */
    float rot_y_matrix[16];
    b3d_get_model_matrix(rot_y_matrix);
    /* Based on debug output, after b3d_rotate_y(0.1f): */
    /* [0]~0.995004 (cos), [2]~0.099833 (sin), [8]~-0.099833 (-sin), */
    /* [10]~0.995004 (cos) */
    ASSERT(fabsf(rot_y_matrix[0] - 0.995004f) <
           0.01f); /* m[0][0] should be ~cos(0.1) */
    ASSERT(fabsf(rot_y_matrix[2] - 0.099833f) <
           0.01f); /* m[0][2] should be ~sin(0.1) */
    ASSERT(fabsf(rot_y_matrix[8] + 0.099833f) <
           0.01f); /* m[2][0] should be ~-sin(0.1) */
    ASSERT(fabsf(rot_y_matrix[10] - 0.995004f) <
           0.01f); /* m[2][2] should be ~cos(0.1) */

    /* Test b3d_rotate_z - use smaller angle for more predictable results */
    b3d_reset();
    b3d_rotate_z(0.1f); /* Small angle */
    float rot_z_matrix[16];
    b3d_get_model_matrix(rot_z_matrix);
    /* Based on debug output, after b3d_rotate_z(0.1f): */
    /* [0]~0.995004 (cos), [1]~0.099833 (sin), [4]~-0.099833 (-sin), */
    /* [5]~0.995004 (cos) */
    ASSERT(fabsf(rot_z_matrix[0] - 0.995004f) <
           0.01f); /* m[0][0] should be ~cos(0.1) */
    ASSERT(fabsf(rot_z_matrix[1] - 0.099833f) <
           0.01f); /* m[0][1] should be ~sin(0.1) */
    ASSERT(fabsf(rot_z_matrix[4] + 0.099833f) <
           0.01f); /* m[1][0] should be ~-sin(0.1) */
    ASSERT(fabsf(rot_z_matrix[5] - 0.995004f) <
           0.01f); /* m[1][1] should be ~cos(0.1) */

    /* Test b3d_scale */
    b3d_reset();
    b3d_scale(2.0f, 3.0f, 4.0f);
    float scale_matrix[16];
    b3d_get_model_matrix(scale_matrix);
    ASSERT(fabsf(scale_matrix[0] - 2.0f) < 0.01f);  /* m[0][0] should be 2 */
    ASSERT(fabsf(scale_matrix[5] - 3.0f) < 0.01f);  /* m[1][1] should be 3 */
    ASSERT(fabsf(scale_matrix[10] - 4.0f) < 0.01f); /* m[2][2] should be 4 */

    free(pixels);
    free(depth);
    return 1;
}

/* Test camera APIs */
TEST(api_camera)
{
    const int width = 32, height = 32;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    b3d_init(pixels, depth, width, height, 65.0f);
    b3d_clear();

    /* Test b3d_set_camera and b3d_get_camera roundtrip */
    b3d_camera_t cam_in = {1.0f, 2.0f, 3.0f, 0.1f, 0.2f, 0.3f};
    b3d_set_camera(&cam_in);

    b3d_camera_t cam_out;
    b3d_get_camera(&cam_out);

    /* Verify all fields match exactly */
    ASSERT(fabsf(cam_out.x - cam_in.x) < 0.0001f);
    ASSERT(fabsf(cam_out.y - cam_in.y) < 0.0001f);
    ASSERT(fabsf(cam_out.z - cam_in.z) < 0.0001f);
    ASSERT(fabsf(cam_out.yaw - cam_in.yaw) < 0.0001f);
    ASSERT(fabsf(cam_out.pitch - cam_in.pitch) < 0.0001f);
    ASSERT(fabsf(cam_out.roll - cam_in.roll) < 0.0001f);

    /* Test b3d_set_fov and b3d_get_fov roundtrip */
    b3d_set_fov(90.0f);
    float fov = b3d_get_fov();
    ASSERT(fabsf(fov - 90.0f) < 0.0001f);

    /* Test b3d_look_at (note: invalidates camera orientation) */
    b3d_look_at(5.0f, 6.0f, 7.0f);
    /* Position should remain unchanged, orientation becomes stale */
    b3d_get_camera(&cam_out);
    ASSERT(fabsf(cam_out.x - cam_in.x) < 0.0001f);
    ASSERT(fabsf(cam_out.y - cam_in.y) < 0.0001f);
    ASSERT(fabsf(cam_out.z - cam_in.z) < 0.0001f);

    /* Test b3d_get_view_matrix returns valid data */
    float view[16];
    b3d_get_view_matrix(view);
    /* View matrix should be non-zero (not identity after look_at) */
    int non_zero = 0;
    for (int i = 0; i < 16; i++) {
        if (fabsf(view[i]) > 0.0001f)
            non_zero++;
    }
    ASSERT(non_zero > 0);

    /* Test b3d_get_proj_matrix returns valid data */
    float proj[16];
    b3d_get_proj_matrix(proj);
    /* Projection matrix should have non-zero elements */
    non_zero = 0;
    for (int i = 0; i < 16; i++) {
        if (fabsf(proj[i]) > 0.0001f)
            non_zero++;
    }
    ASSERT(non_zero > 0);

    free(pixels);
    free(depth);
    return 1;
}

/* Test state query APIs */
TEST(api_state_queries)
{
    const int width = 64, height = 48;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    ASSERT(pixels);
    ASSERT(depth);

    /* Test b3d_is_initialized before init */
    /* Note: After a failed init, state is cleared, so this depends on prior
     * state */

    /* Initialize and verify state */
    bool init_ok = b3d_init(pixels, depth, width, height, 75.0f);
    ASSERT(init_ok);

    /* Test b3d_is_initialized after successful init */
    ASSERT(b3d_is_initialized() == true);

    /* Test b3d_get_width and b3d_get_height */
    ASSERT(b3d_get_width() == width);
    ASSERT(b3d_get_height() == height);

    /* Test b3d_get_fov returns init value */
    float fov = b3d_get_fov();
    ASSERT(fabsf(fov - 75.0f) < 0.0001f);

    /* Test after re-init with different dimensions */
    const int width2 = 128, height2 = 96;
    uint32_t *pixels2 =
        malloc((size_t) width2 * (size_t) height2 * sizeof(uint32_t));
    b3d_depth_t *depth2 =
        malloc((size_t) width2 * (size_t) height2 * sizeof(b3d_depth_t));

    ASSERT(pixels2);
    ASSERT(depth2);

    init_ok = b3d_init(pixels2, depth2, width2, height2, 60.0f);
    ASSERT(init_ok);
    ASSERT(b3d_get_width() == width2);
    ASSERT(b3d_get_height() == height2);
    ASSERT(fabsf(b3d_get_fov() - 60.0f) < 0.0001f);

    free(pixels);
    free(depth);
    free(pixels2);
    free(depth2);
    return 1;
}

/* Test matrix stack APIs */
TEST(api_matrix_stack)
{
    const int width = 32, height = 32;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    b3d_init(pixels, depth, width, height, 65.0f);
    b3d_clear();

    /* Test b3d_push_matrix and b3d_pop_matrix */
    b3d_reset();
    b3d_translate(1.0f, 0.0f, 0.0f);
    float matrix_before[16];
    b3d_get_model_matrix(matrix_before);
    ASSERT(matrix_before[12] == 1.0f); /* tx should be 1 */

    /* Push current matrix */
    bool push_result = b3d_push_matrix();
    ASSERT(push_result == true); /* Should succeed */

    /* Modify matrix */
    b3d_translate(2.0f, 0.0f, 0.0f);
    float matrix_after[16];
    b3d_get_model_matrix(matrix_after);
    ASSERT(matrix_after[12] == 3.0f); /* tx should be 3 (1+2) */

    /* Pop matrix - should restore previous state */
    bool pop_result = b3d_pop_matrix();
    ASSERT(pop_result == true); /* Should succeed */

    float matrix_restored[16];
    b3d_get_model_matrix(matrix_restored);
    ASSERT(fabsf(matrix_restored[12] - 1.0f) <
           0.01f); /* tx should be back to 1 */

    /* Test stack overflow protection */
    for (int i = 0; i < B3D_MATRIX_STACK_SIZE + 5; i++) {
        if (b3d_push_matrix() == false)
            break; /* Stack is full, this is expected */
    }

    /* Test stack underflow protection */
    for (int i = 0; i < B3D_MATRIX_STACK_SIZE + 5; i++) {
        if (b3d_pop_matrix() == false)
            break; /* Stack is empty, this is expected */
    }

    free(pixels);
    free(depth);
    return 1;
}

/* Test rendering API with ASCII output comparison */
TEST(api_render_ascii)
{
    const int width = 32, height = 16;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    ASSERT(pixels);
    ASSERT(depth);

    /* Render a simple scene */
    bool init_ok = b3d_init(pixels, depth, width, height, 65.0f);
    ASSERT(init_ok);
    b3d_set_camera(&(b3d_camera_t) {0.0f, 0.0f, -2.3f, 0.0f, 0.0f, 0.0f});
    b3d_clear();

    b3d_reset();
    b3d_rotate_y(0.5f);
    b3d_rotate_x(0.3f);

    /* Render a simple cube */
    b3d_triangle(
        &(b3d_tri_t) {
            {{-0.5f, -0.5f, -0.5f}, {-0.5f, 0.5f, -0.5f}, {0.5f, 0.5f, -0.5f}}},
        0xfcd0a1);
    b3d_triangle(
        &(b3d_tri_t) {
            {{-0.5f, -0.5f, -0.5f}, {0.5f, 0.5f, -0.5f}, {0.5f, -0.5f, -0.5f}}},
        0xb1b695);
    b3d_triangle(
        &(b3d_tri_t) {
            {{0.5f, -0.5f, -0.5f}, {0.5f, 0.5f, -0.5f}, {0.5f, 0.5f, 0.5f}}},
        0x53917e);
    b3d_triangle(
        &(b3d_tri_t) {
            {{0.5f, -0.5f, -0.5f}, {0.5f, 0.5f, 0.5f}, {0.5f, -0.5f, 0.5f}}},
        0x63535b);
    b3d_triangle(
        &(b3d_tri_t) {
            {{0.5f, -0.5f, 0.5f}, {0.5f, 0.5f, 0.5f}, {-0.5f, 0.5f, 0.5f}}},
        0x6d1a36);
    b3d_triangle(
        &(b3d_tri_t) {
            {{0.5f, -0.5f, 0.5f}, {-0.5f, 0.5f, 0.5f}, {-0.5f, -0.5f, 0.5f}}},
        0xd4e09b);
    b3d_triangle(
        &(b3d_tri_t) {
            {{-0.5f, -0.5f, 0.5f}, {-0.5f, 0.5f, 0.5f}, {-0.5f, 0.5f, -0.5f}}},
        0xf6f4d2);
    b3d_triangle(&(b3d_tri_t) {{{-0.5f, -0.5f, 0.5f},
                                {-0.5f, 0.5f, -0.5f},
                                {-0.5f, -0.5f, -0.5f}}},
                 0xcbdfbd);
    b3d_triangle(
        &(b3d_tri_t) {
            {{-0.5f, 0.5f, -0.5f}, {-0.5f, 0.5f, 0.5f}, {0.5f, 0.5f, 0.5f}}},
        0xf19c79);
    b3d_triangle(
        &(b3d_tri_t) {
            {{-0.5f, 0.5f, -0.5f}, {0.5f, 0.5f, 0.5f}, {0.5f, 0.5f, -0.5f}}},
        0xa44a3f);
    b3d_triangle(
        &(b3d_tri_t) {
            {{0.5f, -0.5f, 0.5f}, {-0.5f, -0.5f, 0.5f}, {-0.5f, -0.5f, -0.5f}}},
        0x5465ff);
    b3d_triangle(
        &(b3d_tri_t) {
            {{0.5f, -0.5f, 0.5f}, {-0.5f, -0.5f, -0.5f}, {0.5f, -0.5f, -0.5f}}},
        0x788bff);

    /* Convert to ASCII and verify we have some non-empty output */
    size_t non_empty_chars = 0;
    for (size_t y = 0; y < (size_t) height; ++y) {
        for (size_t x = 0; x < (size_t) width; ++x) {
            uint32_t c = pixels[(size_t) y * (size_t) width + (size_t) x];
            char ch = c ? color_to_char(c) : ' ';
            if (ch != ' ')
                non_empty_chars++;
        }
    }

    /* We should have rendered some pixels (non-space characters) */
    ASSERT(non_empty_chars > 0);

    free(pixels);
    free(depth);
    return 1;
}

/* Test b3d_to_screen API */
TEST(api_to_screen)
{
    const int width = 64, height = 64;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    b3d_init(pixels, depth, width, height, 65.0f);
    b3d_clear();

    /* Test b3d_to_screen - should project world coordinates to screen */
    int sx, sy;
    /* Set up camera first */
    b3d_set_camera(
        &(b3d_camera_t) {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}); /* At origin */
    bool result = b3d_to_screen(0.0f, 0.0f, 1.0f, &sx,
                                &sy); /* Point in front of camera */
    /* This should succeed (return true) since (0,0,1) is in front of camera */
    ASSERT(result == true);
    /* The point should be projected to screen bounds */
    ASSERT(sx >= 0 && sx < width);
    ASSERT(sy >= 0 && sy < height);

    /* Test with NULL pointers - should return false */
    result = b3d_to_screen(0.0f, 0.0f, 0.0f, NULL, &sy);
    ASSERT(result == false);
    result = b3d_to_screen(0.0f, 0.0f, 0.0f, &sx, NULL);
    ASSERT(result == false);

    free(pixels);
    free(depth);
    return 1;
}

/* Test b3d_get_clip_drop_count API */
TEST(api_clip_drop_count)
{
    const int width = 32, height = 32;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    b3d_init(pixels, depth, width, height, 65.0f);
    b3d_clear();

    /* Get initial count */
    size_t initial_count = b3d_get_clip_drop_count();

    /* Render a valid triangle */
    b3d_triangle(
        &(b3d_tri_t) {
            {{0.0f, 0.0f, 0.5f}, {0.5f, 0.5f, 0.5f}, {-0.5f, 0.5f, 0.5f}}},
        0xffffff);

    /* Count should still be the same (no clipping for valid triangle) */
    size_t after_count = b3d_get_clip_drop_count();
    ASSERT(initial_count == after_count);

    /* Clear should reset count */
    b3d_clear();
    size_t count = b3d_get_clip_drop_count();
    ASSERT(count == 0);

    free(pixels);
    free(depth);
    return 1;
}

/* Test buffer size calculation API */
TEST(api_buffer_size)
{
    size_t size;

    /* Test normal cases */
    size = b3d_buffer_size(10, 10, sizeof(uint32_t));
    ASSERT(size == 10 * 10 * sizeof(uint32_t));
    size = b3d_buffer_size(100, 50, sizeof(b3d_depth_t));
    ASSERT(size == 100 * 50 * sizeof(b3d_depth_t));

    /* Test edge cases */
    size = b3d_buffer_size(0, 10, sizeof(uint32_t));
    ASSERT(size == 0);
    size = b3d_buffer_size(10, 0, sizeof(uint32_t));
    ASSERT(size == 0);
    size = b3d_buffer_size(10, 10, 0);
    ASSERT(size == 0);

    /* Test negative values */
    size = b3d_buffer_size(-1, 10, sizeof(uint32_t));
    ASSERT(size == 0);
    size = b3d_buffer_size(10, -1, sizeof(uint32_t));
    ASSERT(size == 0);

    return 1;
}

/* Test b3d_set_model_matrix API */
TEST(api_set_model_matrix)
{
    const int width = 32, height = 32;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    ASSERT(pixels);
    ASSERT(depth);

    bool init_ok = b3d_init(pixels, depth, width, height, 65.0f);
    ASSERT(init_ok);
    b3d_clear();

    /* Create a custom matrix (translation by 5, 10, 15) */
    float custom_matrix[16] = {
        1.0f, 0.0f,  0.0f,  0.0f, /* column 0 */
        0.0f, 1.0f,  0.0f,  0.0f, /* column 1 */
        0.0f, 0.0f,  1.0f,  0.0f, /* column 2 */
        5.0f, 10.0f, 15.0f, 1.0f, /* column 3 (translation) */
    };

    /* Set the custom matrix */
    b3d_set_model_matrix(custom_matrix);

    /* Retrieve and verify */
    float retrieved[16];
    b3d_get_model_matrix(retrieved);
    ASSERT(fabsf(retrieved[12] - 5.0f) < 0.01f);
    ASSERT(fabsf(retrieved[13] - 10.0f) < 0.01f);
    ASSERT(fabsf(retrieved[14] - 15.0f) < 0.01f);

    free(pixels);
    free(depth);
    return 1;
}

/* Test triangle rendering return value */
TEST(api_triangle_return)
{
    const int width = 64, height = 64;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    ASSERT(pixels);
    ASSERT(depth);

    bool init_ok = b3d_init(pixels, depth, width, height, 65.0f);
    ASSERT(init_ok);
    b3d_set_camera(&(b3d_camera_t) {0.0f, 0.0f, -2.3f, 0.0f, 0.0f, 0.0f});
    b3d_clear();
    b3d_reset();

    /* Render a visible front-facing triangle - should return true */
    /* This is the same setup as api_render_ascii which we know works */
    b3d_rotate_y(0.5f);
    b3d_rotate_x(0.3f);
    bool result = b3d_triangle(
        &(b3d_tri_t) {
            {{-0.5f, -0.5f, -0.5f}, {-0.5f, 0.5f, -0.5f}, {0.5f, 0.5f, -0.5f}}},
        0xfcd0a1);
    ASSERT(result == true);

    free(pixels);
    free(depth);
    return 1;
}

/* Test degenerate triangle handling */
TEST(api_degenerate_triangles)
{
    const int width = 32, height = 32;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    ASSERT(pixels);
    ASSERT(depth);

    bool init_ok = b3d_init(pixels, depth, width, height, 65.0f);
    ASSERT(init_ok);
    b3d_set_camera(&(b3d_camera_t) {0.0f, 0.0f, -3.0f, 0.0f, 0.0f, 0.0f});
    b3d_clear();
    b3d_reset();

    /* Test that degenerate triangles are handled gracefully (no crash).
     * The specific return value depends on pipeline stages - we just verify
     * the function doesn't crash and produces consistent behavior.
     */

    /* Collinear points (line) - minimal area, handled by rasterizer */
    (void) b3d_triangle(
        &(b3d_tri_t) {
            {{0.0f, 0.0f, 0.5f}, {0.5f, 0.0f, 0.5f}, {1.0f, 0.0f, 0.5f}}},
        0xffffff);

    /* All same point - zero area triangle */
    (void) b3d_triangle(
        &(b3d_tri_t) {
            {{0.0f, 0.0f, 0.5f}, {0.0f, 0.0f, 0.5f}, {0.0f, 0.0f, 0.5f}}},
        0xffffff);

    /* Verify no pixels were written - check all pixels for robustness */
    const size_t pixel_count = (size_t) width * (size_t) height;
    size_t non_zero = 0;
    for (size_t i = 0; i < pixel_count; i++) {
        if (pixels[i] != 0)
            non_zero++;
    }
    /* Degenerate triangles should produce no pixel output */
    ASSERT(non_zero == 0);

    free(pixels);
    free(depth);
    return 1;
}

/* Test near plane clipping */
TEST(api_near_plane_clip)
{
    const int width = 64, height = 64;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    ASSERT(pixels);
    ASSERT(depth);

    bool init_ok = b3d_init(pixels, depth, width, height, 65.0f);
    ASSERT(init_ok);
    b3d_set_camera(&(b3d_camera_t) {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f});
    b3d_clear();
    b3d_reset();

    const size_t pixel_count = (size_t) width * (size_t) height;

    /* Test that triangles partially behind near plane are clipped gracefully.
     * The near plane in view space is at z = 0.1.
     * We render a triangle and verify no crash occurs.
     */

    /* First, render a fully visible triangle as baseline */
    bool result1 = b3d_triangle(
        &(b3d_tri_t) {
            {{-0.5f, -0.5f, 0.5f}, {0.0f, 0.5f, 0.5f}, {0.5f, -0.5f, 0.5f}}},
        0x00ff00);

    /* Count pixels after first triangle */
    size_t pixels_after_first = 0;
    for (size_t i = 0; i < pixel_count; i++) {
        if (pixels[i] != 0)
            pixels_after_first++;
    }

    ASSERT(result1 == true);
    /* Verify at least some pixels were rendered */
    ASSERT(pixels_after_first > 0);

    /* Now render a triangle fully behind the near plane: it should be culled
     * and produce no pixel output.
     */
    b3d_clear();
    bool clipped = b3d_triangle(
        &(b3d_tri_t) {
            {{-0.5f, -0.5f, 0.05f}, {0.0f, 0.5f, 0.05f}, {0.5f, -0.5f, 0.05f}}},
        0xff0000);
    size_t pixels_after_clipped = 0;
    for (size_t i = 0; i < pixel_count; i++) {
        if (pixels[i] != 0)
            pixels_after_clipped++;
    }
    ASSERT(clipped == false);
    ASSERT(pixels_after_clipped == 0);

    /* Triangle straddling near plane: should be clipped but still render. */
    b3d_clear();
    bool clipped_visible = b3d_triangle(
        &(b3d_tri_t) {
            {{-0.5f, -0.5f, 0.05f}, {0.0f, 0.5f, 0.2f}, {0.5f, -0.5f, 0.2f}}},
        0x0000ff);
    size_t pixels_after_clip = 0;
    for (size_t i = 0; i < pixel_count; i++) {
        if (pixels[i] != 0)
            pixels_after_clip++;
    }
    ASSERT(clipped_visible == true);
    ASSERT(pixels_after_clip > 0);

    free(pixels);
    free(depth);
    return 1;
}

/* Test screen projection with various positions */
TEST(api_to_screen_extended)
{
    const int width = 64, height = 64;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    ASSERT(pixels);
    ASSERT(depth);

    b3d_init(pixels, depth, width, height, 65.0f);
    b3d_set_camera(&(b3d_camera_t) {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f});
    b3d_clear();

    int sx, sy;

    /* Point directly in front - should project to center */
    bool result = b3d_to_screen(0.0f, 0.0f, 1.0f, &sx, &sy);
    ASSERT(result == true);
    ASSERT(sx >= width / 2 - 2 && sx <= width / 2 + 2);
    ASSERT(sy >= height / 2 - 2 && sy <= height / 2 + 2);

    /* Point to the right - should project to right side */
    result = b3d_to_screen(1.0f, 0.0f, 1.0f, &sx, &sy);
    ASSERT(result == true);
    ASSERT(sx > width / 2);

    /* Point above - should project to upper half (y inverted in screen coords)
     */
    result = b3d_to_screen(0.0f, 1.0f, 1.0f, &sx, &sy);
    ASSERT(result == true);
    ASSERT(sy < height / 2);

    /* Point behind camera - should fail */
    result = b3d_to_screen(0.0f, 0.0f, -1.0f, &sx, &sy);
    ASSERT(result == false);

    free(pixels);
    free(depth);
    return 1;
}

/* Test b3d_init validation */
TEST(api_init_validation)
{
    const int width = 32, height = 32;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    /* Valid initialization */
    bool result = b3d_init(pixels, depth, width, height, 65.0f);
    ASSERT(result == true);

    /* NULL pixel buffer should fail */
    result = b3d_init(NULL, depth, width, height, 65.0f);
    ASSERT(result == false);

    /* NULL depth buffer should fail */
    result = b3d_init(pixels, NULL, width, height, 65.0f);
    ASSERT(result == false);

    /* Zero width should fail */
    result = b3d_init(pixels, depth, 0, height, 65.0f);
    ASSERT(result == false);

    /* Zero height should fail */
    result = b3d_init(pixels, depth, width, 0, 65.0f);
    ASSERT(result == false);

    /* Negative dimensions should fail */
    result = b3d_init(pixels, depth, -1, height, 65.0f);
    ASSERT(result == false);

    free(pixels);
    free(depth);
    return 1;
}

/* Test lighting API - direction setting and getting */
TEST(api_lighting_direction)
{
    const int width = 32, height = 32;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    ASSERT(pixels);
    ASSERT(depth);

    bool init_ok = b3d_init(pixels, depth, width, height, 65.0f);
    ASSERT(init_ok);

    float x, y, z;

    /* Test default light direction (0, 0, 1) toward +Z */
    b3d_get_light_direction(&x, &y, &z);
    ASSERT(fabsf(x - 0.0f) < 0.0001f);
    ASSERT(fabsf(y - 0.0f) < 0.0001f);
    ASSERT(fabsf(z - 1.0f) < 0.0001f);

    /* Test setting a new direction (auto-normalized) */
    b3d_set_light_direction(1.0f, 1.0f, 1.0f);
    b3d_get_light_direction(&x, &y, &z);
    float expected = 1.0f / sqrtf(3.0f); /* ~0.577 */
    ASSERT(fabsf(x - expected) < 0.01f);
    ASSERT(fabsf(y - expected) < 0.01f);
    ASSERT(fabsf(z - expected) < 0.01f);

    /* Test zero-length vector rejection (keeps previous value) */
    b3d_set_light_direction(0.0f, 0.0f, 0.0f);
    b3d_get_light_direction(&x, &y, &z);
    ASSERT(fabsf(x - expected) < 0.01f); /* Should still be previous value */

    /* Test NULL-safe getters */
    b3d_get_light_direction(&x, NULL, NULL);
    ASSERT(fabsf(x - expected) < 0.01f);
    b3d_get_light_direction(NULL, &y, NULL);
    ASSERT(fabsf(y - expected) < 0.01f);
    b3d_get_light_direction(NULL, NULL, &z);
    ASSERT(fabsf(z - expected) < 0.01f);

    free(pixels);
    free(depth);
    return 1;
}

/* Test lighting API - ambient setting and getting */
TEST(api_lighting_ambient)
{
    const int width = 32, height = 32;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    ASSERT(pixels);
    ASSERT(depth);

    bool init_ok = b3d_init(pixels, depth, width, height, 65.0f);
    ASSERT(init_ok);

    /* Test default ambient (0.2) */
    float ambient = b3d_get_ambient();
    ASSERT(fabsf(ambient - 0.2f) < 0.0001f);

    /* Test setting ambient within range */
    b3d_set_ambient(0.5f);
    ambient = b3d_get_ambient();
    ASSERT(fabsf(ambient - 0.5f) < 0.0001f);

    /* Test clamping to 0 */
    b3d_set_ambient(-0.5f);
    ambient = b3d_get_ambient();
    ASSERT(fabsf(ambient - 0.0f) < 0.0001f);

    /* Test clamping to 1 */
    b3d_set_ambient(1.5f);
    ambient = b3d_get_ambient();
    ASSERT(fabsf(ambient - 1.0f) < 0.0001f);

    /* Test edge values */
    b3d_set_ambient(0.0f);
    ASSERT(fabsf(b3d_get_ambient() - 0.0f) < 0.0001f);
    b3d_set_ambient(1.0f);
    ASSERT(fabsf(b3d_get_ambient() - 1.0f) < 0.0001f);

    free(pixels);
    free(depth);
    return 1;
}

/* Test lit triangle rendering */
TEST(api_triangle_lit)
{
    const int width = 64, height = 64;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    ASSERT(pixels);
    ASSERT(depth);

    bool init_ok = b3d_init(pixels, depth, width, height, 65.0f);
    ASSERT(init_ok);
    b3d_set_camera(&(b3d_camera_t) {0.0f, 0.0f, -2.0f, 0.0f, 0.0f, 0.0f});
    b3d_clear();
    b3d_reset();

    /* Set light direction toward +Z (facing the triangle) */
    b3d_set_light_direction(0.0f, 0.0f, 1.0f);
    b3d_set_ambient(0.2f);

    /* Render a lit triangle facing the camera with normal toward -Z */
    bool result =
        b3d_triangle_lit(&(b3d_tri_t) {{
                             {-0.5f, -0.5f, 0.0f},
                             {0.0f, 0.5f, 0.0f},
                             {0.5f, -0.5f, 0.0f},
                         }},
                         0.0f, 0.0f, -1.0f, /* Normal pointing toward camera */
                         0xffffff);         /* White base color */
    ASSERT(result == true);

    /* Verify pixels were rendered */
    const size_t pixel_count = (size_t) width * (size_t) height;
    size_t non_zero = 0;
    for (size_t i = 0; i < pixel_count; i++) {
        if (pixels[i] != 0)
            non_zero++;
    }
    ASSERT(non_zero > 0);

    /* Verify shading occurred (pixels should not be pure white due to lighting)
     * With two-sided lighting: |dot(-Z, +Z)| = 1, so intensity = 0.2 + 0.8*1 =
     * 1 Actually full intensity, so pixels should be white (0xffffff)
     */
    size_t center =
        (size_t) (height / 2) * (size_t) width + (size_t) (width / 2);
    uint32_t center_color = pixels[center];
    ASSERT(center_color == 0xffffff); /* Full intensity */

    /* Test with perpendicular normal (should get ambient-only shading) */
    b3d_clear();
    b3d_set_light_direction(1.0f, 0.0f, 0.0f); /* Light from +X */
    result = b3d_triangle_lit(
        &(b3d_tri_t) {{
            {-0.5f, -0.5f, 0.0f},
            {0.0f, 0.5f, 0.0f},
            {0.5f, -0.5f, 0.0f},
        }},
        0.0f, 0.0f, -1.0f, /* Normal toward -Z, perpendicular to light */
        0xffffff);
    ASSERT(result == true);

    /* With perpendicular normal: dot = 0, intensity = ambient = 0.2 */
    center_color = pixels[center];
    /* 0.2 * 255 = 51 = 0x33, so expect ~0x333333 */
    unsigned r = (center_color >> 16) & 0xff;
    unsigned g = (center_color >> 8) & 0xff;
    unsigned b = center_color & 0xff;
    ASSERT(r >= 48 && r <= 56); /* ~0x33 with tolerance */
    ASSERT(g >= 48 && g <= 56);
    ASSERT(b >= 48 && b <= 56);

    free(pixels);
    free(depth);
    return 1;
}

/* Test depth buffer functionality */
TEST(api_depth_buffer)
{
    const int width = 32, height = 32;
    uint32_t *pixels =
        malloc((size_t) width * (size_t) height * sizeof(uint32_t));
    b3d_depth_t *depth =
        malloc((size_t) width * (size_t) height * sizeof(b3d_depth_t));

    ASSERT(pixels);
    ASSERT(depth);

    bool init_ok = b3d_init(pixels, depth, width, height, 65.0f);
    ASSERT(init_ok);
    b3d_set_camera(&(b3d_camera_t) {0.0f, 0.0f, -5.0f, 0.0f, 0.0f, 0.0f});
    b3d_clear();
    b3d_reset();

    /* With camera at z=-5, objects at higher z are FARTHER away:
     * - z=0.5 is 5.5 units from camera (FAR)
     * - z=0.0 is 5.0 units from camera (NEAR)
     * Use CCW winding (v0->v2->v1) for front-facing triangles from camera view.
     */

    /* Render a large green triangle far from camera (z=0.5) */
    uint32_t far_color = 0x00ff00;
    b3d_triangle(&(b3d_tri_t) {{
                     {-1.0f, -1.0f, 0.5f},
                     {0.0f, 1.0f, 0.5f},
                     {1.0f, -1.0f, 0.5f},
                 }},
                 far_color);

    /* Check that center pixel is the far color */
    size_t center =
        (size_t) (height / 2) * (size_t) width + (size_t) (width / 2);
    uint32_t center_before = pixels[center];
    ASSERT(center_before == far_color);

    /* Render a smaller red triangle closer to camera (z=0.0) */
    uint32_t near_color = 0xff0000;
    b3d_triangle(&(b3d_tri_t) {{
                     {-0.5f, -0.5f, 0.0f},
                     {0.0f, 0.5f, 0.0f},
                     {0.5f, -0.5f, 0.0f},
                 }},
                 near_color);

    /* Center pixel should now be the near color (closer triangle wins) */
    uint32_t center_after = pixels[center];
    ASSERT(center_after == near_color);

    /* Render in reverse order: near then far */
    b3d_clear();

    /* Render small red triangle first (closer, z=0.0) */
    b3d_triangle(&(b3d_tri_t) {{
                     {-0.5f, -0.5f, 0.0f},
                     {0.0f, 0.5f, 0.0f},
                     {0.5f, -0.5f, 0.0f},
                 }},
                 near_color);

    /* Render large green triangle second (farther, z=0.5) */
    b3d_triangle(&(b3d_tri_t) {{
                     {-1.0f, -1.0f, 0.5f},
                     {0.0f, 1.0f, 0.5f},
                     {1.0f, -1.0f, 0.5f},
                 }},
                 far_color);

    /* Center pixel should STILL be the near color due to depth test */
    uint32_t center_reverse = pixels[center];
    ASSERT(center_reverse == near_color);

    free(pixels);
    free(depth);
    return 1;
}

int main(void)
{
    printf(ANSI_BOLD "B3D API Validation Tests\n" ANSI_RESET);
    printf("======================\n");

    SECTION_BEGIN("API Initialization");
    RUN_TEST(api_init);
    RUN_TEST(api_init_validation);
    SECTION_END();

    SECTION_BEGIN("API Transformations");
    RUN_TEST(api_transform);
    RUN_TEST(api_set_model_matrix);
    SECTION_END();

    SECTION_BEGIN("API Camera");
    RUN_TEST(api_camera);
    SECTION_END();

    SECTION_BEGIN("API State Queries");
    RUN_TEST(api_state_queries);
    SECTION_END();

    SECTION_BEGIN("API Matrix Stack");
    RUN_TEST(api_matrix_stack);
    SECTION_END();

    SECTION_BEGIN("API Rendering");
    RUN_TEST(api_render_ascii);
    RUN_TEST(api_triangle_return);
    RUN_TEST(api_degenerate_triangles);
    RUN_TEST(api_depth_buffer);
    SECTION_END();

    SECTION_BEGIN("API Lighting");
    RUN_TEST(api_lighting_direction);
    RUN_TEST(api_lighting_ambient);
    RUN_TEST(api_triangle_lit);
    SECTION_END();

    SECTION_BEGIN("API Clipping");
    RUN_TEST(api_near_plane_clip);
    RUN_TEST(api_clip_drop_count);
    SECTION_END();

    SECTION_BEGIN("API Screen Projection");
    RUN_TEST(api_to_screen);
    RUN_TEST(api_to_screen_extended);
    SECTION_END();

    SECTION_BEGIN("API Buffer Size");
    RUN_TEST(api_buffer_size);
    SECTION_END();

    printf("======================\n");
    if (tests_passed == tests_run)
        printf(ANSI_GREEN "All %d tests passed" ANSI_RESET "\n", tests_run);
    else
        printf(ANSI_RED "%d/%d tests failed" ANSI_RESET "\n",
               tests_run - tests_passed, tests_run);

    return tests_passed == tests_run ? 0 : 1;
}
