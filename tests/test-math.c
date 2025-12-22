/*
 * B3D is freely redistributable under the MIT License. See the file
 * "LICENSE" for information on usage and redistribution of this file.
 */

/* Unit tests for math-toolkit.h */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../src/math-toolkit.h"

static int tests_run = 0;
static int tests_passed = 0;
static int section_passed = 0;
static int section_total = 0;

/* ANSI color codes for terminal output */
#define ANSI_GREEN "\033[32m"
#define ANSI_RED "\033[31m"
#define ANSI_BOLD "\033[1m"
#define ANSI_RESET "\033[0m"

#define TEST(name) \
    static int test_##name(void)

#define RUN_TEST(name)                                   \
    do {                                                 \
        tests_run++;                                     \
        section_total++;                                 \
        if (test_##name()) {                             \
            tests_passed++;                              \
            section_passed++;                            \
        }                                                \
    } while (0)

#define SECTION_BEGIN(name)                              \
    do {                                                 \
        section_passed = 0;                              \
        section_total = 0;                               \
        printf("  %-44s", name);                         \
        fflush(stdout);                                  \
    } while (0)

#define SECTION_END()                                                     \
    do {                                                                  \
        if (section_passed == section_total)                              \
            printf("[ " ANSI_GREEN "OK" ANSI_RESET " ]\n");               \
        else                                                              \
            printf("[ " ANSI_RED "FAIL %d/%d" ANSI_RESET " ]\n",          \
                   section_passed, section_total);                        \
    } while (0)

#define ASSERT(cond)       \
    do {                   \
        if (!(cond))       \
            return 0;      \
    } while (0)

#define ASSERT_NEAR(a, b, eps) \
    ASSERT(fabsf((float)(a) - (float)(b)) < (eps))

/*
 * Fixed-point conversion tests
 */
TEST(fp_int_conversion)
{
    /* Integer to fixed-point and back (positive values) */
    ASSERT(B3D_FP_TO_INT(B3D_INT_TO_FP(0)) == 0);
    ASSERT(B3D_FP_TO_INT(B3D_INT_TO_FP(1)) == 1);
    ASSERT(B3D_FP_TO_INT(B3D_INT_TO_FP(100)) == 100);
    ASSERT(B3D_FP_TO_INT(B3D_INT_TO_FP(1000)) == 1000);

    /* Negative values via float conversion (avoids shift UB) */
    ASSERT(B3D_FP_TO_INT(B3D_FLOAT_TO_FP(-1.0f)) == -1);
    ASSERT(B3D_FP_TO_INT(B3D_FLOAT_TO_FP(-100.0f)) == -100);
    return 1;
}

TEST(fp_float_conversion)
{
    /* Float to fixed-point and back */
    ASSERT_NEAR(B3D_FP_TO_FLOAT(B3D_FLOAT_TO_FP(0.0f)), 0.0f, 0.001f);
    ASSERT_NEAR(B3D_FP_TO_FLOAT(B3D_FLOAT_TO_FP(1.0f)), 1.0f, 0.001f);
    ASSERT_NEAR(B3D_FP_TO_FLOAT(B3D_FLOAT_TO_FP(0.5f)), 0.5f, 0.001f);
    ASSERT_NEAR(B3D_FP_TO_FLOAT(B3D_FLOAT_TO_FP(-1.0f)), -1.0f, 0.001f);
    ASSERT_NEAR(B3D_FP_TO_FLOAT(B3D_FLOAT_TO_FP(3.14159f)), 3.14159f, 0.001f);
    return 1;
}

/*
 * Fixed-point arithmetic tests
 */
TEST(fp_multiplication)
{
    b3d_scalar_t a = B3D_FLOAT_TO_FP(2.0f);
    b3d_scalar_t b = B3D_FLOAT_TO_FP(3.0f);
    b3d_scalar_t result = B3D_FP_MUL(a, b);
    ASSERT_NEAR(B3D_FP_TO_FLOAT(result), 6.0f, 0.001f);

    /* Test with fractional values */
    a = B3D_FLOAT_TO_FP(0.5f);
    b = B3D_FLOAT_TO_FP(0.5f);
    result = B3D_FP_MUL(a, b);
    ASSERT_NEAR(B3D_FP_TO_FLOAT(result), 0.25f, 0.001f);

    /* Test negative multiplication */
    a = B3D_FLOAT_TO_FP(-2.0f);
    b = B3D_FLOAT_TO_FP(3.0f);
    result = B3D_FP_MUL(a, b);
    ASSERT_NEAR(B3D_FP_TO_FLOAT(result), -6.0f, 0.001f);

    return 1;
}

TEST(fp_division)
{
    b3d_scalar_t a = B3D_FLOAT_TO_FP(6.0f);
    b3d_scalar_t b = B3D_FLOAT_TO_FP(2.0f);
    b3d_scalar_t result = B3D_FP_DIV(a, b);
    ASSERT_NEAR(B3D_FP_TO_FLOAT(result), 3.0f, 0.001f);

    /* Test division by zero returns 0 */
    result = B3D_FP_DIV(a, 0);
    ASSERT(result == 0);

    /* Test negative division */
    a = B3D_FLOAT_TO_FP(-6.0f);
    b = B3D_FLOAT_TO_FP(2.0f);
    result = B3D_FP_DIV(a, b);
    ASSERT_NEAR(B3D_FP_TO_FLOAT(result), -3.0f, 0.001f);

    return 1;
}

TEST(fp_floor)
{
    ASSERT(B3D_FP_TO_INT(B3D_FP_FLOOR(B3D_FLOAT_TO_FP(1.9f))) == 1);
    ASSERT(B3D_FP_TO_INT(B3D_FP_FLOOR(B3D_FLOAT_TO_FP(1.1f))) == 1);
    ASSERT(B3D_FP_TO_INT(B3D_FP_FLOOR(B3D_FLOAT_TO_FP(1.0f))) == 1);
    ASSERT(B3D_FP_TO_INT(B3D_FP_FLOOR(B3D_FLOAT_TO_FP(0.9f))) == 0);
    return 1;
}

/*
 * Trigonometric function tests
 */
TEST(fp_sin_basic)
{
    /* sin(0) = 0 */
    ASSERT_NEAR(B3D_FP_TO_FLOAT(b3d_fp_sin(0)), 0.0f, 0.01f);

    /* sin(pi/2) = 1 */
    ASSERT_NEAR(B3D_FP_TO_FLOAT(b3d_fp_sin(B3D_FP_PI_HALF)), 1.0f, 0.01f);

    /* sin(pi) = 0 */
    ASSERT_NEAR(B3D_FP_TO_FLOAT(b3d_fp_sin(B3D_FP_PI)), 0.0f, 0.01f);

    /* sin(3*pi/2) = -1 */
    b3d_scalar_t three_pi_half = B3D_FP_PI + B3D_FP_PI_HALF;
    ASSERT_NEAR(B3D_FP_TO_FLOAT(b3d_fp_sin(three_pi_half)), -1.0f, 0.01f);

    return 1;
}

TEST(fp_sin_negative)
{
    /* sin(-x) = -sin(x) */
    b3d_scalar_t x = B3D_FP_PI_HALF;
    float sin_pos = B3D_FP_TO_FLOAT(b3d_fp_sin(x));
    float sin_neg = B3D_FP_TO_FLOAT(b3d_fp_sin(-x));
    ASSERT_NEAR(sin_neg, -sin_pos, 0.01f);

    return 1;
}

TEST(fp_cos_basic)
{
    /* cos(0) = 1 */
    ASSERT_NEAR(B3D_FP_TO_FLOAT(b3d_fp_cos(0)), 1.0f, 0.01f);

    /* cos(pi/2) = 0 */
    ASSERT_NEAR(B3D_FP_TO_FLOAT(b3d_fp_cos(B3D_FP_PI_HALF)), 0.0f, 0.01f);

    /* cos(pi) = -1 */
    ASSERT_NEAR(B3D_FP_TO_FLOAT(b3d_fp_cos(B3D_FP_PI)), -1.0f, 0.01f);

    return 1;
}

TEST(fp_sqrt_basic)
{
    /* sqrt(0) = 0 */
    ASSERT(b3d_fp_sqrt(0) == 0);

    /* sqrt(1) = 1 */
    ASSERT_NEAR(B3D_FP_TO_FLOAT(b3d_fp_sqrt(B3D_FLOAT_TO_FP(1.0f))), 1.0f,
                0.01f);

    /* sqrt(4) = 2 */
    ASSERT_NEAR(B3D_FP_TO_FLOAT(b3d_fp_sqrt(B3D_FLOAT_TO_FP(4.0f))), 2.0f,
                0.01f);

    /* sqrt(2) ~= 1.414 */
    ASSERT_NEAR(B3D_FP_TO_FLOAT(b3d_fp_sqrt(B3D_FLOAT_TO_FP(2.0f))), 1.414f,
                0.01f);

    /* sqrt(negative) = 0 */
    ASSERT(b3d_fp_sqrt(B3D_FLOAT_TO_FP(-1.0f)) == 0);

    return 1;
}

TEST(fp_abs_basic)
{
    ASSERT(b3d_fp_abs(B3D_FLOAT_TO_FP(5.0f)) == B3D_FLOAT_TO_FP(5.0f));
    ASSERT(b3d_fp_abs(B3D_FLOAT_TO_FP(-5.0f)) == B3D_FLOAT_TO_FP(5.0f));
    ASSERT(b3d_fp_abs(0) == 0);
    return 1;
}

/*
 * Vector operation tests
 */
TEST(vec_dot)
{
    b3d_vec_t a = {1, 0, 0, 1};
    b3d_vec_t b = {1, 0, 0, 1};
    ASSERT_NEAR(b3d_vec_dot(a, b), 1.0f, 0.001f);

    b3d_vec_t c = {0, 1, 0, 1};
    ASSERT_NEAR(b3d_vec_dot(a, c), 0.0f, 0.001f);

    b3d_vec_t d = {1, 2, 3, 1};
    b3d_vec_t e = {4, 5, 6, 1};
    ASSERT_NEAR(b3d_vec_dot(d, e), 32.0f, 0.001f); /* 1*4 + 2*5 + 3*6 = 32 */

    return 1;
}

TEST(vec_cross)
{
    b3d_vec_t x = {1, 0, 0, 1};
    b3d_vec_t y = {0, 1, 0, 1};
    b3d_vec_t z = b3d_vec_cross(x, y);
    ASSERT_NEAR(z.x, 0.0f, 0.001f);
    ASSERT_NEAR(z.y, 0.0f, 0.001f);
    ASSERT_NEAR(z.z, 1.0f, 0.001f);

    return 1;
}

TEST(vec_length)
{
    b3d_vec_t v = {3, 4, 0, 1};
    ASSERT_NEAR(b3d_vec_length(v), 5.0f, 0.001f);

    b3d_vec_t unit = {1, 0, 0, 1};
    ASSERT_NEAR(b3d_vec_length(unit), 1.0f, 0.001f);

    return 1;
}

TEST(vec_norm)
{
    b3d_vec_t v = {3, 4, 0, 1};
    b3d_vec_t n = b3d_vec_norm(v);
    ASSERT_NEAR(b3d_vec_length(n), 1.0f, 0.001f);
    ASSERT_NEAR(n.x, 0.6f, 0.001f);
    ASSERT_NEAR(n.y, 0.8f, 0.001f);

    return 1;
}

/*
 * Matrix operation tests
 */
TEST(mat_identity)
{
    b3d_mat_t m = b3d_mat_ident();
    ASSERT(b3d_mat_is_identity(m));

    b3d_vec_t v = {1, 2, 3, 1};
    b3d_vec_t result = b3d_mat_mul_vec(m, v);
    ASSERT_NEAR(result.x, v.x, 0.001f);
    ASSERT_NEAR(result.y, v.y, 0.001f);
    ASSERT_NEAR(result.z, v.z, 0.001f);

    return 1;
}

TEST(mat_translation)
{
    b3d_mat_t m = b3d_mat_trans(10, 20, 30);
    b3d_vec_t v = {0, 0, 0, 1};
    b3d_vec_t result = b3d_mat_mul_vec(m, v);
    ASSERT_NEAR(result.x, 10.0f, 0.001f);
    ASSERT_NEAR(result.y, 20.0f, 0.001f);
    ASSERT_NEAR(result.z, 30.0f, 0.001f);

    return 1;
}

TEST(mat_scale)
{
    b3d_mat_t m = b3d_mat_scale(2, 3, 4);
    b3d_vec_t v = {1, 1, 1, 1};
    b3d_vec_t result = b3d_mat_mul_vec(m, v);
    ASSERT_NEAR(result.x, 2.0f, 0.001f);
    ASSERT_NEAR(result.y, 3.0f, 0.001f);
    ASSERT_NEAR(result.z, 4.0f, 0.001f);

    return 1;
}

TEST(mat_rotation_x)
{
    /* Rotate 90 degrees around X axis: Y becomes Z, Z becomes -Y */
    b3d_mat_t m = b3d_mat_rot_x(B3D_PI / 2.0f);
    b3d_vec_t v = {0, 1, 0, 1};
    b3d_vec_t result = b3d_mat_mul_vec(m, v);
    ASSERT_NEAR(result.x, 0.0f, 0.01f);
    ASSERT_NEAR(result.y, 0.0f, 0.01f);
    ASSERT_NEAR(result.z, 1.0f, 0.01f);

    return 1;
}

/*
 * Depth buffer conversion tests
 */
TEST(depth_conversion)
{
    float d = 0.5f;
    b3d_depth_t stored = b3d_depth_from_float(d);
    float restored = b3d_depth_to_float(stored);
    ASSERT_NEAR(restored, d, 0.001f);

    /* Test boundaries */
    stored = b3d_depth_from_float(0.0f);
    restored = b3d_depth_to_float(stored);
    ASSERT_NEAR(restored, 0.0f, 0.001f);

    stored = b3d_depth_from_float(1.0f);
    restored = b3d_depth_to_float(stored);
    ASSERT_NEAR(restored, 1.0f, 0.001f);

    return 1;
}

int main(void)
{
    printf(ANSI_BOLD "B3D Math Toolkit Tests\n" ANSI_RESET);
    printf("======================\n");

    SECTION_BEGIN("Fixed-point conversions");
    RUN_TEST(fp_int_conversion);
    RUN_TEST(fp_float_conversion);
    SECTION_END();

    SECTION_BEGIN("Fixed-point arithmetic");
    RUN_TEST(fp_multiplication);
    RUN_TEST(fp_division);
    RUN_TEST(fp_floor);
    SECTION_END();

    SECTION_BEGIN("Trigonometric functions");
    RUN_TEST(fp_sin_basic);
    RUN_TEST(fp_sin_negative);
    RUN_TEST(fp_cos_basic);
    RUN_TEST(fp_sqrt_basic);
    RUN_TEST(fp_abs_basic);
    SECTION_END();

    SECTION_BEGIN("Vector operations");
    RUN_TEST(vec_dot);
    RUN_TEST(vec_cross);
    RUN_TEST(vec_length);
    RUN_TEST(vec_norm);
    SECTION_END();

    SECTION_BEGIN("Matrix operations");
    RUN_TEST(mat_identity);
    RUN_TEST(mat_translation);
    RUN_TEST(mat_scale);
    RUN_TEST(mat_rotation_x);
    SECTION_END();

    SECTION_BEGIN("Depth buffer");
    RUN_TEST(depth_conversion);
    SECTION_END();

    printf("======================\n");
    if (tests_passed == tests_run)
        printf(ANSI_GREEN "All %d tests passed" ANSI_RESET "\n", tests_run);
    else
        printf(ANSI_RED "%d/%d tests failed" ANSI_RESET "\n",
               tests_run - tests_passed, tests_run);

    return tests_passed == tests_run ? 0 : 1;
}
