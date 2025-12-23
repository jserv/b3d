/*
 * B3D is freely redistributable under the MIT License. See the file
 * "LICENSE" for information on usage and redistribution of this file.
 */

/* Unit tests for DSL-generated math functions (math_gen.h) */

#include <math.h>
#include <stdio.h>

#include "../src/math-toolkit.h"
#include "../src/math_gen_test.h"

static int tests_run = 0;
static int tests_passed = 0;

#define ANSI_GREEN "\033[32m"
#define ANSI_RED "\033[31m"
#define ANSI_BOLD "\033[1m"
#define ANSI_RESET "\033[0m"

#define TEST(name) static int test_##name(void)
#define RUN_TEST(name)                                      \
    do {                                                    \
        tests_run++;                                        \
        printf("  %-40s", #name);                           \
        if (test_##name()) {                                \
            tests_passed++;                                 \
            printf("[ " ANSI_GREEN "OK" ANSI_RESET " ]\n"); \
        } else {                                            \
            printf("[ " ANSI_RED "FAIL" ANSI_RESET " ]\n"); \
        }                                                   \
    } while (0)

#define ASSERT(cond)  \
    do {              \
        if (!(cond))  \
            return 0; \
    } while (0)

#define ASSERT_NEAR(a, b, eps) ASSERT(fabsf((float) (a) - (float) (b)) < (eps))

/* Test generated vec_dot matches existing implementation */
TEST(gen_vec_dot)
{
    b3d_vec_t a = {1, 2, 3, 1};
    b3d_vec_t b = {4, 5, 6, 1};

    /* Both should compute 1*4 + 2*5 + 3*6 = 32 */
    float existing = b3d_vec_dot(a, b);
    float generated = b3d_vec_dot_gen(a, b);

    ASSERT_NEAR(existing, 32.0f, 0.001f);
    ASSERT_NEAR(generated, 32.0f, 0.001f);
    return 1;
}

/* Test generated vec_cross matches existing implementation */
TEST(gen_vec_cross)
{
    b3d_vec_t x = {1, 0, 0, 1};
    b3d_vec_t y = {0, 1, 0, 1};

    b3d_vec_t existing = b3d_vec_cross(x, y);
    b3d_vec_t generated = b3d_vec_cross_gen(x, y);

    /* x cross y = z */
    ASSERT_NEAR(existing.z, 1.0f, 0.001f);
    ASSERT_NEAR(generated.z, 1.0f, 0.001f);
    return 1;
}

/* Test generated vec_length */
TEST(gen_vec_length)
{
    b3d_vec_t v = {3, 4, 0, 1};

    float existing = b3d_vec_length(v);
    float generated = b3d_vec_length_gen(v);

    ASSERT_NEAR(existing, 5.0f, 0.001f);
    ASSERT_NEAR(generated, 5.0f, 0.001f);
    return 1;
}

/* Test generated vec_norm */
TEST(gen_vec_norm)
{
    b3d_vec_t v = {3, 4, 0, 1};

    b3d_vec_t existing = b3d_vec_norm(v);
    b3d_vec_t generated = b3d_vec_norm_gen(v);

    ASSERT_NEAR(existing.x, 0.6f, 0.001f);
    ASSERT_NEAR(existing.y, 0.8f, 0.001f);
    ASSERT_NEAR(generated.x, 0.6f, 0.001f);
    ASSERT_NEAR(generated.y, 0.8f, 0.001f);
    return 1;
}

/* Test generated vec_add */
TEST(gen_vec_add)
{
    b3d_vec_t a = {1, 2, 3, 1};
    b3d_vec_t b = {4, 5, 6, 1};

    b3d_vec_t existing = b3d_vec_add(a, b);
    b3d_vec_t generated = b3d_vec_add_gen(a, b);

    ASSERT_NEAR(existing.x, 5.0f, 0.001f);
    ASSERT_NEAR(existing.y, 7.0f, 0.001f);
    ASSERT_NEAR(existing.z, 9.0f, 0.001f);
    ASSERT_NEAR(generated.x, 5.0f, 0.001f);
    ASSERT_NEAR(generated.y, 7.0f, 0.001f);
    ASSERT_NEAR(generated.z, 9.0f, 0.001f);
    return 1;
}

/* Test generated vec_sub */
TEST(gen_vec_sub)
{
    b3d_vec_t a = {4, 5, 6, 1};
    b3d_vec_t b = {1, 2, 3, 1};

    b3d_vec_t existing = b3d_vec_sub(a, b);
    b3d_vec_t generated = b3d_vec_sub_gen(a, b);

    ASSERT_NEAR(existing.x, 3.0f, 0.001f);
    ASSERT_NEAR(existing.y, 3.0f, 0.001f);
    ASSERT_NEAR(existing.z, 3.0f, 0.001f);
    ASSERT_NEAR(generated.x, 3.0f, 0.001f);
    ASSERT_NEAR(generated.y, 3.0f, 0.001f);
    ASSERT_NEAR(generated.z, 3.0f, 0.001f);
    return 1;
}

/* Test generated vec_mul */
TEST(gen_vec_mul)
{
    b3d_vec_t v = {1, 2, 3, 1};
    float s = 2.0f;

    b3d_vec_t generated = b3d_vec_mul_gen(v, s);

    ASSERT_NEAR(generated.x, 2.0f, 0.01f);
    ASSERT_NEAR(generated.y, 4.0f, 0.01f);
    ASSERT_NEAR(generated.z, 6.0f, 0.01f);
    return 1;
}

/* Test generated vec_neg */
TEST(gen_vec_neg)
{
    b3d_vec_t v = {1, 2, 3, 1};
    b3d_vec_t result = b3d_vec_neg_gen(v);

    ASSERT_NEAR(result.x, -1.0f, 0.001f);
    ASSERT_NEAR(result.y, -2.0f, 0.001f);
    ASSERT_NEAR(result.z, -3.0f, 0.001f);
    ASSERT_NEAR(result.w, 1.0f, 0.001f); /* w preserved */
    return 1;
}

/* Test zero-length vector normalization */
TEST(gen_vec_norm_zero)
{
    b3d_vec_t v = {0, 0, 0, 1};
    b3d_vec_t result = b3d_vec_norm_gen(v);

    /* Should return zero vector for degenerate input */
    ASSERT_NEAR(result.x, 0.0f, 0.001f);
    ASSERT_NEAR(result.y, 0.0f, 0.001f);
    ASSERT_NEAR(result.z, 0.0f, 0.001f);
    return 1;
}

int main(void)
{
    printf(ANSI_BOLD "B3D Generated Code Tests (math_gen.h)\n" ANSI_RESET);
    printf("======================================\n");

    RUN_TEST(gen_vec_dot);
    RUN_TEST(gen_vec_cross);
    RUN_TEST(gen_vec_length);
    RUN_TEST(gen_vec_norm);
    RUN_TEST(gen_vec_add);
    RUN_TEST(gen_vec_sub);
    RUN_TEST(gen_vec_mul);
    RUN_TEST(gen_vec_neg);
    RUN_TEST(gen_vec_norm_zero);

    printf("======================================\n");
    if (tests_passed == tests_run) {
        printf(ANSI_GREEN "All %d tests passed" ANSI_RESET "\n", tests_run);
    } else {
        printf(ANSI_RED "%d/%d tests failed" ANSI_RESET "\n",
               tests_run - tests_passed, tests_run);
    }

    return tests_passed == tests_run ? 0 : 1;
}
