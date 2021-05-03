/*
 * Adapted Delaunator bench.js program.
 * Note: 1 million degenerate does not work well with our floating point
 * precision. It results in a lot of colinear points, whereas Delaunator
 * being JavaScript is able to maintain uniqueness with double precision.
 * This is the reason we lag waaaay behind in that test, because we're
 * running into a completely different issue on those data sets. However
 * since I'm creating this library mainly for use in game development,
 * I'm not interested at all in double precision. :)
 */

#include "delaunay/delaunay.h"
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define DELC_PI 3.14159265358979323846f

static inline float frand(void);
static inline float pseudoNormal(void);

/* Domain generating functions, all populate pt with npt points */
static inline void degenerate(float *pt, size_t npt);
static inline void gaussian(float *pt, size_t npt);
static inline void grid(float *pt, size_t npt);
static inline void uniform(float *pt, size_t npt);

int
main(void)
{
    /* Seed our random number generator */
    srand((unsigned)time(NULL));

    size_t npt = 100000;
    float  *pt  = calloc(2 * npt, sizeof *pt);
    size_t *del = calloc(DELAUNAY_SZ(npt), sizeof *del);
    if (pt == NULL || del == NULL) {
        perror("Error allocating point buffer");
        return EXIT_FAILURE;
    }

#define BENCHMARK(FN) do {                                   \
        FN(pt, npt);                                         \
        if (clock_gettime(CLOCK_MONOTONIC, &t0) != 0) {      \
            perror("Error measuring time");                  \
            goto error_triangulating;                        \
        }                                                    \
        if (triangulate(del, pt, npt) != 0) {               \
            perror("Error triangulating " #FN);              \
            goto error_triangulating;                        \
        }                                                    \
        if (clock_gettime(CLOCK_MONOTONIC, &t1) != 0) {      \
            perror("Error measuring time");                  \
            goto error_triangulating;                        \
        }                                                    \
        float tD = 1000.0f * t1.tv_sec + 1e-6 * t1.tv_nsec   \
                - (1000.0f * t0.tv_sec + 1e-6 * t0.tv_nsec); \
        if (printf(#FN " took %.2fms\n", tD) < 0) {          \
            perror("Error reporting to stdout");             \
            goto error_triangulating;                        \
        }                                                    \
    } while (0)

    struct timespec t0;
    struct timespec t1;
    BENCHMARK(uniform);
    BENCHMARK(gaussian);
    BENCHMARK(grid);
    BENCHMARK(degenerate);

    free(del);
    free(pt);
    return EXIT_SUCCESS;

error_triangulating:
    free(pt);
    return EXIT_FAILURE;
}

static inline float
frand(void)
{
    return (float)rand() / (float)(RAND_MAX - 1);
}

static inline float
pseudoNormal(void)
{
    float v = frand() + frand() + frand() + frand() + frand() + frand();
    float n = 0.5f * (v - 3.0f) / 3.0f;
    return (n < 1) ? n : 1;
}

static inline void
degenerate(float *pt, size_t npt)
{
    pt[0] = 0;
    pt[1] = 0;
    for (size_t i = 1; i < npt; ++ i) {
        float angle = 2.0f * DELC_PI * (float)i / (float)npt;
        pt[i*2+0] = 1e10 * sinf(angle);
        pt[i*2+1] = 1e10 * cosf(angle);
    }
}

static inline void
gaussian(float *pt, size_t npt)
{
    for (size_t i = 0; i < npt; ++ i) {
        pt[i*2+0] = pseudoNormal() * 1e3;
        pt[i*2+1] = pseudoNormal() * 1e3;
    }
}

static inline void
grid(float *pt, size_t npt)
{
    size_t size = llroundf(floorf(sqrtf(npt)));
    for (size_t y = 0; y < size; ++ y)
    for (size_t x = 0; x < size; ++ x) {
        size_t i = y * size + x;
        pt[i*2+0] = x;
        pt[i*2+1] = y;
    }
    /* Don't leave uninitialized if unsquare */
    for (size_t i = size*size; i < npt; ++ i) {
        pt[i*2+0] = i;
        pt[i*2+1] = i;
    }
}

static inline void
uniform(float *pt, size_t npt)
{
    for (size_t i = 0; i < npt; ++ i) {
        pt[i*2+0] = frand() * 1e3;
        pt[i*2+1] = frand() * 1e3;
    }
}

