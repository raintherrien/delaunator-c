/*
 * Simple delaunay triangulation and voronoi partition of Poisson
 * distributed points.
 * TODO: Actual triangulation!
 */

#include "delaunay/delaunator.h"
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define DELC_PI 3.14159265358979323846f

/*
 * Returns an array of two dimensional points pt (and length npt)
 * distributed within a domain defined by the dimensions (w,h), and
 * return zero on success. Returns errno on failure, pt and npt are
 * undefined.
 *
 * This is a crude implementation of:
 *  Robert Bridson, 2007. "Fast Poisson Disk Sampling in Arbitrary
 *  Dimensions." SIGGRAPH '07
 *
 * This implementation makes heap allocations and is non-deterministic.
 * Neither great qualities for serious applications.
 */
static int poisson(float **pt, size_t *ptsz, float r, float w, float h);

/*
 * Returns a random floating point value between [0,1) (excluding 1)
 */
static inline float frand(void);

int
main(void)
{
    float    radius = 16.0f;
    unsigned width  = 960;
    unsigned height = 540;

    srand((unsigned)time(NULL));

    float *pt = NULL;
    size_t ptsz = 0;
    int rc = poisson(&pt, &ptsz, radius, width, height);
    if (rc != 0) {
        errno = rc;
        perror("Error creating Poisson distribution");
        return EXIT_FAILURE;
    }

    int *rgb = calloc(3 * width * height, sizeof *rgb);
    if (!rgb) {
        perror("Error creating PPM RGB buffer");
        goto free_poisson;
    }

    for (size_t i = 0; i < ptsz; ++ i) {
        float *p = pt + i*2;
        long x = lroundf(p[0]);
        long y = lroundf(p[1]);
        assert(x >= 0 && x < width);
        assert(y >= 0 && y < height);
        int *c = rgb + 3*(y * width + x);
        c[0] = 255;
        c[1] = 255;
        c[2] = 255;
    }

    if (printf("P6 %d %d 255 ", width, height) < 0) goto printf_error;
    for (int y = height; y --; )
    for (int x = width;  x --; ) {
        size_t i = 3 * (y * width + x);
        int r = rgb[i+0];
        int g = rgb[i+1];
        int b = rgb[i+2];
        if (printf("%c%c%c", r, g, b) < 0) goto printf_error;
    }

    free(rgb);
    free(pt);
    return EXIT_SUCCESS;

printf_error:
    perror("Error writing to stdout");
    free(rgb);
free_poisson:
    free(pt);
    return EXIT_FAILURE;
}

static int
poisson(float **pt, size_t *ptsz, float r, float w, float h)
{
    /*
     * Backing grid to accelerate spacial searches stores the index into
     * (*pt) of the point at that location (position floored).
     *
     * Backing grid cells have length radius/sqrt(2) such that two
     * points farther than radius apart will never occupy the same cell.
     */
    float cl = r / 1.414213562f;
    long gw = lroundf(ceilf(w / cl));
    long gh = lroundf(ceilf(h / cl));

    /* Buffer for points, with reasonable maximum count of area */
    size_t mxnpt = gw * gh;
    *pt = malloc(2 * mxnpt * sizeof **pt);
    size_t npt = 0;

    /* Index into pt, queried with grid coordinate */
    size_t *lookup = malloc(gw * gh * sizeof *lookup);

    /* Active list of points that have not yet spawned children */
    size_t *active = malloc(gw * gh * sizeof *active);
    size_t nactive = 0;

    if (!*pt || !lookup || !active) {
        perror("Error allocating buffers");
        goto free_buffers;
    }

    /* Assign initial values to SIZE_MAX, which means no point */
    memset(lookup, 0xFF, gw * gh * sizeof *lookup);

    /* Prime with one random point */
    {
        float strtx = frand() * w;
        float strty = frand() * h;
        long gridx = lroundf(floorf(strtx / cl));
        long gridy = lroundf(floorf(strty / cl));
        (*pt)[0] = strtx;
        (*pt)[1] = strty;
        lookup[gridy * gw + gridx] = 0;
        active[0] = 0;
        npt = 1;
        nactive = 1;
    }

    do {
        float *parent = (*pt) + 2*active[-- nactive];

        /* Maximum attempts made to spawn a child point */
        size_t max_child_samples = 30;
        for (size_t i = 0; i < max_child_samples && npt < mxnpt; ++ i) {
            /* Calculate a new point between radius and 2*radius */
            float dist = r * (1 + frand());
            float angle = DELC_PI * 2.0f * frand();
            float childx = parent[0] + cosf(angle) * dist;
            float childy = parent[1] + sinf(angle) * dist;
            long gridx = lroundf(floorf(childx / cl));
            long gridy = lroundf(floorf(childy / cl));

            /* Check for points within radius, out of bounds */
            if (gridx < 2 || gridx > gw-3 ||
                gridy < 2 || gridy > gh-3) {
                continue; /* Neighbors would be out of bounds */
            }
            int gridi = gridy * gw + gridx;
            int intersect = 0;
            /* TODO: Doesn't need to check corner neighbors */
            for (int ny = -2; ny <= 2 && !intersect; ++ ny)
            for (int nx = -2; nx <= 2 && !intersect; ++ nx) {
                size_t ni = lookup[gridi + nx + ny * gw];
                if (ni == SIZE_MAX) continue;
                float *n = (*pt) + 2*ni;
                intersect = hypotf(childx - n[0], childy - n[1]) < r;
            }

            /* This point is within radius of a neighbor */
            if (intersect) continue;

            /* This is a valid point */
            (*pt)[npt*2+0] = childx;
            (*pt)[npt*2+1] = childy;
            active[nactive ++] = npt;
            lookup[gridi] = npt;
            ++ npt;
        }
    } while (nactive && npt < mxnpt);

    /* Reclaimed unused points */
    float *newpt = realloc(*pt, 2 * npt * sizeof **pt);
    if (!newpt) {
        perror("Error shrinking Poisson points");
        goto free_buffers;
    }
    *pt = newpt;
    *ptsz = npt;

    free(active);
    free(lookup);
    return 0;

free_buffers:
    free(active);
    free(lookup);
    free(*pt);
    return errno;
}

static inline float
frand(void)
{
    return (float)rand() / (float)(RAND_MAX - 1);
}
