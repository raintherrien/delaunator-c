/*
 * Simple delaunay triangulation and voronoi partition of Poisson
 * distributed points.
 */

#include "delaunay/delaunator.h"
#include "delaunay/helper.h"
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define DELC_PI 3.14159265358979323846f

float    radius = 16.0f;
unsigned width  = 960;
unsigned height = 540;

/*
 * Returns a random floating point value between [0,1) (excluding 1)
 */
static inline float frand(void);

/*
 * Returns whether a 2D point (p) is within triangle (a,b,c)
 */
static int
pt_in_tri(const float * restrict p,
          const float * restrict a,
          const float * restrict b,
          const float * restrict c);

/*
 * Returns whether a 2D point (p) is within a polygon, defined by points
 * pol.
 */
static int
pt_in_pol(const float * restrict p,
          const float * restrict pol,
          size_t                 polsz);

/*
 * Calculates an RGB value from the first 12 bits of a size_t index.
 */
static void swizzle(size_t i, int *rgb);

/*
 * Rasterizes a two dimensional line, uses Bresenham's line algorithm.
 */
static void putline(int *rgb, unsigned width, unsigned height,
                    int r, int g, int b,
                    long x0, long y0, long x1, long y1);

/*
 * Rasterizes a triangle, suitable callback for foreach_tri()
 */
static void puttri(void *xrgb, size_t t, const float * restrict a,
                                         const float * restrict b,
                                         const float * restrict c);

/*
 * Rasterizes a voronoi cell, suitable callback for foreach_cell_dup()
 */
static void putcell(void *xrgb, size_t c, const float * restrict pt, size_t npt);

/*
 * Returns an array of two dimensional points pt (and length npt)
 * distributed within a domain defined by the dimensions (w,h), and
 * return zero on success. Sets and returns the value in errno on
 * failure, pt and npt are undefined.
 *
 * This is a crude implementation of:
 *  Robert Bridson, 2007. "Fast Poisson Disk Sampling in Arbitrary
 *  Dimensions." SIGGRAPH '07
 *
 * This implementation makes heap allocations and is non-deterministic.
 * Neither great qualities for serious applications.
 */
static int poisson(float **pt, size_t *ptsz, float r, float w, float h);

int
main(int argc, char **argv)
{
    int mode; /* 0 == triangulation, 1 == voronoi */
    if (argc < 2) goto print_help;
    else if (strcmp(argv[1], "-t") == 0) mode = 0;
    else if (strcmp(argv[1], "-v") == 0) mode = 1;
    else goto print_help;

    /* Seed our random number generator */
    srand((unsigned)time(NULL));

    /* Construct the Poisson distribution */
    float *pt = NULL;
    size_t ptsz = 0;
    if (poisson(&pt, &ptsz, radius, width, height) != 0) {
        perror("Error creating Poisson distribution");
        return EXIT_FAILURE;
    }

    /* Triangulate the distribution */
    size_t *del = calloc(DELAUNAY_SZ(ptsz), sizeof *del);
    if (del == NULL || triangulate(del, pt, ptsz) != 0) {
        perror("Error triangulating Poisson distribution");
        goto error_triangulating;
    }
    size_t *halfedge =  DELAUNAY_HALFEDGE(del, ptsz);
    size_t *triverts =  DELAUNAY_TRIVERTS(del, ptsz);
    size_t  ntrivert = *DELAUNAY_NTRIVERT(del, ptsz);

    /* Allocate an RGB buffer to draw our magic into */
    int *rgb = calloc(3 * width * height, sizeof *rgb);
    if (!rgb) {
        perror("Error creating PPM RGB buffer");
        goto error_calloc_rgb;
    }

    if (mode == 0) {
        foreach_tri(del, pt, ptsz, puttri, rgb);
    } else {
        foreach_cell_dup(del, pt, ptsz, putcell, rgb);
    }

    /* Draw points in white */
    for (size_t i = 0; i < ptsz&&0; ++ i) {
        float *p = pt + i*2;
        long x = lroundf(p[0]);
        long y = lroundf(p[1]);
        if (x >= 0 && x < width && y >= 0 && y < height) {
            int *c = rgb + 3*(y * width + x);
            c[0] = 255;
            c[1] = 255;
            c[2] = 255;
        }
    }

    /* Print a PPM image file to stdout */
    if (printf("P6 %d %d 255 ", width, height) < 0) goto printf_error;
    for (int y = height; y --; )
    for (int x = width;  x --; ) {
        size_t i = 3 * (y * width + x);
        int r = rgb[i+0];
        int g = rgb[i+1];
        int b = rgb[i+2];
        if (printf("%c%c%c", r, g, b) < 0) goto printf_error;
    }

    free(del);
    free(rgb);
    free(pt);
    return EXIT_SUCCESS;

printf_error:
    perror("Error writing to stdout");
    free(rgb);
error_calloc_rgb:
error_triangulating:
    free(del);
    free(pt);
    return EXIT_FAILURE;

print_help:
    fprintf(stderr, "Usage: ./poisson [-t|-v] > filename.ppm\n"
                    "\t-t  Rasterize Delaunay triangles\n"
                    "\t-v  Rasterize Voronoi diagram\n");
    return EXIT_SUCCESS;
}

static inline float
frand(void)
{
    return (float)rand() / (float)(RAND_MAX - 1);
}

static int
pt_in_tri(const float * restrict p,
          const float * restrict a,
          const float * restrict b,
          const float * restrict c)
{
#define SGN(A,B,C) (A[0]-C[0])*(B[1]-C[1])-(B[0]-C[0])*(A[1]-C[1])
    float d0 = SGN(p, a, b);
    float d1 = SGN(p, b, c);
    float d2 = SGN(p, c, a);
    return !( ((d0 < 0) || (d1 < 0) || (d2 < 0)) &&
              ((d0 > 0) || (d1 > 0) || (d2 > 0)) );
#undef  SGN
}

static int
pt_in_pol(const float * restrict p,
          const float * restrict pol,
          size_t                 polsz)
{
    /* I don't understand this function at all! */
    /* https://stackoverflow.com/questions/11716268/point-in-polygon-algorithm */
    int c = 0;
    size_t i = 0;
    size_t j = polsz - 1;
    for (; i < polsz; j = i ++) {
        if( ( (pol[i*2+1] >= p[1] ) != (pol[j*2+1] >= p[1]) ) &&
            (p[0] <= (pol[j*2+0] - pol[i*2+0]) * (p[1] - pol[i*2+1]) / (pol[j*2+1] - pol[i*2+1]) + pol[i*2+0])
          ) {
            c = !c;
        }
    }
  return c;
}

static void
swizzle(size_t i, int *rgb)
{
    /*
     * Very wonky RGB construction from the lower 12 bits of index
     */
    float r = (((unsigned)i & (1<< 0))>> 0) +
              (((unsigned)i & (1<< 3))>> 3) +
              (((unsigned)i & (1<< 6))>> 6) +
              (((unsigned)i & (1<< 9))>> 9);

    float g = (((unsigned)i & (1<< 1))>> 1) +
              (((unsigned)i & (1<< 4))>> 4) +
              (((unsigned)i & (1<< 7))>> 7) +
              (((unsigned)i & (1<<10))>>10);

    float b = (((unsigned)i & (1<< 2))>> 2) +
              (((unsigned)i & (1<< 5))>> 5) +
              (((unsigned)i & (1<< 8))>> 8) +
              (((unsigned)i & (1<<11))>>11);
    rgb[0] = (int)lroundf(255 * r / 4);
    rgb[1] = (int)lroundf(255 * g / 4);
    rgb[2] = (int)lroundf(255 * b / 4);
}

static void
putline(int *rgb, unsigned width, unsigned height, int r, int g, int b,
        long x0, long y0, long x1, long y1)
{
    long dx = +labs(x1 - x0);
    long dy = -labs(y1 - y0);
    long sx = x0 < x1 ? 1 : -1;
    long sy = y0 < y1 ? 1 : -1;
    long de = dx + dy;
    long e2 = 0; /* Accumulate error */

    for(;;) {
        if (x0 >= 0 && x0 < width && y0 >= 0 && y0 < height) {
            size_t i = 3 * (y0 * width + x0);
            rgb[i+0] = r;
            rgb[i+1] = g;
            rgb[i+2] = b;
        }
        if (x0 == x1 && y0 == y1) break;
        e2 = 2 * de;
        if (e2 >= dy) { de += dy; x0 += sx; }
        if (e2 <= dx) { de += dx; y0 += sy; }
    }
}

static void
puttri(void *xrgb, size_t t, const float * restrict a,
                             const float * restrict b,
                             const float * restrict c)
{
    int *rgb = xrgb;
    /* This could not be any lazier */
    for (size_t h = 0; h < height; ++ h)
    for (size_t w = 0; w < width; ++ w) {
        size_t i = 3 * (h * width + w);
        float p[2] = { w, h };
        if (pt_in_tri(p, a, b, c)) {
            swizzle(t, rgb + i);
        }
    }
}


static void
putcell(void *xrgb, size_t c, const float * restrict pt, size_t npt)
{
    int *rgb = xrgb;
    /* This could not be any lazier */
    for (size_t h = 0; h < height; ++ h)
    for (size_t w = 0; w < width; ++ w) {
        size_t i = 3 * (h * width + w);
        float p[2] = { w, h };
        if (pt_in_pol(p, pt, npt)) {
            swizzle(c, rgb + i);
        }
    }
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

