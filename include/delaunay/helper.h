/*
 * Fast 2D Delaunay triangulation in C. Yet another port of Delaunator.
 * https://github.com/mapbox/delaunator
 *
 * Copyright (c) 2017, Mapbox
 * Copyright (c) 2021, Rain Therrien
 * Distributed under the ISC license.
 * See accompanying LICENSE
 */

#ifndef DELAUNAY_HELPER_H_
#define DELAUNAY_HELPER_H_

#include <stddef.h>

/* Pointers into buffer; see delaunator.h for layout */
#define DELAUNAY_HALFEDGE(D,NPT) ((D))
#define DELAUNAY_HULLHASH(D,NPT) ((D) + (size_t)(NPT) * 6)
#define DELAUNAY_HULLNEXT(D,NPT) ((D) + (size_t)(NPT) * 7)
#define DELAUNAY_HULLPREV(D,NPT) ((D) + (size_t)(NPT) * 8)
#define DELAUNAY_HULLTRIS(D,NPT) ((D) + (size_t)(NPT) * 9)
#define DELAUNAY_TRIVERTS(D,NPT) ((D) + (size_t)(NPT) * 10)
#define DELAUNAY_NTRIVERT(D,NPT) ((D) + (size_t)(NPT) * 16)
#define DELAUNAY_HULLSIZE(D,NPT) ((D) + (size_t)(NPT) * 16 + 1)
#define DELAUNAY_HULLSTRT(D,NPT) ((D) + (size_t)(NPT) * 16 + 2)
/* Max number of points before overflow */
#define DELAUNAY_MAXNPT          ((SIZE_MAX - 3) / 16 / sizeof(uint32_t))

static inline void
circumcenter(float ax, float ay, float bx, float by, float cx, float cy, float *c)
{
    float dx = bx - ax;
    float dy = by - ay;
    float ex = cx - ax;
    float ey = cy - ay;

    float bl = dx * dx + dy * dy;
    float cl = ex * ex + ey * ey;
    float d  = 0.5f / (dx * ey - dy * ex);

    c[0] = ax + (ey * bl - dy * cl) * d;
    c[1] = ay + (dx * cl - ex * bl) * d;
}

/*
 * The following are helper functions, descriptions and code copied from
 * https://mapbox.github.io/delaunator/
 * I've also reproduced the descriptions there for clarity.
 */

/*
 * tri_edges() stores the half-edge indices of triangle t in edges. The
 * half-edges of triangle t are 3*t, 3*t+1, and 3*t+2.
 */
static inline void
tri_edges(uint32_t t, uint32_t edges[3])
{
    edges[0] = 3 * t;
    edges[1] = 3 * t + 1;
    edges[2] = 3 * t + 2;
}

/*
 * tri_of_edges() returns the triangle id of the half-edge e. The
 * triangle of half-edge e is floor(e/3).
 */
static inline uint32_t
tri_of_edge(uint32_t e)
{
    return e / 3;
}

/*
 * tri_next_edge() and tri_prev_edge() return the next or previous edge
 * in a triangle.
 */
static inline uint32_t
tri_next_edge(uint32_t e)
{
    return (e % 3 == 2) ? e - 2 : e + 1;
}
static inline uint32_t
tri_prev_edge(uint32_t e)
{
    return (e % 3 == 0) ? e + 2 : e - 1;
}

static inline void
tri_center(uint32_t *triverts, const float * restrict pt, uint32_t t, float *q)
{
    uint32_t te[3];
    tri_edges(t, te);
    uint32_t tp[3] = { triverts[te[0]], triverts[te[1]], triverts[te[2]] };
    circumcenter(pt[tp[0]*2],pt[tp[0]*2+1],
                 pt[tp[1]*2],pt[tp[1]*2+1],
                 pt[tp[2]*2],pt[tp[2]*2+1], q);
}

/*
 * foreach_tri() and foreach_cell_dup() are more-or-less direct ports
 * from the Delaunator guide, but with even less of a bent towards
 * efficiency. foreach_cell_dup() specifically does not check for
 * duplicate cells, so that is left to the client (which they can do by
 * bitmasking the cell index c). These functions are more informational
 * than useful.
 *
 * Oh, and foreach_cell_dup might trigger a stack overflow on you if a
 * cell has too many hull vertices!
 */
static inline void
foreach_tri(uint32_t *delaunay, const float *pt, uint32_t npt,
            void(*fn)(void *, uint32_t t, const float * restrict a,
                                          const float * restrict b,
                                          const float * restrict c),
            void *arg)
{
    uint32_t *triverts =  DELAUNAY_TRIVERTS(delaunay, npt);
    uint32_t ntris     = *DELAUNAY_NTRIVERT(delaunay, npt) / 3;
    uint32_t e[3];
    for (uint32_t t = 0; t < ntris; ++ t) {
        tri_edges(t, e);
        fn(arg, t, pt+2*triverts[e[0]],
                   pt+2*triverts[e[1]],
                   pt+2*triverts[e[2]]);
    }
}

static inline void
foreach_cell_dup(uint32_t *delaunay, const float *pt, uint32_t npt,
                 void(*fn)(void *, uint32_t c, const float * restrict pt, uint32_t npt),
                 void *arg)
{
    uint32_t *halfedge =  DELAUNAY_HALFEDGE(delaunay, npt);
    uint32_t *triverts =  DELAUNAY_TRIVERTS(delaunay, npt);
    uint32_t ntriverts = *DELAUNAY_NTRIVERT(delaunay, npt);
    for (uint32_t incoming = 0; incoming < ntriverts; ++ incoming) {
        uint32_t endpoint = triverts[tri_next_edge(incoming)];
        if (halfedge[endpoint] == -1) continue;

        /* Count the number of Voronoi cell vertices */
        uint32_t cell_npt = 0;
        uint32_t current = incoming;
        do {
            ++ cell_npt;
            current = halfedge[tri_next_edge(current)];
        } while (current != -1 && current != incoming);

        /* Determine cell vertices */
        /* XXX stack overflow right here :) */
        float cell_pt[cell_npt * 2];
        current = incoming;
        for (uint32_t i = 0; i < cell_npt; ++ i) {
            tri_center(triverts, pt, tri_of_edge(current), cell_pt+i*2);
            current = halfedge[tri_next_edge(current)];
        }

        fn(arg, endpoint, cell_pt, cell_npt);
    }
}

#endif /* DELAUNAY_HELPER_H_ */

