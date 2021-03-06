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

/* Pointers into buffer; see delaunator.h for layout */
#define DELAUNAY_HALFEDGE(D,NPT) ((D))
#define DELAUNAY_HULLHASH(D,NPT) ((D) + (NPT) * 6)
#define DELAUNAY_HULLNEXT(D,NPT) ((D) + (NPT) * 7)
#define DELAUNAY_HULLPREV(D,NPT) ((D) + (NPT) * 8)
#define DELAUNAY_HULLTRIS(D,NPT) ((D) + (NPT) * 9)
#define DELAUNAY_TRIVERTS(D,NPT) ((D) + (NPT) * 10)
#define DELAUNAY_NTRIVERT(D,NPT) ((D) + (NPT) * 16)
#define DELAUNAY_HULLSIZE(D,NPT) ((D) + (NPT) * 16 + 1)
#define DELAUNAY_HULLSTRT(D,NPT) ((D) + (NPT) * 16 + 2)
/* Max number of points before overflow */
#define DELAUNAY_MAXNPT          ((SIZE_MAX - 3) / 16 / sizeof(size_t))

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
tri_edges(size_t t, size_t edges[3])
{
    edges[0] = 3 * t;
    edges[1] = 3 * t + 1;
    edges[2] = 3 * t + 2;
}

/*
 * tri_of_edges() returns the triangle id of the half-edge e. The
 * triangle of half-edge e is floor(e/3).
 */
static inline size_t
tri_of_edge(size_t e)
{
    return e / 3;
}

/*
 * tri_next_edge() and tri_prev_edge() return the next or previous edge
 * in a triangle.
 */
static inline size_t
tri_next_edge(size_t e)
{
    return (e % 3 == 2) ? e - 2 : e + 1;
}
static inline size_t
tri_prev_edge(size_t e)
{
    return (e % 3 == 0) ? e + 2 : e - 1;
}

static inline void
tri_center(size_t *triverts, const float * restrict pt, size_t t, float *q)
{
    size_t te[3];
    tri_edges(t, te);
    size_t tp[3] = { triverts[te[0]], triverts[te[1]], triverts[te[2]] };
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
foreach_tri(size_t *delaunay, const float *pt, size_t npt,
            void(*fn)(void *, size_t t, const float * restrict a,
                                        const float * restrict b,
                                        const float * restrict c),
            void *arg)
{
    size_t *triverts =  DELAUNAY_TRIVERTS(delaunay, npt);
    size_t ntris     = *DELAUNAY_NTRIVERT(delaunay, npt) / 3;
    size_t e[3];
    for (size_t t = 0; t < ntris; ++ t) {
        tri_edges(t, e);
        fn(arg, t, pt+2*triverts[e[0]],
                   pt+2*triverts[e[1]],
                   pt+2*triverts[e[2]]);
    }
}

static inline void
foreach_cell_dup(size_t *delaunay, const float *pt, size_t npt,
                 void(*fn)(void *, size_t c, const float * restrict pt, size_t npt),
                 void *arg)
{
    size_t *halfedge =  DELAUNAY_HALFEDGE(delaunay, npt);
    size_t *triverts =  DELAUNAY_TRIVERTS(delaunay, npt);
    size_t ntriverts = *DELAUNAY_NTRIVERT(delaunay, npt);
    for (size_t incoming = 0; incoming < ntriverts; ++ incoming) {
        size_t endpoint = triverts[tri_next_edge(incoming)];
        if (halfedge[endpoint] == -1) continue;

        /* Count the number of Voronoi cell vertices */
        size_t cell_npt = 0;
        size_t current = incoming;
        do {
            ++ cell_npt;
            current = halfedge[tri_next_edge(current)];
        } while (current != -1 && current != incoming);

        /* Determine cell vertices */
        /* XXX stack overflow right here :) */
        float cell_pt[cell_npt * 2];
        current = incoming;
        for (size_t i = 0; i < cell_npt; ++ i) {
            tri_center(triverts, pt, tri_of_edge(current), cell_pt+i*2);
            current = halfedge[tri_next_edge(current)];
        }

        fn(arg, endpoint, cell_pt, cell_npt);
    }
}

#endif /* DELAUNAY_HELPER_H_ */

