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

/*
 * ---------------------------------------------------------------------
 * The following are helper functions, descriptions and code copied from
 * https://mapbox.github.io/delaunator/
 * ---------------------------------------------------------------------
 */

/*
 * Triangle ids and half-edge ids are related.
 *  - The half-edges of triangle t are 3*t, 3*t + 1, and 3*t + 2.
 *  - The triangle of half-edge id e is floor(e/3).
 */

static inline void
triangle_edges(size_t t, size_t *halfedgeids)
{
    halfedgeids[0] = 3 * t;
    halfedgeids[1] = 3 * t + 1;
    halfedgeids[2] = 3 * t + 2;
}

static inline size_t
triangle_of_edge(size_t edge)
{
    return edge / 3;
}

static inline size_t
next_halfedge(size_t e)
{
    return (e % 3 == 2) ? e - 2 : e + 1;
}

static inline size_t
prev_halfedge(size_t e)
{
    return (e % 3 == 0) ? e + 2 : e - 1;
}

void triangle_center(size_t *triverts, float *pt, size_t t, float *q);

#endif /* DELAUNAY_HELPER_H_ */

