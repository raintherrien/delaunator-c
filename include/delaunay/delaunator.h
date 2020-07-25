/*
 * Fast 2D Delaunay triangulation in C. Yet another port of Delaunator.
 * https://github.com/mapbox/delaunator
 *
 * Copyright (c) 2017, Mapbox
 * Copyright (c) 2020, Rain Therrien
 * Distributed under the ISC license.
 * See accompanying LICENSE
 */

#ifndef DELAUNAY_DELAUNATOR_H_
#define DELAUNAY_DELAUNATOR_H_

#include <stddef.h>

/*
 * Layout depends on point count and shouldn't matter externally anyway.
 *
 * size_t halfedge[npt * 6]; // Edge to adjacent triangle
 * size_t hullhash[npt];     // Hull verts in order of pseudo-angle
 * size_t hullnext[npt];     // Edge to next edge
 * size_t hullprev[npt];     // Edge to prev edge
 * size_t hulltris[npt];     // Edge to adjacent triangle
 * size_t triverts[npt * 6]; // Triangle vertices
 * size_t ntrivert;
 * size_t hullsize;
 * size_t hullstrt;
 *
 * Total size: sizeof(size_t) * (npt * 16 + 3)
 *
 * helper.h defines macros to retrieve these pointers.
 */
typedef size_t *delaunay;

/*
 * Calculates the delaunay triangulation of points and returns zero on
 * success, otherwise sets and returns the value in errno.
 *
 * Internally performs allocation of delaunay array. delaunay could be
 * passed as a pre-allocated array, but the user would need a function
 * to determine the required size and that is probably unecessary.
 */
int triangulate(delaunay *, float *pt, size_t npt);

/*
 * Frees data allocated by triangulate and destroys a delaunay state.
 */
void delaunay_free(delaunay *);

#endif /* DELAUNAY_DELAUNATOR_H_ */
