/*
 * Fast 2D Delaunay triangulation in C. Yet another port of Delaunator.
 * https://github.com/mapbox/delaunator
 *
 * Copyright (c) 2017, Mapbox
 * Copyright (c) 2021, Rain Therrien
 * Distributed under the ISC license.
 * See accompanying LICENSE
 */

#ifndef DELAUNAY_DELAUNAY_H_
#define DELAUNAY_DELAUNAY_H_

#include <stddef.h>

/*
 * triangulate() calculates the delaunay triangulation of points (pt).
 * npt is the length or number of points, represented by two contiguous
 * floating point coordinates, in pt. The results are stored in
 * delaunay, which helper.h provides macros for accessing.
 * Zero is returned on success, otherwise the state of delaunay is
 * unspecified and errno is returned as:
 * EINVAL shall be returned if all points are colinear;
 * ERANGE shall be returned if less than three points are passed, or if
 * the number of points would cause a numeric overflow;
 * ENOMEM shall be returned if triangulate is unable to allocate scrach
 * buffers (See TODO below).
 *
 * delaunay array layout depends on point count:
 *   size_t halfedge[npt * 6]; // Edge to adjacent triangle
 *   size_t hullhash[npt];     // Hull verts in order of pseudo-angle
 *   size_t hullnext[npt];     // Edge to next edge
 *   size_t hullprev[npt];     // Edge to prev edge
 *   size_t hulltris[npt];     // Edge to adjacent triangle
 *   size_t triverts[npt * 6]; // Triangle vertices
 *   size_t ntrivert;
 *   size_t hullsize;
 *   size_t hullstrt;
 * helper.h defines macros to retrieve these pointers.
 * Total size: sizeof(size_t) * (npt * 16 + 3)
 *
 * Use the DELAUNAY_SZ macro to allocate this array.
 *
 * TODO: This algorithm DOES perform two allocations which can fail:
 * one for a fixed sized flip stack and one for a point distance field.
 * Both of these could be hoisted out into a "scratch" buffer allocated
 * by the user.
 */

#define DELAUNAY_SZ(NPT) ((NPT) * 16 + 3)

int triangulate(size_t *delaunay, float *pt, size_t npt);

#endif /* DELAUNAY_DELAUNAY_H_ */

