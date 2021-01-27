#include "delaunay/delaunator.h"
#include "delaunay/helper.h"

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/* Internal aliases for some type clarity */
typedef size_t vid;
typedef size_t tid;

/*
 * Creates a new triangle from three verts and links to neighbors
 */
static tid addtri(vid *triverts, size_t *ntrivert, tid *halfedge,
                  vid v0, vid v1, vid v2, tid t0, tid t1, tid t2);

/*
 * Returns whether the three points are counter-clockwise in order.
 */
static int ccw(float, float, float, float, float, float);

/*
 * Stores the center of the circumcircle of three points in c.
 */
static void circc(float, float, float, float, float, float, float *c);

/*
 * Returns the circumradius between three points, or FLT_MAX if points
 * are colinear.
 */
static float circr(float, float, float, float, float, float);

/*
 * Returns a pseudo-angle between two points which monotonically
 * increases with real angle but doesn't need expensive trigonometry.
 */
static size_t hashk(float, float, float, float, size_t mask);

/*
 * Returns whether the point p is within the circumcircle formed by the
 * three vertices.
 */
static int incircc(float ax, float ay, float bx, float by,
                   float cx, float cy, float px, float py);

/*
 * Recursively flips neighboring triangles until they (locally) satisfy
 * the Delaunay condition.
 *
 * Stack is a pre-allocated fixed sized stack used to avoid actual
 * recursion, the size of which must be at least npt.
 */
static tid legalize(tid *halfedge, vid *hullprev, tid *hulltris,
                    vid *triverts, size_t *hullstrt,
                    float *pt, size_t npt, tid t, tid *stack);

/*
 * Links two triangles with half edges
 */
static void link(tid *halfedge, tid, tid);

/*
 * Seed stores the vertices of a triangle close to the center of points
 * in the first three elements of s, and returns zero on success.
 *
 * If all points are colinear errno is set to EINVAL, which is returned,
 * and s is undefined.
 */
static int seed(float *pt, size_t npt, vid *s);

/*
 * Structure and comparator of a (index, distance) tuple used to order
 * points by distance from domain seed.
 */
struct pointdist { vid i; float d; };
static int pointdistcmp(const void *xa, const void *xb);

int
triangulate(size_t *delaunay, float *pt, size_t npt)
{
    if (npt < 3) {
        errno = ERANGE;
        return errno;
    }

    if (npt > DELAUNAY_MAXNPT) {
        errno = ERANGE;
        return errno;
    }

    /* Define pointers into buffer; see delaunator.h for layout */
    tid    *halfedge = DELAUNAY_HALFEDGE(delaunay, npt);
    vid    *hullhash = DELAUNAY_HULLHASH(delaunay, npt);
    vid    *hullnext = DELAUNAY_HULLNEXT(delaunay, npt);
    vid    *hullprev = DELAUNAY_HULLPREV(delaunay, npt);
    tid    *hulltris = DELAUNAY_HULLTRIS(delaunay, npt);
    vid    *triverts = DELAUNAY_TRIVERTS(delaunay, npt);
    size_t *ntrivert = DELAUNAY_NTRIVERT(delaunay, npt);
    size_t *hullsize = DELAUNAY_HULLSIZE(delaunay, npt);
    size_t *hullstrt = DELAUNAY_HULLSTRT(delaunay, npt);

    /*
     * Array of point indices to sort by distance from seed and a fixed
     * size stack for legalize.
     *
     * Static asserts here that: if the buffer allocation is safe from
     * overflow, then we're definitely safe with these temp allocations.
     */
    _Static_assert(sizeof(struct pointdist) <= 16, "");
    _Static_assert(sizeof(tid)              <= 16, "");
    struct pointdist *pds = malloc(npt * sizeof *pds);
    /* Fixed sized stack for legalize */
    tid *stack = malloc(npt * sizeof *stack);

    if (pds == NULL || stack == NULL)
    {
        goto free_buffers;
    }

#define ADDTRI(...) addtri(triverts, ntrivert, halfedge, __VA_ARGS__);
#define LEGALIZE(T) legalize(halfedge, hullprev, hulltris, triverts, \
                             hullstrt, pt, npt, T, stack);

    /* Assign initial values to SIZE_MAX, which means no vert/tri */
    memset(delaunay, 0xFF, DELAUNAY_SZ(npt) * sizeof *delaunay);

    /* Find seed */
    vid s[3] = { 0 };
    if (seed(pt, npt, s) != 0) {
        goto free_buffers;
    }

    float c[2];
    circc(pt[s[0]*2], pt[s[0]*2+1],
          pt[s[1]*2], pt[s[1]*2+1],
          pt[s[2]*2], pt[s[2]*2+1], c);

    for (size_t i = 0; i < npt; ++ i) {
        pds[i] = (struct pointdist) {
            .i = i,
            .d = hypotf(pt[i*2]-c[0], pt[i*2+1]-c[1])
        };
    }

    /* Sort indices by distance from "center" of pts */
    qsort(pds, npt, sizeof *pds, pointdistcmp);

    /* Create the initial hull around seed triangle */
    size_t hashsz = llroundf(ceilf(sqrtf(npt)));
    *hullstrt = s[0];
    *hullsize = 3;
    hullnext[s[0]] = hullprev[s[2]] = s[1];
    hullnext[s[1]] = hullprev[s[0]] = s[2];
    hullnext[s[2]] = hullprev[s[1]] = s[0];
    hulltris[s[0]] = 0;
    hulltris[s[1]] = 1;
    hulltris[s[2]] = 2;
    hullhash[hashk(c[0],c[1], pt[s[0]*2],pt[s[0]*2+1], hashsz)] = s[0];
    hullhash[hashk(c[0],c[1], pt[s[1]*2],pt[s[1]*2+1], hashsz)] = s[1];
    hullhash[hashk(c[0],c[1], pt[s[2]*2],pt[s[2]*2+1], hashsz)] = s[2];
    *ntrivert = 0;
    assert(*ntrivert + 3 < npt * 6);
    ADDTRI(s[0], s[1], s[2], SIZE_MAX, SIZE_MAX, SIZE_MAX);

    float p[2] = { FLT_MAX };
    for (size_t idx = 0; idx < npt; ++ idx) {
        vid i = pds[idx].i;
        float *v = pt + i*2;

        /* Skip near-duplicate points */
        if (fabsf(v[0]-p[0]) <= FLT_EPSILON &&
            fabsf(v[1]-p[1]) <= FLT_EPSILON)
        {
            continue;
        }
        p[0] = v[0];
        p[1] = v[1];

        /* Skip seed triangle points */
        if (i == s[0] || i == s[1] || i == s[2]) {
            continue;
        }

        /* Find a visible edge on the convex hull using edge hash */
        vid start = 0;
        size_t key = hashk(c[0],c[1], v[0],v[1], hashsz);
        for (size_t j = 0; j < hashsz; ++ j) {
            start = hullhash[(key + j) % hashsz];
            if (start != SIZE_MAX && start != hullnext[start]) {
                break;
            }
        }
        start = hullprev[start];
        vid e = start;
        vid q = hullnext[e];
        while (ccw(v[0],v[1], pt[e*2],pt[e*2+1], pt[q*2],pt[q*2+1])) {
            e = q;
            if (e == start) {
                e = SIZE_MAX;
                break;
            }
            q = hullnext[e];
        }
        /* Likely a near-duplicate point; skip it */
        if (e == SIZE_MAX) continue;

        /* Add the first triangle from the point */
        assert(*ntrivert + 3 < npt * 6);
        tid t = ADDTRI(e, i, hullnext[e], SIZE_MAX, SIZE_MAX, hulltris[e]);
        /* Flip triangles until they satisfy the Delaunay condition */
        hulltris[i] = LEGALIZE(t + 2);
        hulltris[e] = t; /* Keep track of boundary (hull) triangles */
        ++ hullsize;

        /* Walk through the hull adding more triangles and flipping */
        vid n = hullnext[e];
        q = hullnext[n];
        while (!ccw(v[0],v[1], pt[n*2],pt[n*2+1], pt[q*2],pt[q*2+1])) {
            assert(*ntrivert + 3 < npt * 6);
            t = ADDTRI(n, i, q, hulltris[i], SIZE_MAX, hulltris[n]);
            hulltris[i] = LEGALIZE(t + 2);
            hullnext[n] = n; /* mark as removed */
            -- hullsize;
            n = q;
            q = hullnext[n];
        }

        /* Walk backwards adding more tris and flipping */
        if (e == start) {
            q = hullprev[e];
            while (!ccw(v[0],v[1], pt[q*2],pt[q*2+1], pt[e*2],pt[e*2+1])) {
                assert(*ntrivert + 3 < npt * 6);
                t = ADDTRI(q, i, e, SIZE_MAX, hulltris[e], hulltris[q]);
                LEGALIZE(t + 2);
                hulltris[q] = t;
                hullnext[e] = e; /* mark as removed */
                -- hullsize;
                e = q;
                q = hullprev[e];
            }
        }

        /* Update the hull indices */
        *hullstrt   = hullprev[i] = e;
        hullnext[e] = hullprev[n] = i;
        hullnext[i] = n;

        /* Save the two new edges in the hash table */
        hullhash[hashk(c[0],c[1], v[0],v[1], hashsz)] = i;
        hullhash[hashk(c[0],c[1], pt[e*2],pt[e*2+1], hashsz)] = e;
    }

    /* XXX Shrink halfedges to only valid triangles */
    /* this.triangles = this._triangles.subarray(0, this.trianglesLen);
       this.halfedge = this._halfedges.subarray(0, this.trianglesLen); */

    free(stack);
    free(pds);

    return 0;

free_buffers:
    free(stack);
    free(pds);
    return errno;
#undef LEGALIZE
#undef ADDTRI
}

void
triangle_center(vid *triverts, float *pt, size_t t, float *q)
{
    size_t te[3];
    triangle_edges(t, te);
    vid tp[3] = { triverts[te[0]], triverts[te[1]], triverts[te[2]] };
    circc(pt[tp[0]*2],pt[tp[0]*2+1],
          pt[tp[1]*2],pt[tp[1]*2+1],
          pt[tp[2]*2],pt[tp[2]*2+1], q);
}

static int
seed(float *pt, size_t npt, vid *s)
{
    /* Calculate pt bounding box */
    float maxx = FLT_MIN;
    float maxy = FLT_MIN;
    float minx = FLT_MAX;
    float miny = FLT_MAX;
    for (vid i = 0; i < npt; ++ i) {
        float x = pt[i*2];
        float y = pt[i*2+1];
        if (x > maxx) maxx = x;
        if (x < minx) minx = x;
        if (y > maxy) maxy = y;
        if (y < miny) miny = y;
    }
    /* Pick a seed close to the center */
    float ctrx = (minx + maxx) / 2.0f;
    float ctry = (miny + maxy) / 2.0f;
    vid s0 = SIZE_MAX;
    vid s1 = SIZE_MAX;
    vid s2 = SIZE_MAX;
    float mind = FLT_MAX;
    for (vid i = 0; i < npt; ++ i) {
        float d = hypotf(pt[i*2]-ctrx, pt[i*2+1]-ctry);
        if (d < mind) {
            s0 = i;
            mind = d;
        }
    }
    /* Pick the Vertex closest to seed */
    mind = FLT_MAX;
    for (vid i = 0; i < npt; ++ i) {
        if (s0 == i) continue;
        float d = hypotf(pt[s0*2] - pt[ i*2], pt[s0*2+1] - pt[ i*2+1]);
        if (d < mind) {
            s1 = i;
            mind = d;
        }
    }
    /* Pick the third vertex which forms the smallest circumcircle */
    float minr = FLT_MAX;
    for (vid i = 0; i < npt; ++ i) {
        if (s0 == i || s1 == i) continue;
        float r = circr(pt[s0*2], pt[s0*2+1],
                        pt[s1*2], pt[s1*2+1],
                        pt[ i*2], pt[ i*2+1]);
        if (r < minr) {
            s2 = i;
            minr = r;
        }
    }

    /*
     * Note: original JS implementation handles colinear domains. I just
     * don't want the extra case.
     */
    if (minr == FLT_MAX) {
        errno = EINVAL;
        return errno;
    }

    /* Swap the order of our seed triangle until it is ccw */
    if (!ccw(pt[s0*2], pt[s0*2+1],
             pt[s1*2], pt[s1*2+1],
             pt[s2*2], pt[s2*2+1]))
    {
        vid t = s1;
        s1 = s2;
        s2 = t;
    }

    s[0] = s0;
    s[1] = s1;
    s[2] = s2;

    return 0;
}

static inline tid
addtri(vid *triverts, size_t *ntrivert, tid *halfedge,
       vid v0, vid v1, vid v2, tid t0, tid t1, tid t2)
{
    size_t t = *ntrivert;
    triverts[t+0] = v0;
    triverts[t+1] = v1;
    triverts[t+2] = v2;
    *ntrivert += 3;
    link(halfedge, t+0, t0);
    link(halfedge, t+1, t1);
    link(halfedge, t+2, t2);
    return t;
}

static inline int
ccw(float ax, float ay, float bx, float by, float cx, float cy)
{
    return (by - ay) * (cx - bx) - (bx - ax) * (cy - by) >= 0;
}

static inline void
circc(float ax, float ay, float bx, float by, float cx, float cy, float *c)
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

static inline float
circr(float ax, float ay, float bx, float by, float cx, float cy)
{
    float dx = bx - ax;
    float dy = by - ay;
    float ex = cx - ax;
    float ey = cy - ay;
    float bl = dx * dx + dy * dy;
    float cl = ex * ex + ey * ey;
    float d  = dx * ey - dy * ex;

    float x = (ey * bl - dy * cl) * 0.5f / d;
    float y = (dx * cl - ex * bl) * 0.5f / d;

    if (bl != 0.0f && cl != 0.0f && d  != 0.0f)
    {
        return x * x + y * y;
    } else {
        return FLT_MAX;
    }
}

static inline size_t
hashk(float cx, float cy, float vx, float vy, size_t mask)
{
    float d0 = vx - cx;
    float d1 = vy - cy;
    /*
     * Pseudo angle monotonically increases with real angle, but doesn't
     * need expensive trigonometry
     */
    float p = d0 / (fabsf(d0) + fabsf(d1));
    float angle = (d1 > 0 ? 3 - p : 1 + p) / 4.0f; /* [0..1] */
    return llroundf(floorf(angle * (float)mask)) % mask;
}

static inline int
incircc(float ax, float ay, float bx, float by,
        float cx, float cy, float px, float py)
{
    float dx = ax - px;
    float dy = ay - py;
    float ex = bx - px;
    float ey = by - py;
    float fx = cx - px;
    float fy = cy - py;

    float ap = dx * dx + dy * dy;
    float bp = ex * ex + ey * ey;
    float cp = fx * fx + fy * fy;

    return ( dx * (ey * cp - bp * fy) -
             dy * (ex * cp - bp * fx) +
             ap * (ex * fy - ey * fx) ) < 0;
}

static tid
legalize(tid *halfedge, vid *hullprev, tid *hulltris, vid *triverts,
         size_t *hullstrt, float *pt, size_t npt, tid a, tid *stack)
{
    size_t stacksz = 0;

    tid ar = 0;

    /* Recursion eliminated with a fixed-size stack. */
    for (;;) {
        (void) npt;
        assert(stacksz < npt);

        tid b = halfedge[a];

        /*
         * If the pair of triangles doesn't satisfy the Delaunay
         * condition (p1 is inside the circumcircle of [p0, pl, pr]),
         * flip them, then do the same check/flip recursively for the
         * new pair of triangles.
         *
         *           pl                    pl
         *          /||\                  /  \
         *       al/ || \bl            al/    \a
         *        /  ||  \              /      \
         *       /  a||b  \    flip    /___ar___\
         *     p0\   ||   /p1   =>   p0\---bl---/p1
         *        \  ||  /              \      /
         *       ar\ || /br             b\    /br
         *          \||/                  \  /
         *           pr                    pr
         *
         */

        tid a0 = a - a % 3;
        ar = a0 + (a + 2) % 3;

        /* Convex hull edge */
        if (b == SIZE_MAX) {
            if (stacksz == 0) break;
            a = stack[-- stacksz];
            continue;
        }

        tid b0 = b - b % 3;
        tid al = a0 + (a + 1) % 3;
        tid bl = b0 + (b + 2) % 3;

        vid p0 = triverts[ar];
        vid pr = triverts[a];
        vid pl = triverts[al];
        vid p1 = triverts[bl];

        int flip = incircc(pt[p0*2], pt[p0*2+1], pt[pr*2], pt[pr*2+1],
                           pt[pl*2], pt[pl*2+1], pt[p1*2], pt[p1*2+1]);
        if (flip) {
            triverts[a] = p1;
            triverts[b] = p0;

            tid hbl = halfedge[bl];

            /*
             * Edge swapped on the other side of the hull (rare); fix
             * the halfedge reference
             */
            if (hbl == SIZE_MAX) {
                vid e = *hullstrt;
                do {
                    if (hulltris[e] == bl) {
                        hulltris[e] = a;
                        break;
                    }
                    e = hullprev[e];
                } while (e != *hullstrt);
            }
            link(halfedge, a, hbl);
            link(halfedge, b, halfedge[ar]);
            link(halfedge, ar, bl);

            tid br = b0 + (b + 1) % 3;

            stack[stacksz ++] = br;
        } else {
            if (stacksz == 0) break;
            a = stack[-- stacksz];
        }
    }

    return ar;
}

static inline void
link(tid *halfedge, tid a, tid b)
{
    halfedge[a] = b;
    if (b != SIZE_MAX) halfedge[b] = a;
}

static inline int
pointdistcmp(const void *xa, const void *xb)
{
    const struct pointdist *a = xa;
    const struct pointdist *b = xb;
    if (a->d < b->d) return -1;
    if (a->d > b->d) return 1;
    return 0;
}

