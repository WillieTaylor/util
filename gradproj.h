#ifndef __GRADPROJ_H__
#define __GRADPROJ_H__

/* ==== HEADER gradproj.h ==== */

/* Gradual Projection algorithm. */

/* ANSI C, IRIX 5.3, 21. Nov. 1995. */

/* ---- STANDARD HEADERS ---- */

/* #include <stdlib.h> */
#include <stdio.h>
#include <math.h>

/* ---- INCLUDE FILES ---- */

#include "util/aa/incl/matrix.h"	/* matrix handler */
#include "util/aa/incl/ql.h"		/* QL eigenvalue functions */

/* ---- PROTOTYPES ---- */

/* grad_proj(): performs a gradual projection iteration down to
 * 3D. The initial squared distances are in Dist, the corresponding
 * strictness values in Strict (both are triangular matrices with
 * size Size x Size). In each iteration, the distances in Dist
 * are applied to the actual distance matrix with appropriate
 * strictnesses. The final 3D embedding is obtained in Xyz
 * (should be an Rno x Rno square matrix).
 * Return value: the number of iterations done.
int drag_proj(const Trimat_ Dist, const Trimat_ Strict, int Rno, Sqmat_ Xyz);
 */
static double sq_dist(double Vec1[], double Vec2[], int Dim);

/* ==== END OF HEADER ==== */

#endif
