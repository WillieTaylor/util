#ifndef __SVD_H__
#define __SVD_H__

/* ==== HEADER svd.h ==== */

/* Singular Value Decomposition, based on the Numerical Recipes
 * routine 'svdcmp'. The following modifications are made:-
 *   + double precision throughout
 *   + elimination of macros
 *   + ANSI C prototyping
 *   + a wrapper routine to implement [0..N-1] style indexing
 *   + elimination of the functions vector(), nrerror() and free_vector()
 * 
 * Reference:- Press, W.H., Flannery, B.P., Teukolsky, S.A. and Vetterling, W.T.:
 * Numerical Recipes in C: The Art of Scientific Computing.
 * Cambridge University Press, 1992.
 */

/* ANSI C, IRIX 5.2, 26. July 1994. Modifications by Andras Aszodi. */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

/* ---- DEFINITIONS ---- */

#ifdef DBL_EPSILON
#define SVD_EPSILON (10.0*DBL_EPSILON)
#else
#define SVD_EPSILON (1.0e-10)
#endif

/* ---- PROTOTYPES ---- */

int svd_setup(int Row, int Col, double ***U, double **W, double ***V);
int svd_comp(const double **A, int Row, int Col, 
	double **U, double *W, double **V);
int rank_cond(double W[], int N, double Eps, double *Cond);
void svd_solve(const double **U, const double W[], const double **V, 
	int Row, int Col, const double B[], double X[]);
void free_svd(int Row, int Col, double **U, double *W, double **V);

/* ==== END OF HEADER svd.h ==== */
#endif
