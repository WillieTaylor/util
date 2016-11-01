#ifndef __CDF_H__
#define __CDF_H__

/* ==== HEADER cdf.h ==== */

/* Routine collection for estimating cumulative distribution functions. */

/* ANSI C, IRIX 4.0.5, 1. Oct. 1993. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* ---- GLOBAL TYPES ---- */

/* Cdfrec_ : this structure holds the variables for estimating a CDF.
 * X is the indep. var, Y is the CDF, Yint holds the integer counts
 * of the values falling below a given X. X, Y, Yint are Bno long arrays.
 * N is the total number of values processed.
 */
typedef struct
{
    double *X, *Y;
    int *Yint;
    int N, Bno;
}
Cdfrec_ ;

/* ---- PROTOTYPES ---- */

int init_cdf(int Binno, double Low, double Up, Cdfrec_ *Cdf);
int add_cdf(double Value, Cdfrec_ *Cdf);
int eval_cdf(Cdfrec_ *Cdf);
void remove_cdf(Cdfrec_ *Cdf);

/* ==== END OF HEADER cdf.h ==== */
#endif
