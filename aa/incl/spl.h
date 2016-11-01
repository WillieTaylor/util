#ifndef __SPL_H__
#define __SPL_H__

/* ==== HEADER spl.h ==== */

/* For the fit and evaluation of third-order splines.
 * This routine collection is different from its predecessor
 * "spline.c" in that it is more 'object-oriented'. A spline
 * structure Splrec_ holds the necessary auxiliary arrays
 * and evaluation functions operate on this record.
 */

/* Output from p2c, the Pascal-to-C translator */
/* From input file "spline.p" */

/* ANSI C, IRIX 4.0.5, 1. Oct. 1993. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

/* ---- GLOBAL TYPES ---- */

/* Splrec_ : data structure for spline variables. X[] and Y[] are
 * the dependent and independent variables. Y2[] holds the second
 * derivatives, Yin[] holds the integrals. All four arrays are
 * N long.
 */
typedef struct
{
    double *X, *Y, *Y2, *Yin;
    long N;
}
Splrec_ ;

/* ---- PROTOTYPES ---- */

int init_spl(long n, Splrec_ *Spl);
int fit_spl(double yp1, double ypn, Splrec_ *Spl);
double eval_spl(const Splrec_ *Spl, 
		double xi, double *Der1, double *Der2,
		double *Der3, double *Integ);
double integ_spl(const Splrec_ *Spl, double Low, double Up);
void remove_spl(Splrec_ *Spl);

/* ==== END OF HEADER spline.h ==== */
#endif
