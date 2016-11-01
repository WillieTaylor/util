#ifndef __SPLINE_H__
#define __SPLINE_H__

/* ==== HEADER spline.h ==== */

/* For the fit and evaluation of third-order splines. */

/* Output from p2c, the Pascal-to-C translator */
/* From input file "spline.p" */

/* ANSI C, IRIX 4.0.5, 20. Sep. 1993. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <stdlib.h>

/* ---- PROTOTYPES ---- */

void spline_fitter(long n, double *x, double *y, double *y2, double *yin,
		  double yp1, double ypn);
double spline_eval(long n, double *x, double *y, double *y2, double *yin,
		double xi, double *Der1, double *Der2,
		double *Der3, double *Integ);
long pos_finder(long n, double *x, double xi, double *a, double *h);
double spline_intpol(long n, long kl, double a, double h, double *y,
		    double *y2);
double spline_der1(long n, long kl, double a, double h, double *y, double *y2);
double spline_der2(long n, long kl, double a, double h, double *y2);
double spline_der3(long n, long kl, double h, double *y2);
double spline_integ(long n, double *x, double *y, double *y2, double *yin,
		   double Low, double Up);

/* ==== END OF HEADER spline.h ==== */
#endif
