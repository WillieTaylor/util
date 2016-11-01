#ifndef __PARAMEST_H__
#define __PARAMEST_H__

/* ==== HEADER paramest.h ==== */

/* Output from p2c, the Pascal-to-C translator */
/* From input file "regr.p" (translated 7-Sep-93) */

/* This module is based on my parameter estimation library
 * written in Hungary in 1990. The algorithms come from
 * Valko+Vajda. The p2c raw output is after-edited.
 */

/* ANSI C, IRIX 4.0.5, 14. Sep. 1993. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ---- INCLUDE FILES ---- */

#include "matrix.h"

/* ---- GLOBAL TYPES ---- */

/* Userfunct_ : the type of the function which describes
 * the (X, P)->Y relationship for the nonlinear regression
 * routine as (X, P, Y). X is the vector of independent variables, P is
 * the vector of parameters and Y is the vector of the 
 * dependent variables. This general forms enables the
 * parameter estimation for vector-vector functions.
 */
typedef void (*Userfunct_)(const double *, const double *, double *);

/* Userfunct11_ : the type of the function which describes
 * the (X, P)->Y relationship for the 1:1 nonlinear regression
 * routine where X, Y are scalars,  and P is a vector of parameters.
 */
typedef double (*Userfunct11_)(double, const double *);

/* ---- PROTOTYPES ---- */

long lin_reg(long Nx, long Nm, double *Xmeas[], double Ymeas[], double W[], 
	    double *Sres, double P[], double Sdev[], Trimat_ Correl,
	    double *Dstat);
double lms_fit(double *x, double *y, long N, double *b, double *a, double *Sigb,
	    double *Siga, double *Sigab);
long nonlin_reg(Userfunct_ Funct, long Nx, long Ny, long Np, long Nm,
	       long ItMax, double StepLim, double **Xmeas, double **Ymeas,
	       double **W, long *ItNo, double *Sres, double P[],
	       double Sdev[], Trimat_ Correl);
long nonlin11_reg(Userfunct11_ Funct, long Np, long Nm,
	       long ItMax, double StepLim, double *Xmeas, double *Ymeas,
	       double *W, long *ItNo, double *Sres, double P[],
	       double Sdev[], Trimat_ Correl);
double tcrit_95(long Nf);

/* ==== END OF HEADER paramest.h ==== */
#endif
