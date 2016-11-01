#ifndef __LRSVD_H__
#define __LRSVD_H__

/* ==== HEADER lrsvd.h ==== */

/* Weighted linear regression routines based on SVD. */

/* ANSI C, IRIX 4.0.5, 8. July 1994. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* ---- INCLUDE FILES ---- */

#include "matrix.h"
#include "svd.h"

/* ---- TYPEDEFS ---- */

/* Xfunc1_: a pointer to a function that transforms a
 * double into another. Useful in linear regression problems
 * with indep. variable transformations:
 * y(x)=a1*f1(x)+a2*f2(x)+...+aN*fN(x)
 * where x is scalar. Cf. lr_xtr().
 */
typedef double (*Xfunc1_)(double);

/* Xfuncn_: a pointer to a function that transforms a
 * double vector into a double scalar. To be used in linear
 * regression problems of the form
 * y(x1,...,xq)=a1*f1(x1,...,xq)+a2*f2(x1,...,xq)+...+aN*fN(x1,...,xq)
 * The first parameter is the double array that stores the
 * vector, the second is the length of the array. Cf. lr_xvectr().
 */
typedef double (*Xfuncn_)(const double *, int);

/* ---- PROTOTYPES ---- */

long lr_svd(long Nx, long Nm, const double *Xmeas[], const double Ymeas[], const double Wgt[], 
	    double *Sres, double P[], double Sdev[], Trimat_ Correl,
	    double *Cond);
long lr_xtr(long Np, long Nm, const double Xmeas[], const double Ymeas[],
	    const double Wgt[], const Xfunc1_ Xtrs[], double *Sres, double P[], 
	    double Sdev[], Trimat_ Correl, double *Cond);
long lr_xvectr(long Np, long Nx, long Nm, const double *Xmeas[], const double Ymeas[],
	    const double Wgt[], const Xfuncn_ Xtrs[], double *Sres, double P[], 
	    double Sdev[], Trimat_ Correl, double *Cond);

/* ==== END OF HEADER lrsvd.h ==== */

#endif
