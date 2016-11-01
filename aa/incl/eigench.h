#ifndef __EIGENCH_H__
#define __EIGENCH_H__

/* ==== HEADER eigench.h ==== */

/* Calculates a given number of eigenvalues and eigenvectors by
 matrix Chebyshev polynomials. */

/* ANSI C, Iris Indigo IRIX 4.0.5, 23. Nov. 1992. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* ---- INCLUDE ---- */

#include "matrix.h"

/* ---- PROTOTYPES ---- */

int eigen_ch(Trimat_ Matr, int Rno, double Eval[], Sqmat_ Evec,
	int Eno, double Toler, int Iterno);

/* ==== END OF HEADER ==== */

#endif
