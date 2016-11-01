#ifndef __KOLMS_H__
#define __KOLMS_H__

/* ==== HEADER kolms.h ==== */

/* Kolmogorov/Smirnov test for two sampled distributions.
 * Based on Numerical Recipes.
 */

/* ANSI C, IRIX 5.2, 17. Feb. 1995. Andris */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* ---- PROTOTYPES ---- */

/* ks_2(): the Kolmogorov/Smirnov test for two sampled distributions.
 * Cdf1, 2 are discrete cumulative distribution functions supplied as
 * tables containing N1, 2 sample points where the steps occur.
 * Each step is 1/N1 (or 1/N2) high.
 * Cdf1, 2 are assumed to have been sorted into ascending order
 * before the call: not tested.
 * Returns the K/S statistics D (the maximal difference between Cdf1, 2)
 * in *Kstat if Kstat!=NULL.
 * Return value: the probability that the two distributions are the same
 * or -1.0 on error.
 */
double ks_2(const double Cdf1[], int N1, const double Cdf2[], int N2, 
	    double *Kstat);

/* ==== END OF HEADER kolms.h ==== */

#endif
