/* ==== HEADER histo.h ==== */

#ifndef __HISTO_H__
#define __HISTO_H__

/* A simple averaging, SD and histogram-constructing function collection
 for double data vectors. */

/* ANSI C, IRIX 4.0.5, 1. Feb. 1993. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <math.h>

/* ---- PROTOTYPES ---- */

void avg_sd(double Data[], int Datno, double *Avg, double *Sd);
void make_hist(double Data[], int Datno, int Binno,
	double *Low, double *Step, int Bins[]);

#endif

/* ==== END OF HEADER histo.h ==== */
