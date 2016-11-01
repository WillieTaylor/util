#ifndef __PDF2_H__
#define __PDF2_H__

/* ==== HEADER pdf2.h ==== */

/* Routine collection for estimating probability density functions. */

/* ANSI C, IRIX 4.0.5, 9. Nov. 1993. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* ---- GLOBAL TYPES ---- */

/* Pdfrec2_ : this structure holds the variables for estimating a 
 * joint PDF of 2 random variables.
 * X1, 2 are the indep. vars, Y is the PDF, Yint holds the integer counts
 * of the values falling in each bin. The bins are centred around the
 * values in X1[] (Bno1 long) and X2[] (Bno2 long). Y, Yint are Bno1 x Bno2 arrays.
 * Wid1, 2 are the half-bin widths for the 1st and 2nd variables, 
 * N is the total number of values processed.
 */
typedef struct
{
    double *X1, *X2;
    double **Y;
    unsigned long **Yint;
    double Wid1, Wid2;
    unsigned long N, Bno1, Bno2;
}
Pdfrec2_ ;

/* ---- PROTOTYPES ---- */

int init_pdf2(int Binno1, int Binno2, double Binwidth1, double Binwidth2, 
	    double Low1, double Up1, double Low2, double Up2, 
	    Pdfrec2_ *Pdf);
int add_pdf2(double Value1, double Value2, Pdfrec2_ *Pdf);
int eval_pdf2(Pdfrec2_ *Pdf);
void remove_pdf2(Pdfrec2_ *Pdf);

/* ==== END OF HEADER pdf2.h ==== */

#endif
