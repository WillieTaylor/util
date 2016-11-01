#ifndef __PDF_H__
#define __PDF_H__

/* ==== HEADER pdf.h ==== */

/* Routine collection for estimating probability density functions. */

/* ANSI C, IRIX 4.0.5, 28. Sept. 1993. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* ---- GLOBAL TYPES ---- */

/* Pdfrec_ : this structure holds the variables for estimating a PDF.
 * X is the indep. var, Y is the PDF, Yint holds the integer counts
 * of the values falling in each bin. The bins are centred around the
 * values in X[]. X, Y, Yint are Bno long arrays. Wid is the half-bin width, 
 * N is the total number of values processed.
 */
typedef struct
{
    double *X, *Y;
    int *Yint;
    double Wid;
    int N, Bno;
}
Pdfrec_ ;

/* ---- PROTOTYPES ---- */

int init_pdf(int Binno, double Binwidth, double Low, double Up, 
		Pdfrec_ *Pdf);
int add_pdf(double Value, Pdfrec_ *Pdf);
int eval_pdf(Pdfrec_ *Pdf);
void remove_pdf(Pdfrec_ *Pdf);

/* ==== END OF HEADER pdf.h ==== */
#endif
