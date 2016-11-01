#ifndef __PLOTDISTMAT_H__
#define __PLOTDISTMAT_H__

/* ==== HEADER plotdistmat.h ==== */

/* Plots a triangular distance matrix using the Silicon Graphics
 * GL language. This program is therefore NOT PORTABLE. */

/* ANSI C and SGI GL, SGI IRIX 4.0.5   17. Nov. 1993, Andris Aszodi */

/* ---- INCLUDE FILES ---- */

#include "matrix.h"
#include <math.h>
#include <gl/gl.h>      /* Silicon Graphics GL graphics library headers: NOT PORTABLE!!! */

/* ---- DEFINES ---- */

#define DXORIG 100
#define DYORIG 460

/* ---- TYPES ---- */

typedef struct		/* color scale record */
{
    float Threshold;	/* an entry below Threshold */
    long Color;		/* will have this colour */
}
Colscale_ ;

/* ---- GLOBAL VARS ---- */

/* the default colour scale */
#define SCALENO 5
extern Colscale_ Colscale[SCALENO];

/* ---- PROTOTYPES ---- */

long init_plot(int Mno, char *Title, long Xorig, long Yorig);
void plot_distmat(long Gid, Trimat_ Mat, int Mno,
		    Colscale_ Colscale[], int Scaleno);
void plot_lumat(long Gid, Trimat_ Lm, Trimat_ Um, int Mno,
		    Colscale_ Lsc[], Colscale_ Usc[], 
		    int Lscno, int Uscno);

/* ==== END OF HEADER plotdistmat.h ==== */

#endif
