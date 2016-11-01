#ifndef __PDBDRAW_H__
#define __PDBDRAW_H__

/* ==== HEADER pdbdraw.h ==== */

/* Draws a colour-coded PDB chain which can
 be interactively rotated/moved on the screen by the mouse. */

/* PORTABILITY WARNING: Highly Silicon Graphics-specific code! */
/* Link in the SGI shared graphics library (-lgl_s) into executable */

/* ANSI C (otherwise...), Silicon Graphics Iris Indigo, IRIX 4.0.5 */

/* 31. Jan. 1994, Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>

/* ---- GRAPHICS (NONSTANDARD) HEADERS ---- */

#include <gl.h>		/* SGI graphics library */
#include <device.h>	/* SGI event handler library */

/* ---- INCLUDE FILES ---- */

#include "pdb.h"	/* PDB reader and types */

/* ---- DEFINITIONS ---- */

/* 32 bit colour constants, full brightness. The alpha byte is set to
 FF in all cases (no transparency). The byte order is ABGR. */
#define RGB_WHITE 0xFFFFFFFF
#define RGB_YELLOW 0xFF00FFFF
#define RGB_MAGENTA 0xFFFF00FF
#define RGB_CYAN 0xFFFFFF00
#define RGB_RED 0xFF0000FF
#define RGB_GREEN 0xFF00FF00
#define RGB_BLUE 0xFFFF0000
#define RGB_BLACK 0xFF000000

/* default window position */
#define CXORIG 100
#define CYORIG 10

/* symbolic constants for manipulating the static variables in
 * "pdbdraw.c" (menu options etc.) See set_pdrvars()
 */
#define PDR_NOP 0	    /* do nothing */
#define PDR_DEPTH 1	    /* depth toggle switch */
#define PDR_CALPHA 2	    /* C-alpha or full main chain switch */
#define PDR_SIDECH 3	    /* side chain draw toggle */
#define PDR_HBOND 4	    /* main-chain H-bond draw toggle */
#define PDR_GOURAUD 5	    /* Gouraud shade toggle */
#define PDR_VIEW 7	    /* viewing angle set */
#define PDR_LCD 8	    /* largest coordinate set */
#define PDR_INV 9	    /* invert background colour */
#define PDR_COL 10	    /* colour scheme adjustment */
#define PDR_LABEL 11	    /* labels on C-alpha atoms */

#define PDR_ON 1	    /* switches */
#define PDR_OFF 0

#define VIEW_ORTHO 1	    /* projection angle constants */
#define VIEW_WIDE 2
#define VIEW_NORM 3
#define VIEW_TELE 4
#define VIEW_DFLT 5

/* colour scheme constants */
#define COLOUR_LEVITT 1	    /* Levitt's hydrophobicity */
#define COLOUR_BINPHOB 2    /* binary hydrophobicity */
#define COLOUR_OCCU 3	    /* occupancy */
#define COLOUR_BFACT 4	    /* B-factor */
#define COLOUR_OCBF 5	    /* occupancy*Bfactor */
#define COLOUR_USER 6	    /* user-defined colour table */

/* ---- TYPES ---- */

typedef struct		/* one atom to be drawn */
{
    float X[3];	/* coords */
    unsigned long Col;	/* colour */
    Str3_ Id;	/* atom ID */
    float Occu, Bfact;	/* occupancy and B-factor */
    int Aaord;	/* diff of AA code from 'A' */
}
Coord_ ;

typedef struct	    /* a pair of H-bonded atoms */
{
    int i, j;	    /* atom indices in Coords */
}
Hpair_ ;

typedef struct		/* the whole chain to be drawn */
{
    Coord_ *Coords;	/* atom coordinates,colour and IDs */
    int Cono;
    Hpair_ *Hpairs;	/* H-bond entries: pairwise */
    int Hpno;
    float Rot[3][3];	/* rotation matrix */
}
Drawchain_ ;

typedef unsigned long Colcode_[26];	/* code table for AA->colour conversion */

/* ---- GLOBAL VARS ---- */

extern Colcode_ Phob2col;
extern Colcode_ Levittcol;

/* ---- PROTOTYPES ---- */

#ifdef __cplusplus
extern "C" {
#endif

Drawchain_ *transfer_pdbcoords(Chain_ *Chain, Colcode_ Col);
void destruct_drawchain(Drawchain_ *Drawchain);
long init_draw(char *Windowtitle, short Perspang, short Linewidth, 
	long Xorig, long Yorig);
void display_pdbchain(long Gid, Drawchain_ *Drawchain);
void draw_pdbchain(long Gid, Drawchain_ *Drawchain);
void set_pdrvars(int Varsym, ...);

#ifdef __cplusplus
}
#endif

/* ==== END OF HEADER pdbdraw.h ==== */

#endif
