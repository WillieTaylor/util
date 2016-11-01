#ifndef __CGR_H__
#define __CGR_H__

/* ==== HEADER cgr.h ==== */

/* Routine collection for processing centre-of-gravity data. */

/* ANSI C, Iris Indigo IRIX 4.0.5, 30. Oct. 1992. Andris Aszodi */

/* ---- INCLUDE FILES ---- */

#include "pdb.h"

/* ---- PROTOTYPES ---- */

extern Atom_ *pdb_cgrav(Chain_ *Chain);
extern Atom_ cgr_centroid(Atom_ Cgravs[], int Aano);

/* ==== END OF HEADER cgr.h ==== */

#endif
