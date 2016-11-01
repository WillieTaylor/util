#ifndef __DSSPCA_H__
#define __DSSPCA_H__

/* ==== HEADER dsspca.h ==== */

/* Little program to get secondary structure info and C-alpha coordinates
   from DSSP files. H-bonds and SS-bridges are not read. 
   Had to be re-written from the original 8-Feb-93 code 
*/

/* ANSI C, IRIX 5.3, 16. May 1995. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>

/* ---- GLOBAL TYPES ---- */

typedef struct        /* main chain entry */
{
    int Resno;   /* residue no. */
    char Chain;  /* chain identifier , <SPACE> if 1 chain only */
    char Res;    /* AA code, 1-letter */
    char Secstruct; /* secondary structure code, see ref. */
    char Turns3, Turns4, Turns5;    /* helix turn labels */
    char Bend;	    /* geometrical bend */
    char Chir;	    /* chirality: '+' or '-' */
    char Bridge1, Bridge2;  /* beta labels: CAPS if antiparallel, lowercase if parallel */
    int Beta1, Beta2;	    /* partner resnums in beta-bridges */
    char Sheet;		    /* beta-sheet label (always caps) */
    int Access;  /* solvent accessibility, see ref. */
    float Xca,Yca,Zca; /* C-alpha coordinates */
} Dssprec_;

/* ---- PROTOTYPES ---- */

/* dssp_ca: reads the (opened) Dssp stream. This is a text file produced
   by DSSP, Version Oct. 1988 (ref: W. Kabsch, C. Sander, Biopolymers
   22:2577-2637 (1983)). Returns Entry, an array of size Nres, containing
   chain ID, residue no, 1-letter AA code, Kabsch/Sander secondary
   structure code and solvent accessibility for each residue. 
   Returns the size of Entry (i.e. no. of AAs) or 0 on error.
*/
int dssp_ca(FILE *Dssp, Dssprec_ **Entry, int *Chainno);

/* dssp_cadist(): returns the CA:CA distance between two DSSP records
 * pointed to by Dp1 and Dp2. Checks if they're NULL.
 */
double dssp_cadist(const Dssprec_ *Dp1, const Dssprec_ *Dp2);

/* ==== END OF HEADER dsspca.h ==== */

#endif
