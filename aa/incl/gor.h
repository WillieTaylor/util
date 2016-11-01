#ifndef __GOR_H__
#define __GOR_H__

/* ==== HEADER gor.h ==== */

/* The GOR prediction method. */

/* ANSI C, IRIX 4.0.5, 11. Oct. 1993. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>

/* ---- GOR TYPES ---- */

/* State_ : holds the GOR parameter table (Table) and the structural
 * state name. The table has as many rows as there are amino acids
 * in the permutation string (see below) and has GOR_POS columns.
 * Column[0] corresponds to the j-GOR_HALFPOS-th entry.
 * Dc is the decision constant which must be subtracted from the
 * predicted information values to avoid periodic structure
 * overprediction.
 */
typedef struct
{
    char *Name;
    int **Table;
    int Dc;
}
State_ ;

/* Gortable_ : stores the amino acid permutation string Aaperm,
 * the number of structural states (Statno), the length of Aaperm
 * (Aano) and a Statno-long array (States[]) whose elements are
 * the State_ records (see above) for each structural state.
 * Descr is a string describing the table (the 1st line in the
 * table file).
 */
typedef struct
{
    char *Descr, *Aaperm;
    int Statno, Aano;
    State_ *States;
}
Gortable_ ;

/* --- PROTOTYPES ---- */

Gortable_ read_gortable(const char *Fname);
char *gor_predict(const char *Seq, Gortable_ Gtbl, int *Info[]);
void gor_pred(const char *Seq, Gortable_ Gtbl, int *Info[]);
char pred_char(Gortable_ Gtbl, const int Infos[]);
void clean_gortable(Gortable_ *Gtbl);

/* ==== END OF HEADER gor.h ==== */
#endif
