#ifndef __STAT2_H__
#define __STAT2_H__

/* ==== HEADER stat2.h ==== */

/* One- and two-variable statistics: mean, SD and correlation
 * for double data.
 */

/* ANSI C, IRIX 5.2, 17. Jan. 1995. Andris */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ---- TYPES ---- */

/* Statrec_ : holds the sum and squared sum of X.
 * N is the number of data points.
 */
typedef struct
{
    double Sx, Sx2;
    unsigned long N;
}
Statrec_ ;

/* Stat2rec_ : holds the sums and squared sums and the mixed
 * sum of X and Y. N is the number of data points.
 */
typedef struct
{
    double Sx, Sy, Sx2, Sy2, Sxy;
    unsigned long N;
}
Stat2rec_ ;

/* ---- DEFINITIONS ---- */

/* The following symbolic constants control the calculation
 * of the avg and SD for the 1-var statistics.
 */
#define STAT_NOOP 0
#define STAT_AVG 1
#define STAT_SD 2

/* The following symbolic constants control the calculation
 * of the avg's and SD-s for the 2-var statistics.
 */
#define STAT2_NOOP 0
#define STAT2_AVGX 1
#define STAT2_SDX 2
#define STAT2_AVGY 3
#define STAT2_SDY 4
#define STAT2_CORR 5

/* ---- PROTOTYPES ---- */

#ifdef __cplusplus
extern "C" {
#endif

void init_stat(Statrec_ *St);
void add_stat(double X, Statrec_ *St);
double eval_stat(const Statrec_ *St, int Op);

void init_stat2(Stat2rec_ *St);
void add_stat2(double X, double Y, Stat2rec_ *St);
double eval_stat2(const Stat2rec_ *St, int Op);

#ifdef __cplusplus
}
#endif

/* ==== END OF HEADER stat2.h ==== */
#endif
