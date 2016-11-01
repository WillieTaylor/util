#ifndef __DATETIME_H__
#define __DATETIME_H__

/* ==== HEADER datetime.h ==== */

/* Routines for producing the time and date in ASCII format. */

/* ANSI C, Iris Indigo IRIX 4.0.5, 6. Nov. 1992. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <stdio.h>
#include <time.h>
#include <strings.h>
#include <math.h> 	/* for floor */

/* ---- PROTOTYPES ---- */

extern time_t time_stamp(char *Asctime, char *Nicetime, char *Shortime);
extern long julian_date(int Day, int Month, int Year);

/* ==== END OF HEADER datetime.h ==== */

#endif
