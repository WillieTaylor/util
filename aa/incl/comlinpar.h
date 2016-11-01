#ifndef __COMLINPAR_H__
#define __COMLINPAR_H__

/* ==== HEADER comlinpar.h ==== */

/* For the uniform processing of command-line parameters. */

/* ANSI C, IRIX 5.2, 26. July 1994. Andris */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* ---- DEFINITIONS ---- */

#define CMDSW_PREFIX '-'  /* UNIX convention; set to '/' for DOS */

/* ---- PROTOTYPES ---- */

int getcmd_switch(int argc, char *argv[], const char *Switch);
int getcmd_fname(int argc, char *argv[], int No);
char **getcmd_fnames(int argc, char *argv[], int *Fnos);

int getcmd_strpar(int argc, char *argv[], const char *Parname, char **Str);
int getcmd_intpar(int argc, char *argv[], const char *Parname, int *Val);
int getcmd_dblpar(int argc, char *argv[], const char *Parname, double *Val);


/* ==== END OF HEADER comlinpar.h ==== */

#endif
