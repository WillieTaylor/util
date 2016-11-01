#ifndef __VECT_H__
#define __VECT_H__

/* ==== HEADER vect.h ==== */

/* Useful collection of 3D vector routines. */

/* ANSI C, Iris Indigo IRIX 5.2, 26. July 1994. Andris */

/* ---- INCLUDE FILE ---- */

#include <math.h>
#include <stdio.h>

/* ---- GLOBAL CONSTANT ---- */

#define VERYSMALL (1.0e-15)

/* ---- GLOBAL TYPE ---- */

typedef struct     /* 3D vector type */
{
  double X,Y,Z;
} Vector_;

/* ---- PROTOTYPES ---- */

extern double vector_length(const Vector_ *V);
extern double unit_vector(const Vector_ *V, Vector_ *U);
extern void sum_vector(const Vector_ *V, const Vector_ *W, Vector_ *S);
extern void diff_vector(const Vector_ *V, const Vector_ *W, Vector_ *D);
extern double diff_length(const Vector_ *V, const Vector_ *W);
extern double scalar_product(const Vector_ *V, const Vector_ *W);
extern void vector_product(const Vector_ *V, const Vector_ *W, Vector_ *P);
extern double dir_cosine(const Vector_ *V, const Vector_ *W);
extern double norm_diff(const Vector_ *V, const Vector_ *W, Vector_ *N);
extern void list_vector(const Vector_ *V, const char *Format);

/* ==== END OF HEADER vect.h ==== */

#endif
