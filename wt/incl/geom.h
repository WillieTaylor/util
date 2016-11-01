/*
*/
#define PI 3.14159265358979323846
#define twoPI PI+PI

/*	STRUCTURES FOR 3D GEOMETRY	*/

typedef struct	{ float	x, y, z; } Vec;
typedef struct  { Vec A, B, C; } Mat;

void  vinit ( Vec* );
void  vcopy ( Vec, Vec* );
void  vunit ( Vec, Vec* );
void  vnorm ( Vec* );
void  vave  ( Vec, Vec, Vec* );
void  vsum  ( Vec, Vec* );
void  vadd  ( Vec, Vec, Vec* );
void  vsub  ( Vec, Vec, Vec*);
void  vmul  ( Vec*, float );
void  vdiv  ( Vec*, float );
void  vmul2 ( Vec, float, Vec* );
void  vdiv2 ( Vec, float, Vec* );
void  vrset ( Vec*, float );
void  vradd ( Vec*, float );
float vdif  ( Vec, Vec );
float vddif ( Vec, Vec );
float vdad  ( Vec, Vec );
float vddad ( Vec, Vec );
void  vat0  ( Vec, Vec, Vec*, float );
void  vatA  ( Vec, Vec, Vec*, float );
float vsqr  ( Vec );
float vmod  ( Vec );
float vdot  ( Vec, Vec );
void  vprod ( Vec, Vec, Vec* );
float vtri  ( Vec, Vec, Vec );
float pdotp ( Vec, Vec, Vec, Vec );
float pvol  ( Vec, Vec, Vec, Vec );
float phand ( Vec, Vec, Vec, Vec );
void identM ( Mat* );
void  VtoM  ( Vec, Vec, Vec, Mat* );
void  MtoM  ( Mat*, Mat* );
void  MtoT  ( Mat*, Mat* );
void Mprint ( Mat* );
void  MmulM ( Mat*, Mat*, Mat* );
void  Minv  ( Mat*, Mat*, float );
float Mdet  ( Mat* );
void  MmulV ( Mat*, Vec, Vec* );
void  VmulM ( Mat*, Vec, Vec* );
void  Mrot  ( char axis, Mat*, float );
void Mframe ( Vec, Vec, Vec, Mat* );
int   line2tri ( Vec, Vec, Vec, Vec, Vec );
int   line2line( Vec, Vec, Vec, Vec, float );
float dist2line( Vec, Vec, Vec, Vec );
int   norm2line( Vec, Vec, Vec, Vec, Vec*, Vec* );
float lineOline( Vec, Vec, Vec, Vec, Vec* );
void  dot2line ( Vec, Vec, Vec, Vec* );
int   dot1pair ( Vec, Vec, Vec, Vec* );
int   dot2pair ( Vec, Vec, Vec, Vec*, Vec* );
int   line2dot ( Vec, Vec, Vec, float );
float endOline ( Vec, Vec, Vec );
float dotOline ( Vec, Vec, Vec, Vec* );
void  rotate   ( Vec, Vec, Vec*, float );
float angle    ( Vec, Vec, Vec );
float torsion  ( Vec, Vec, Vec, Vec );
float angle1pi ( float, float );
float angle2pi ( float, float );
float angdif   ( float, float );
void  separate ( Vec*, Vec*, float, float );
int   fsort4min ( float, float, float, float );
int   fsort4max ( float, float, float, float );
int   isort4min ( int, int, int, int );
int   isort4max ( int, int, int, int );
