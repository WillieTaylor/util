/*
*/
#define PI 3.14159265358979323846
#define twoPI PI+PI

/*	STRUCTURES FOR 3D GEOMETRY	*/

typedef struct	{ float	x, y, z; } Vec;
typedef struct  { Vec A, B, C; } Mat;

float   Mdet();
float	vmod();
float	vsqr();
float	vdif();
float	vddif();
float	vdot();
float	vtri();
float	pvol();
float	pdotp();
float	phand();
float	dotOline();
float	endOline();
float	dist2line();
float	lineOline();
float	angle();
float	torsion();
float	angle1pi();
float	angle2pi();
float	angdif();

int	dot2pair();
