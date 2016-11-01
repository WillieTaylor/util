#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <alloca.h>
#include <time.h>

#define LOWER(c)	((c)<'a'?((c)-'A'+'a'):(c))
#define UPPER(c)	((c)<'a'?(c):((c)-'a'+'A'))
#define DIST3v(v)	sqrt( (v).x*(v).x + (v).y*(v).y + (v).z*(v).z )
#define DDIST3v(v)	(v).x*(v).x + (v).y*(v).y + (v).z*(v).z
#define DIST3(x,y,z)	sqrt((x)*(x)+(y)*(y)+(z)*(z))
#define PRINTi(i)     printf("%d",i);
#define Pi(i) printf(" " #i " = %d",i);
#define Pu(u) printf(" " #u " = %u",u);
#define PRINTr(r)     printf("%f",r);
#define Pr(r) printf(" " #r " = %f",r);
#define PRINTc(x)     printf("%c",x);
#define Pc(x) printf(" " #x " = %c",x);
#define PRINTs(x)     printf("%s",x);
#define Ps(x) printf(" " #x " = %s",x);
#define	PRINTv(v)	printf("%f %f %f", v.x, v.y, v.z); 
#define	Pv(v) printf(" " #v " = %f %f %f", v.x, v.y, v.z); 
#define PRINTt(x) printf("" #x "");
#define Pt(x) printf(" " #x " ");
#define Px printf("here\n");
#define DO(i,n)		for (i=0; i<n; i++)
#define DO1(i,n)	for (i=1; i<=n; i++)
#define SP		printf(" ");
#define SPP		printf("  ");
#define SPPP		printf("   ");
#define TP		printf("\t");
#define TPP		printf("\t\t");
#define TPPP		printf("\t\t\t");
#define NL		printf("\n");
#define NLL		printf("\n\n");
#define NLLL		printf("\n\n\n");
#define TEST(obj) 	if (obj == NULL) { printf("malloc fail for " #obj "\n"); exit(1); }
#define REST(obj) 	if (obj == NULL) { printf("realloc fail for " #obj "\n"); exit(1); }

#define YES 1
#define NO  0
#define SET   1
#define UNSET 0
//# define LIVE 1
#define DEAD 0
#define TRUE  1
#define FALSE 0
#define LEFT -1
#define RIGHT 1
#define BIG 1000000000
#define WEE -BIG
#define NOISE 0.00001
#define BYTE  255

#define CHAR  127
#define CHARu 255
#define SMAX  32767
#define SMAXu 65535
#define IMAX  2147483647
#define IMAXu 4294967295
#define LMAX  9223372036854775807
#define LMAXu 18446744073709551615

unsigned int pack();

//float fmin();
//float fmax();
