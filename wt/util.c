#include <stdio.h>
#include <math.h>

unsigned int pack (i,j,n)	int	i, j; unsigned int n;
/*	references top half (+diag) of a 2D array (N*N) as a vector	*/
{	long long	ip, jp;
	if(n <= 0) { printf("Bad dimension in IPACK\n"); return -1; }
	if(i <= 0) { printf("Bad I in IPACK: i=%d\n",i); return -2; }
	if(j <= 0) { printf("Bad J in IPACK: j=%d\n",j); return -3; }
	ip = max(i,j);
	jp = min(i,j)-1;
	return (unsigned int)(ip-jp+n*jp-(jp*jp-jp)/2);
}

unpack (ip,jp,nx,idx)	int *ip, *jp, nx; unsigned int idx;
/*	references top half (+diag) of a 2D array (N*N) from a vector	*/
{	long long	i, j, n=nx, id=idx;
	double 	b = n+0.5,
		a = b*b-2*id;
	// (92681*92681-92681)/2 = 4294837540 (92681id = 4294930220, 92682id > IMAXu)
	if(n>92681) printf("Error in UNPACK: n=%d is too big\n",nx);
	if(a<=0.0) printf("Error in UNPACK: n=%d, id=%d\n", nx,idx);
	j = (long long)(b-sqrt(a)-0.00001);
	i = j+id-n*j+(j*j-j)/2;
	j++;
	if (j<=0 || j>n ) printf("Bad J in UNPACK: j=%lld (n=%d id=%u)\n",j,nx,idx);
	if (i<=0 || i>n ) printf("Bad I in UNPACK: i=%lld (n=%d id=%u)\n",i,nx,idx);
	*ip = min(i,j);
	*jp = max(i,j);
}

read_line (file,string) FILE *file; char *string;
{	char	c;
	int	i=0;
	*string = 0;
	while(c=getc(file)) {
		/* printf("%d >%c<\n", c,c); */
		if (feof(file)) return -i-1;
		if (c=='\n') return i;
		string[i] = c;
		string[++i] = 0;
	}
}
next_line(file) FILE *file;
{	char	c;
	while(c=getc(file)) {
		if (feof(file)) return 0;
		if (c=='\n') return 1;
	}
}

min(i,j) int i,j; { if(i<j) return i; else return j; }
max(i,j) int i,j; { if(i>j) return i; else return j; }

//float fmin(i,j) float i,j; { if(i<j) return i; else return j; }
//float fmax(i,j) float i,j; { if(i>j) return i; else return j; }
