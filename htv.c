#include "wt/incl/util.h"
main() {
int i,j, in, n = 0;
#define N 2000
char line[N][N];
FILE	*txt;
	txt = fopen("text.dat", "r");
	for (i=0; i<N; i++) for (j=0; j<N; j++) line[i][j] = (char)0;
	while (1) { int io = read_line(txt,line[n]);
		if (io<1) break;
		n++;
	}
	for (j=0; j<N; j++) {
		in = 0;
		for (i=0; i<n; i++) { char c = line[i][j];
			if (!c) { printf(" "); }
			 else   { printf("%c", line[i][j]); in = 1; }
		}
		printf("\n");
		if (!in) exit(1);
	}
}

