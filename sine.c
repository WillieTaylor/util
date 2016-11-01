#include "wt/incl/util.h"
main() {
int i; float x;
int f[100];
	for (i=0; i<100; i++) f[i] = 0;
	for (i=0; i<99999999; i++) {
		x = sin((float)i);
		if (drand48()>0.1) continue;
		f[abs((int)(100.0*x))]++;
	}
	for (i=0; i<100; i++) printf("%5d %4d\n",i,f[i]);
}

