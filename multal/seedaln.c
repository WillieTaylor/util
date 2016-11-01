#include "util/wt/incl/util.h"

main (argc, argv) int argc; char *argv[]; {
int	nblock = 0, seed = 0;
FILE	*aln, *out;
	if (argc > 1) aln = fopen(argv[1],"r");
		else  aln = fopen("final.aln","r");
	out = fopen("seed.aln","w");
	while (1) { char line[555]; int io;
		io = read_line(aln,line);
		if (!io ) continue;
		if (io<0) break;
		if (strstr(line,"SEED:")) seed = 1;
		if (!strncmp(line,"Block ",6)) {
			if (nblock) {
				if (seed) {
					printf("1\n");
					exit(1);
				} else {
					fclose(out);
					out = fopen("seed.aln","w");
				}
			}
			nblock++;
		}
		fprintf(out,"%s\n", line);
	}
        if (seed) printf("1\n"); else printf("0\n");
}
