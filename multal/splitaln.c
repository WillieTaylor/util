#include "util/wt/incl/util.h"

main (argc, argv) int argc; char *argv[]; {
int	nblock = 0, seed = 0;
FILE	*aln, *head, *tail;
	if (argc > 1) aln = fopen(argv[1],"r");
		else  aln = fopen("final.aln","r");
	head = fopen("head.aln","w");
	tail = fopen("tail.aln","w");
	while (1) { char line[555]; int io;
		io = read_line(aln,line);
		if (!io ) continue;
		if (io<0) break;
		if (nblock==1 && strstr(line,"SEED:")) seed = 1;
		if (!strncmp(line,"Block ",6)) nblock++;
		if (nblock==1) fprintf(head,"%s\n", line);
		if (nblock >1) fprintf(tail,"%s\n", line);
	}
	printf("%d\n", seed);
}
