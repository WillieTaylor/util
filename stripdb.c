#include "wt/incl/util.h"

main (argc, argv) int argc; char *argv[]; {
int	chain = 0;
FILE	*aln, *head, *tail;
	if (argc > 1) aln = fopen(argv[1],"r");
		else  aln = fopen("temp.pdb","r");
	head = fopen("head.pdb","w");
	tail = fopen("tail.pdb","w");
	while (1) { char line[555]; int io;
		io = read_line(aln,line);
		if (!io ) continue;
		if (io<0) break;
		if (chain==0) fprintf(head,"%s\n", line);
		if (chain >0) fprintf(tail,"%s\n", line);
		if (!strncmp(line,"ENDMOL",6)) chain++; 
		if (!strncmp(line,"TER",3)) chain++; 
	}
}
