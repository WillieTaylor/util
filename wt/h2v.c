/* cc h2v.c -o h2v util.o -lm
*/
#include "util.def"
/*	converts horrizontal sequences to vertical format	*/
main() {
	int	i, j, ngap, nseqs, maxres;
	char seqs[100][255], text[255];
	FILE *horr, *vert;
	horr = fopen("horr.aln","r");
	vert = fopen("vert.aln","w");
	read_line(horr,text);
	sscanf(text+10,"%d", &nseqs);
	fprintf(vert,"%d sequences\n", nseqs);
	for (i=0; i<nseqs; i++) {
		read_line(horr,text);
		fprintf(vert,"%s\n", text);
	}
	next_line(horr);
	while (!feof(horr)) {
		maxres=0;
		for (i=0; i<nseqs; i++)
		{ int l;
			l = read_line(horr,seqs[i]);
			if (l>maxres) maxres = l;
		}
		next_line(horr);
		for (i=0; i<maxres; i++) {
			ngap = 0;
			for (j=0; j<nseqs; j++) {
				if (!isalpha(seqs[j][i])) {
					ngap++;
					seqs[j][i] = '-';
				}
			}
			if (ngap<nseqs) {
				for (j=0; j<nseqs; j++) fprintf(vert,"%c", seqs[j][i]);
				fprintf(vert,"\n");
			}
		}
	}
}
