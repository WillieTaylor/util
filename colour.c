#include "wt/incl/util.h"
/* cc colour.c -o colour wt/util.o -lm
*/
#include <stdio.h>
#include <alloca.h>
#define NSEQS 500
#define NACID 30
main(argc,argv) int argc; char *argv[];{
FILE	*aln;
char	**seqs;
int	nseqs, length;
char	pct = '%';
float	s;
	aln = fopen("final.aln","r");
	getaln(aln,&seqs,&nseqs,&length);
	printf("%c!PS-Adobe-2.0\n", pct);
	printf("%c%cBoundingBox: 0 0 %d %d\n", pct,pct,length*10,nseqs*15);
	printf("%c Alignment length = %d\n", pct, length);
	printf("%c Number of sequences %d\n", pct, nseqs);
	printf("%c%cEndComments\n", pct,pct);
	s = 58.0/(float)length;
	if (argc>2 && *argv[2]=='p') printf("%5.3f %5.3f scale\n", s,s);
	colour(*argv[1],seqs,nseqs,length);
}

colour (mode,seqs,nseqs,length) 
char	mode, **seqs;
int	nseqs, length;
{
/*
aspartate  (D) ASP = BRR redish-redey-red                (red) = 1.0 0.0 0.0
serine     (S) SER = YRR yellowish-redey-red         (scarlet) = 1.0 0.2 0.0
threonine  (T) THR = RYR redish-yellowey-red      (vermillion) = 1.0 0.4 0.0
glycine    (G) GLY = YRY yellowish-redey-yellow       (orange) = 1.0 0.6 0.0
proline    (P) PRO = RYY redish-yellowey-yellow    (tangerine) = 1.0 0.8 0.0
cystine    (C) CYS = YYY yellowish-yellowey-yellow    (yellow) = 1.0 1.0 0.0
alanine    (A) ALA = GYY greenish-yellowey-yellow      (lemon) = 0.8 1.0 0.0
valine     (V) VAL = YGY yellowish-greeney-yellow (lemon-lime) = 0.6 1.0 0.0
isoleucine (I) ILE = GYG greenish-yellowey-green        (lime) = 0.4 1.0 0.0
leucine    (L) LEU = YGG yellowish-greeney-green       (grass) = 0.2 1.0 0.0
methionine (M) MET = GGG greenish-greeney-green        (green) = 0.0 1.0 0.0
phenylala. (F) PHE = BGG blueish-greeney-green       (emerald) = 0.0 1.0 0.4
tyrosine   (Y) TYR = GBG greenish-blueey-green     (turquoise) = 0.0 1.0 0.8
tryptophan (W) TRP = BGB blueish-greeney-blue          *(cyan) = 0.0 0.8 1.0
histidine  (H) HIS = GBB greenish-blueey-blue        (peacock) = 0.0 0.4 1.0
arginine   (R) ARG = BBB blueish-blueey-blue            (blue) = 0.0 0.0 1.0
lysine     (K) LYS = RBB redish-blueey-blue           (indigo) = 0.4 0.0 1.0
asparagine (N) ASN = BRB blueish-redey-blue           (purple) = 0.8 0.0 1.0
glutamine  (Q) GLN = RBR redish-bluey-red           *(magenta) = 1.0 0.0 0.8
glutamate  (E) GLU = BRR blueish-redey-red            (violet) = 1.0 0.0 0.4
*/
float	rgb[26][3] = {
0.8, 1.0, 0.0,
0.0, 0.0, 0.0,
1.0, 1.0, 0.0,
1.0, 0.0, 0.0,
1.0, 0.0, 0.4,
0.0, 1.0, 0.4,
1.0, 0.6, 0.0,
0.0, 0.4, 1.0,
0.4, 1.0, 0.0,
0.0, 0.0, 0.0,
0.4, 0.0, 1.0,
0.2, 1.0, 0.0,
0.0, 1.0, 0.0,
0.8, 0.0, 1.0,
0.0, 0.0, 0.0,
1.0, 0.8, 0.0,
1.0, 0.0, 0.8,
0.0, 0.0, 1.0,
1.0, 0.2, 0.0,
1.0, 0.4, 0.0,
0.0, 0.0, 0.0,
0.6, 1.0, 0.0,
0.0, 0.8, 1.0,
1.0, 1.0, 1.0,
0.0, 1.0, 0.8,
0.0, 0.0, 0.0 
};
int	i, j, x0=0, y0=0, cw=10, ch=15;
int	naa, aa[NACID];
float	na, ns = (float)nseqs;
char	pct = '%';
	printf("/Times-Roman findfont\n12 scalefont\nsetfont\n");
	printf("newpath\n");
	for (i=0; i<length; i++)
	{ int	x = x0 + i*cw;
	  float r, g, b;
		for (j=0; j<NACID; j++) aa[j] = 0;
		printf("\n%c ", pct);
		r = g = b = 0.0;
		for (j=0; j<nseqs; j++)
		{ char	s = toupper(seqs[j][i]);
		  int	y = y0 + j*ch,
			a = s-'A';
			if (mode=='4' && seqs[j][i]==' ') { r += 1.0; g += 1.0; b += 1.0; }
			if (s<'A' || s>'Z') continue;
			aa[a] = 1;
			printf("%c", s);
			r += rgb[a][0]; g += rgb[a][1]; b += rgb[a][2];
		}
		printf("\n");
		r /= ns; g /= ns; b /= ns;
		naa = 0;
		for (j=0; j<NACID; j++) naa += aa[j];
		if (!naa) naa = 20;
		na = (float)naa;
		na = (na-1.0)/30.0;
		if (mode=='3') { r -= na; g -= na; b -= na; } /* +black */
		if (mode=='4') { r += na; g += na; b += na; } /* +white */
		if (r>1.0) r=1.0; if (g>1.0) g=1.0; if (b>1.0) b=1.0;
		if (r<0.0) r=0.0; if (g<0.0) g=0.0; if (b<0.0) b=0.0;
		for (j=0; j<nseqs; j++)
		{ char	s = toupper(seqs[j][i]);
		  int	y = y0 + j*ch,
			a = s-'A';
			if (mode=='1') {
				/* individual colours */
				if (s<'A' || s>'Z') r = g = b = 0.0;
				else { r = rgb[a][0]; g = rgb[a][1]; b = rgb[a][2]; }
			}
			printf("%d %d moveto\n", x,y);
			printf("0 15 rlineto\n10 0 rlineto\n0 -15 rlineto\nclosepath\n");
			printf("%6.3f %6.3f %6.3f setrgbcolor\nfill\n0 0 0 setrgbcolor\n", r,g,b);
			printf("%d %d moveto\n(%c) show\n", x+1, y+3, seqs[j][i]);
		}
	}
	printf("showpage\n");
}


getaln (aln, seq, nseq, length)
FILE    *aln;
char	***seq;
int	*nseq, *length;
{
char	**seqs;
int     i,j, k, n, len, nseqs, inseq;
int     maxlen, minlen, gapsin, skip[NSEQS];
char    line[2222];
char    **tmp;
int     N, C, seqlen;
        while (read_line(aln,line)>0) if (strstr(line,"Block")) break;
        read_line(aln,line);
        sscanf(line,"%d", &inseq);
        seqs = (char**)malloc(sizeof(char*)*inseq);
        tmp = (char**)alloca(sizeof(char*)*inseq);
        nseqs = 0;
        for (i=0; i<inseq; i++) {
                skip[i] = 1;
                read_line(aln,line);
                if (!strncmp(line,"SKIP",4)) continue;
                if (!strncmp(line,"skip",4)) continue;
                skip[i] = 0;
                tmp[nseqs] = (char*)alloca(sizeof(char)*2222);
                nseqs++;
        }
        gapsin = nseqs/4;
        if (nseqs<4) gapsin = 1;
gapsin = 999;
        n = 0; N = -1;
        while (read_line(aln,line)>0)
        { int   gaps = 0, ns = 0;
                for (i=0; i<inseq; i++) {
                        if (skip[i]) continue;
                        tmp[ns++][n] = line[i];
                        if (line[i] == '-') gaps++;
                }
                if (gaps<gapsin && N<0) N = n;
                n++;
                if (gaps<gapsin) C = n;
        }
        seqlen = len = C - N;
        maxlen = 0;
        minlen = 99999;
        for (i=0; i<nseqs; i++) {
                seqs[i] = (char*)malloc(sizeof(char)*len);
                n = 0;
                for (j=N; j<C; j++)
                { char  r = tmp[i][j], R = toupper(r);
                        if (r=='-') {
                                seqs[i][n] = ' ';
                        } else {
                                if (R<'A' || R>'Z') { seqs[i][n] = 'X'; }
                                        else        { seqs[i][n] = R; }
                        }
                        n++;
                }
                if (seqlen<minlen) minlen = seqlen;
                if (seqlen>maxlen) maxlen = seqlen;
        }
	*seq = seqs;
	*nseq = nseqs;
	*length = maxlen;
        return nseqs;
}

matin (file,mat)
        char    *file;
        int     **mat;
{
        int     i, j, mat_const;
        char    acid[NACID], c;
        FILE    *mat_file;

        mat_file = fopen(file,"r");
        while( c = getc(mat_file), c != '\n' ) putchar(c); printf("\n");
        fscanf(mat_file,"%s\n",acid);
        printf("%s\n",acid);
        fscanf(mat_file,"%d\n",&mat_const);
        printf("matrix constant = %d\n",mat_const);
        for (i=0; i<NACID; i++) for (j=0; j<NACID; j++) mat[i][j] = 0;
        for( i = 0; acid[i]; i++ )
        {       int     ai = acid[i]-'A';
                for( j = 0; acid[j]; j++ )
                {       int aj = acid[j]-'A';
                        fscanf(mat_file,"%d",&mat[ai][aj]);
                        mat[ai][aj] += mat_const;
                }
        }
        mat[1][1] = mat[9][9] = mat[14][14] = mat[20][20] = mat[23][23] = mat[25][25] = -1;
}


