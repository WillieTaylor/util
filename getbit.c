/*
cc getcas.c -o getcas wt/util.o aa/pdbprot.o -lm
*/
#include <stdlib.h>
#include <alloca.h>
#include "wt/incl/util.h"
#include "wt/incl/geom.h"
#include "aa/incl/pdbprot.h"
#include "aa/incl/matplot.h"
#include "aa/incl/matrix.h"
#define NALLOC 1000
#define NACID 30

float cadist();
float dmatch();

main(argc,argv)
int argc; char *argv[];
{
int	i, j, l, n, nout, breaks, anyseq = -1;
char	code[20], line[100], secs[2000], chid = ' ';
Pdbentry_ *prot;
int	start, stop;
	if (argv[1][0]=='>') strcpy(code,argv[1]+1);
		else	     strcpy(code,argv[1]);
	l = strlen(code);
        if (l==5) { chid = code[l-1]; code[l-1] = 0; }
        strcpy(line,"/nimr/mb/data/pdb/brk/");
        strcat(line,code);
        strcat(line,".brk");
	prot = get_pdb(line,1,1);
	code[4] = 0;
	nout = 0;
	sscanf(argv[2],"%d", &start);
	sscanf(argv[3],"%d", &stop );
	Pi(start) Pi(stop) NL
	for (j=0; j<prot->Chainno; j++)
	{ int	nres = prot->Chains[j].Aano,
		nsec = prot->Chains[j].Secsno;
	  float asum, bsum, apct, bpct;
		if (nres==0) continue;
		if (prot->Chains[j].Chid != chid) continue;
		if ((anyseq>-1) && (prot->Chains[j].Chid == chid)) /* the next chain continues previous */
		{ int	difs = 0;			/* unless it is identical */
			for (i=0; i<min(30,nres); i++) { /* check first 30 res */
				if (prot->Chains[j].Atoms[i].Aa == prot->Chains[anyseq].Atoms[i].Aa) 
					continue;
				difs = 1;
				break;
			}
			if (difs) {
				for (i=0; i<nres; i++)
				{ Atom_ *Atoms = prot->Chains[j].Atoms;
					nout++;
					if (nout<start) continue;
					if (nout>stop ) continue;
					printf("ATOM%7d  CA ", nout);
            				printf("%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f",
                				Atoms[i].Alt, aa_code13(Atoms[i].Aa),
                				chid,Atoms[i].Resno,Atoms[i].Rid,
                				Atoms[i].X,Atoms[i].Y,Atoms[i].Z,
                				Atoms[i].Occu,Atoms[i].Bfact);
					printf("\n");
				}
				continue;
			}
		}
		if (anyseq>-1) break;
		n = 0;
		line[n] = ' ';
		for (i=0; i<strlen(prot->Compound); i++) {
			if (n==99) break;
			if (isspace(prot->Compound[i]) && isspace(line[n])) continue;
			n++;
			line[n] = prot->Compound[i];
		}
		line[n] = 0;
		printf("COMPND   %s\n",line);
		n = 0;
		line[n] = ' ';
		for (i=0; i<strlen(prot->Source); i++) {
			if (n==99) break;
			if (isspace(prot->Source[i]) && isspace(line[n])) continue;
			n++;
			line[n] = prot->Source[i];
		}
		line[n] = 0;
		printf("SOURCE   %s\n",line);
		printf("REMARK    ");
		chid = prot->Chains[j].Chid;
		printf("%s %c",prot->Pdbcode, chid);
		printf("%6.2f ",prot->Resol);
		printf(" %d-%d ",  prot->Chains[j].Atoms[0].Resno,
				   	prot->Chains[j].Atoms[nres-1].Resno);
		printf("(%d)",nres);
		printf("  PERCENT SEC.STR ");
		asum = bsum = -1.0;
		if (nsec) {
			asum = bsum = 0.0;
			for (i=0; i<nsec; i++)
			{ Secstr_ *ss = prot->Chains[j].Secs+i;
			  int	len = ss->End - ss->Beg + 1;
				if (ss->Sectype == HELIX) asum += (float)len;
				if (ss->Sectype == SHEET) bsum += (float)len;
			}
			asum = 100.0*asum/(float)nres;
			bsum = 100.0*bsum/(float)nres;
		}
		if (asum>99.9) asum = bsum = -1.0;
		if (bsum>99.9) asum = bsum = -1.0;
		printf("%5.1f%5.1f ", asum,bsum);
		breaks = consec_seen (prot->Chains+j,secs,&apct,&bpct);
		printf("(%5.1f%5.1f) ", apct,bpct);
		printf("\n");
		bsum = 0.0;
		for (i=0; i<nres; i++) bsum += prot->Chains[j].Atoms[i].Bfact;
		printf("REMARK    Mean B factor =%7.2f\n", bsum/(float)nres);
		for (i=0; i<nres; i++)
		{ Atom_ *Atoms = prot->Chains[j].Atoms;
			nout++;
			if (nout<start) continue;
			if (nout>stop ) continue;
			printf("ATOM%7d  CA ", nout);
            		printf("%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f",
                		Atoms[i].Alt, aa_code13(Atoms[i].Aa),
                		chid,Atoms[i].Resno,Atoms[i].Rid,
                		Atoms[i].X,Atoms[i].Y,Atoms[i].Z,
                		Atoms[i].Occu,Atoms[i].Bfact);
			printf("\n");
		}
		anyseq = j;
	}
	if (anyseq>-1) printf("TER\n");
}

consec_seen (chain,secs,pcta,pctb)
Chain_ *chain;
char *secs;
float	*pcta, *pctb;
{
int	breaks = 0;
int	i, j, len = chain->Aano;
float alpha[9] = {0.0,  3.82,  5.41,  5.04,  6.21,  8.66,  9.82, 10.52, 12.37 };
float betas[9] = {0.0,  3.75,  6.89, 10.44, 13.78, 17.29, 20.67, 24.16, 27.56 };
float noise = 0.35;
float     *psec, **mat;
	psec = (float*)alloca(sizeof(float)*(len+2));
	mat = (float**)alloca(sizeof(float*)*(len+2));
        for (i=1; i<=len; i++) {
		mat[i] = (float*)alloca(sizeof(float)*(len+2));
        	for (j=1; j<=len; j++) {
			mat[i][j] = cadist(chain->Atoms,i-1,j-1);
		}
		if (i==1) continue;
		if (mat[i][i-1] > 4.0) breaks++;
	}
	*pcta = *pctb = 0.0;
        for (i=1; i<=len; i++)
        { float rmsa, rmsb;
                rmsa = 1.25 * dmatch(mat,i,alpha,3,len);
                rmsb = 1.75 * dmatch(mat,i,betas,3,len);
                secs[i] = '-';
                psec[i] = -1.0;
                if (rmsa > noise) {
                        secs[i] = 'H';
                        psec[i] = rmsa;
                }
                if (rmsb > noise) {
                        secs[i] = 'E';
                        psec[i] = rmsb;
                }
                if (secs[i]=='E' && rmsa>rmsb) {
                        secs[i] = 'H';
                        psec[i] = rmsa;
                }
		if (secs[i]=='H') *pcta += 1.0;
		if (secs[i]=='E') *pctb += 1.0;
        }
	*pcta = *pcta * 100.0/(float)len;
	*pctb = *pctb * 100.0/(float)len;
	return breaks;
}

float cadist (a,i,j)
Atom_ *a;
int	i, j;
{
float	x = a[i].X - a[j].X,
	y = a[i].Y - a[j].Y,
	z = a[i].Z - a[j].Z;
	return sqrt(x*x + y*y + z*z);
}
	
float dmatch (mat, m, str, n, l)
float     **mat;
int     m, n, l;
float   *str;
{
float   score, sum, d, in, nn = (float)(n+n+1);
int     i, j;
        nn = nn*nn-nn;
        in = sum = 0.0;
        for (i=-n; i<=n; i++) {
                for (j=-n; j<=n; j++)
                { int   mi = m+i, mj = m+j;
                        if (i==j) continue;
                        if (mi<1 || mj<1) continue;
                        if (mi>l || mj>l) continue;
                        d = mat[mi][mj] - str[abs(i-j)];
                        sum += d*d;
                        in += 1.0;
                }
        }
        score = nn/(nn+sum);
        score = score*in/nn;
        return score;
}
