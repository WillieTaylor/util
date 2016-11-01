/*
cc -O super.c -o super util/aa/bestrot.o util/wt/util.o util/wt/geom.o util/wt/sort.o util/aa/pdbprot.o util/aa/matrix.o util/aa/siva.o util/aa/ql.o util/eigen.o -lm -m32

*/
#include <alloca.h>
#include "util/wt/incl/util.h"
#include "util/wt/incl/geom.h"
#include "util/aa/incl/pdbprot.h"
#include "util/aa/incl/matrix.h"

double	supermac();

main(argc,argv)
int argc; char *argv[];
{
char	file1[255], file2[255];
int	i, j, len;
float	rms;
Pdbentry_ *prot1, *prot2;
double  *ww, *ac, *bc;
double  **va, **vb;
Sqmat_  rot;
	if (argc < 3) {
		printf("usage:  sap file1.pdb file2.pdb\n");
		exit(1);
	}
	prot1 = get_pdb(argv[1],1,1);
	Ps(prot1->Compound) Pi(prot1->Chains[0].Aano) NL
	prot2 = get_pdb(argv[2],1,1);
	Ps(prot2->Compound) Pi(prot2->Chains[0].Aano) NL
	len = prot1->Chains[0].Aano;
	if (prot2->Chains[0].Aano < prot2->Chains[0].Aano) len = prot2->Chains[0].Aano;
        ww = (double*)alloca(sizeof(double)*len);
        va = (double**)alloca(sizeof(double*)*len);
        vb = (double**)alloca(sizeof(double*)*len);
	for (i=0; i<len; i++) {
		ww[i] = 1.0;
                va[i] = (double*)alloca(sizeof(double)*3);
                va[i][0] = prot1->Chains[0].Atoms[i].X;
                va[i][1] = prot1->Chains[0].Atoms[i].Y;
                va[i][2] = prot1->Chains[0].Atoms[i].Z;
                vb[i] = (double*)alloca(sizeof(double)*3);
                vb[i][0] = prot2->Chains[0].Atoms[i].X;
                vb[i][1] = prot2->Chains[0].Atoms[i].Y;
                vb[i][2] = prot2->Chains[0].Atoms[i].Z;
	}
	rot = alloc_sqmat(3);
	ac = (double*)alloca(sizeof(double)*3);
	bc = (double*)alloca(sizeof(double)*3);
	rms = supermac(ww,va,vb,len,ac,bc,rot);
	Pr(rms) NL
}
