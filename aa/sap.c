#include <alloca.h>
#include "code/util.def"
#include "code/geom.def"
#include "aaincl/pdbprot.h"
#include "aaincl/matplot.h"
#define NALLOC 1000
#define NACID 30

typedef struct {
	Vec	v;
	float	d;
} Tri;

typedef struct {
	int	a, b;
	float	c;
} Pairs;

typedef struct {
	char	*res;
	float	*acc;
	Vec	*xyz;
	int	len;
} Seq;

float	compare();
float	recycle();
float	score();
float	match_list();
float	add_path();
float	get_path();
float	local_rms();

#define CYCLES 10
#define BIAS_WT 0.01
#define BIAS_DAMP 0.5
#define PATHSUM_WT 0.0
#define PATH_WT 1.0
#define N 8

int	seqmat[NACID][NACID];

main(argc,argv)
int argc; char *argv[];
{
Tri	**a[N];
Seq	seq[N];
char	id[N*4+1], file1[255], file2[255];
float	**scores;
int	i, j, half;
float	z = 1.0;
Pdbentry_ *prot1, *prot2;
	if (argc < 3) {
		strcpy(file1,"test1.pdb");
		strcpy(file2,"test2.pdb");
	} else {
		strcpy(file1,argv[1]);
		strcpy(file2,argv[2]);
	}
	Ps(file1) NL
	prot1 = get_pdb(file1,1,1);
	Ps(prot1->Compound) NL
	cones(prot1);
	Ps(file2) NL
	prot2 = get_pdb(file2,1,1);
	Ps(prot2->Compound) NL
	cones(prot2);
	matin("matrix/md.mat",seqmat);
	for (i=0; i<N/2; i++)
	{ int	ii = i*2,
		ij = i*2+1,
		rev = 0;
		if (i>1) rev = 1;
		protin(prot1,seq+ii,i,a+ii,z,rev),
		protin(prot2,seq+ij,i,a+ij,z,rev);
		z *= -1.0;
	}
	NLLL
	compare(a[1],a[0],seq+1,seq+0,1);
	strcpy(id,"A>+ B>+ A>- B>- A<+ B<+ A<- B<- ");
	scores = (float**)alloca(sizeof(float*)*N);
	NLL
	half = 0;
	for (i=0; i<N; i++) {
		scores[i] = (float*)alloca(sizeof(float)*N);
		id[i*4+3] = 0;
		printf("   %s   ", id+i*4);
		for (j=0; j<N; j++)
		{ float s; int k;
			if (half && i<j) continue;
			k = i+j+1;
			if (k-(2*(k/2))) {
				printf("         ");
				continue;
			}
			scores[i][j] = s = compare(a[i],a[j],seq+i,seq+j,0);
			printf("%9.1f",s);
		} NLL
	}
	printf("           ");
	for (i=0; i<N; i++) printf("   %s   ", id+i*4);
	NLLL
	stats(half,scores,N);
}

stats (half, data, n) float **data; int half, n; {
double	fn1, fn2, ave1, ave2, var1, var2;
float	sig1, sig2, dmax, dmin, smax, clear, score, noise;
int     i, j, k, n1, n2;
	n1 = n2 = 0;
        ave1 = ave2 = 0.0;
	smax = dmax = 0.0;
        for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			if (half && i<j) continue;
			k = i+j+1; 
			if (k-(2*(k/2))) continue;
			k = i+j-1;
			if ((i-j)*(i-j)==1 && !(k-(4*(k/4))) ) {
				if (dmax<data[i][j]) dmax = data[i][j];
				ave1 += data[i][j];
				n1++;
			} else {
				if (smax<data[i][j]) smax = data[i][j];
				ave2 += data[i][j];
				n2++;
			}
		}
	}
	fn1 = (double)n1; fn2 = (double)n2;
        ave1 /= fn1;
        ave2 /= fn2;
        var1 = var2 = 0.0;
	dmin = dmax;
	noise = 0.0;
        for (i=0; i<n; i++) { float d;
		for (j=0; j<n; j++) {
			if (half && i<j) continue;
			k = i+j+1;
			if (k-(2*(k/2))) continue;
			k = i+j-1;
			if ((i-j)*(i-j)==1 && !(k-(4*(k/4))) ) {
				if (dmin>data[i][j]) dmin = data[i][j];
				d = data[i][j] - ave1;
				var1 += d*d;
			} else {
				d = data[i][j] - ave2;
				var2 += d*d;
			}
			if (half || i<j) continue;
			d = data[i][j] - data[j][i];
			noise += d*d;
		}
	}
        var1 /= fn1; var2 /= fn2;
	stutest(ave1,ave2,var1,var2,n1,n2);
	sig1 = sqrt(var1); sig2 = sqrt(var2);
	noise = sqrt(2.0*noise/(fn1+fn2));
	score = (data[1][0] - ave2)/sig2;
	clear = (dmin-smax)/sig2;
	dmax = (dmax - ave2)/sig2;
	dmin = (dmin - ave2)/sig2;
	NLL
	printf("Alignment score = %5.3f StD(2) above mean controls\n", score);
	printf("Worst score = %5.3f StD(2) above best control\n", clear);
	printf("Max score = %5.3f, Min score = %5.3f\n", dmax, dmin);
	printf("StD(1) = %5.3f, StD(2) = %5.3f\n", sig1,sig2);
	if (half) return;
	printf("RMS alignment order noise = %5.3f\n\n", noise);
}

selsort (ac,bc) const void *ac, *bc;
{
Pairs   *a = (Pairs*)ac, *b = (Pairs*)bc;
        if (a->c < b->c) return 1;
        if (a->c > b->c) return -1;
        return 0;
}

float compare (a,b,seqa,seqb,print)
Tri	**a, **b; 
Seq	*seqa, *seqb;
int	print;
{
float	**bias, **sim, score;
int	i, j, m, n, cycle, nsel;
int	lena, lenb;
Pairs	*sel;
	lena = seqa->len;
	lenb = seqb->len;
        sel = (Pairs*)alloca(sizeof(Pairs)*lena*lenb);
        sim = (float**)alloca(sizeof(float*)*(lenb+2));
        bias = (float**)alloca(sizeof(float*)*(lenb+2));
        for (i=0; i<lenb+2; i++) {
                sim[i] = (float*)alloca(sizeof(float)*(lena+2));
                bias[i] = (float*)alloca(sizeof(float)*(lena+2));
                for (j=0; j<lena+2; j++) bias[i][j] = 0.0;
        }
        for (i=2; i<lenb-1; i++) {
                for (j=2; j<lena-1; j++)
		{ float acc, rms, seq, aa, ab, dab;
		  int	ra, rb;
			aa = seqa->acc[j];
			ab = seqb->acc[i];
			acc = fabs(aa-ab);
			acc = 1.0 * sqrt(acc);
			rms = 1.0 * local_rms(j,i,a,b);
			ra = seqa->res[j] - 'A';
			rb = seqb->res[i] - 'A';
			seq = 0.5*(float)seqmat[ra][rb];
			bias[i][j] = (10.0+seq)/(1.0+acc+rms);
		}
        }
	nsel = cycle = 0;
	score = recycle(cycle,seqa,seqb,bias,sel,bias,a,b,&nsel,0);
	for (cycle=1; cycle<=CYCLES; cycle++) {
		score = recycle(cycle,seqa,seqb,bias,sel,sim,a,b,&nsel,0);
	}
	score = recycle(10,seqa,seqb,bias,sel,sim,a,b,&nsel,print);
	return score;
}
/*
Pi(lena) Pi(lenb) NL
print_mat(10.0,bias,lena,lenb);
draw_mat("bias",1.0,bias,lena,lenb);
*/

float local_rms(m,n,a,b)
Tri     **a, **b;
{
float	sum, d;
int	i, j;
	sum = 0.0;
	for (i=-2; i<=2; i++) {
		for (j=2; j<=2; j++) {
			d = a[m+i][m+j].d - b[n+i][n+j].d;
			sum += d*d;
		}
	}
	return sqrt(sum/25.0);
}

float recycle (cycle,seqa,seqb,bias,sel,sim,a,b,nsel,print)
int	cycle;
Seq	*seqa, *seqb;
float	**bias, **sim;
Tri	**a, **b; 
Pairs	*sel;
int	*nsel, print;
{
int	**aln, len, i, j, n;
int     lena, lenb, onaln;
float	score, acut, cyc_no = (float)cycle;
	lena = seqa->len;
	lenb = seqb->len;
	if (*nsel) score_pair(cyc_no,bias,sel,sim,a,b,lena,lenb,*nsel);
	aln = (int**)alloca(sizeof(int*)*2); TEST(aln)
	for (i=0; i < 2; i++) {
		aln[i] = (int*)alloca(sizeof(int)*(lena+lenb)); TEST(aln[i])
	}
	score = get_path(aln,sim,lena,lenb,&len);
	if (print==1) { NL Pr(score) NL NL }
	if (print==2) print_mat(1.0,sim,lena,lenb);
	if (print==3) draw_mat("sim",100.0,sim,lena,lenb);
	for (i=1; i<=lenb; i++) {
		for (j=1; j<=lena; j++) {
			bias[i][j] *= BIAS_DAMP;
		}
	}
	for (i=len; i>0; i--) 
	{ int	a = aln[0][i],
		b = aln[1][i];
		bias[b][a] += log(1.0+sim[b][a]) * BIAS_WT;
	}
	for (i=0; i<10; i++)  normn(3.0,bias,lena,lenb);
	if (print==3) draw_mat("bias",1.0,bias,lena,lenb);
	n = 0;
	acut = 0.5 - 0.1*cyc_no;
	for (i=1; i<=lenb; i++) {
		if (seqb->acc[i] < acut) continue;
		for (j=1; j<=lena; j++) {
			if (seqa->acc[j] < acut) continue;
			sel[n].a = j;
			sel[n].b = i;
			sel[n].c = bias[i][j];
			n++;
		}
        }
	qsort(sel,n-1,sizeof(Pairs),selsort);
	if (*nsel) *nsel = (int)(0.05*sqrt((float)(lena*lenb))*cyc_no);
	*nsel += 20;
	if (!print) return score;
	onaln = check_sel(aln,sel,*nsel,len);
	printf("Percent sel on aln = %7.2f\n", 100.0*(float)onaln/(float)*nsel);
	printf("Percent aln in sel = %7.2f\n", 100.0*(float)onaln/(float)len);
	return score;
}
/*
printf("\nSEL %d\n", nsel);
for (i=0; i<nsel; i++) printf("%4d%4d%7.3f\n", sel[i].a,sel[i].b,sel[i].c);
*/

check_sel (aln,sel,nsel,naln) Pairs *sel; int **aln, nsel, naln; {
int	i, j, n = 0;
	for (i=0; i<nsel; i++)
	{ int	sela = sel[i].a,
		selb = sel[i].b;
		for (j=1; j<=naln; j++)
		{ int	alna = aln[0][j],
			alnb = aln[1][j];
			if (sela==alna && selb==alnb) {
				n++;
				break;
			}
		}
	}
	return n;
}

draw_mat (title, scale, a,m,n)
char *title; float scale, **a; int m, n;
{
long	id;
double	**b;
int	i, j;
	b = (double**)alloca(sizeof(double*)*(n+2));
	for (i=0; i<=n; i++) b[i] = (double*)alloca(sizeof(double)*(m+2));
	for (i=1; i<=n; i++) for (j=1; j<=m; j++) b[i-1][j-1] = a[i][j];
	id = init_matplot(n-1,m-1,title,DXORIG+200,DYORIG);
	display_mat(id,b,n-1,m-1,0.0,scale);
}

trace_mat (mat,lena,lenb) int **mat, lena, lenb; {
int	i, j;
	NL Pi(lena) Pi(lenb) NL
	for (i=1; i<=lenb; i++) { 
		for (j=1; j<=lena; j++) printf("%3d", mat[i][j]); NL
	} NL
}

print_mat (scale, mat,lena,lenb) float scale, **mat; int lena, lenb; {
int	i, j;
	NL Pi(lena) Pi(lenb) NL
	for (i=1; i<=lenb; i++) { 
		for (j=1; j<=lena; j++)
		{ char	c;
	  	  float sij = mat[i][j]*scale;
			if (sij<10.0) c = '0'+(int)sij;
				else c = 'A'+(int)(sij*0.1)-1;
				if (sij > 260.0) c = '*';
				if (c=='0') c = '.';
			printf("%c", c);
		} NL
	} NL
}

score_pair (dif_wt,bias,sel,sim,a,b,la,lb,nsel)
float	dif_wt;
float	**bias, **sim;
Tri	**a, **b; 
Pairs	*sel;
int	la, lb, nsel;
{
int	i, j, k, l;
        for (i=0; i<lb+2; i++) {
                for (j=0; j<la+2; j++) sim[i][j] = 0.0;
        }
	for (i=0; i<nsel; i++) {
		score(dif_wt,bias,sim,a,b,sel[i].a,sel[i].b,la,lb);
	}
}

float score (dif_wt,bias,sim,a,b,m,n,la,lb)
float	dif_wt;
float	**bias, **sim;
Tri	**a, **b; 
int	m, n, la, lb;
{
int	i, j, k, na, nb, minlen;
float   **smn, path_score;
        smn = (float**)alloca(sizeof(float*)*(lb+2));
        for (i=0; i<lb+2; i++) {
                smn[i] = (float*)alloca(sizeof(float)*(la+2));
                for (j=0; j<la+2; j++) smn[i][j] = 0.0;
        }
	for (i=1; i<=lb; i++) {
		if (n==i) continue;
                for (j=1; j<=la; j++) { float d, v;
			if (m==j) continue;
			d = vddif(a[m][j].v,b[n][i].v);
			smn[i][j] = 10.0/(1.0+d)
				  + 0.2*((10.0-dif_wt)*bias[i][j]);
		}
	}
	path_score = add_path(sim,smn,la,lb,m,n);
	sim[n][m] +=  PATHSUM_WT * path_score;
}
/*
			d = a[m][j].d - b[n][i].d;
			v = vdot(a[m][j].v,b[n][i].v)+1.0;
			smn[i][j] = 10.0*(1.0-exp(-v*v))/(1.0+d*d)
				  + (10.0-dif_wt)*bias[i][j];
*/

float add_path (sim,smn,na,nb,m,n) 
float	**sim, **smn;
int	na, nb, m, n;
{
int	**aln, len, i;
float	**s, score = 0.0;
	aln = (int**)alloca(sizeof(int*)*2); TEST(aln)
	for (i=0; i < 2; i++) {
		aln[i] = (int*)alloca(sizeof(int)*(na+nb)); TEST(aln[i])
	}
	if (m>1 && n> 1) {
		score += get_path(aln,smn,m-1,n-1,&len);
		for (i=len; i>0; i--) 
		{ int	a = aln[0][i],
			b = aln[1][i];
			sim[b][a] += smn[b][a] * PATH_WT;
		}
	}
	if (m<nb && n<na) {
		s = (float**)alloca(sizeof(float*)*(nb+2)); TEST(s)
		for (i=n; i<nb+2; i++) s[i-n] = smn[i]+m;
		score += get_path(aln,s,na-m,nb-n,&len);
		for (i=len; i>0; i--) 
		{ int	a = aln[0][i],
			b = aln[1][i];
			sim[b+n][a+m] += s[b][a] * PATH_WT;
		}
	}
	return score;
}

float get_path (aln,sim,na,nb,length) 
float	**aln, **sim;
int	na, nb, *length;
{
float	**mat;
int	**ptr, i, j, k;
float	score, *colmax, rowmax;
int	*maxcol, maxrow, maxi, maxj, len;
int	naa = na+2, nbb = nb+2, now;
	mat = (float**)alloca(sizeof(float*)*2); TEST(mat)
	for (i=0; i<2; i++) {
		 mat[i] = (float*)alloca(sizeof(float)*naa); TEST(mat[i])
	}
	ptr = (int**)alloca(sizeof(int*)*nbb); TEST(ptr)
	for (i=0; i<nbb; i++) {
		ptr[i] = (int*)alloca(sizeof(int)*naa); TEST(ptr[i])
	}
	colmax = (float*)alloca(sizeof(float)*(naa)); TEST(colmax)
	maxcol = (int*)alloca(sizeof(int)*(naa)); TEST(maxcol)
	for (i=0; i<naa; i++) {
		maxcol[i] = 0;
		colmax[i] = mat[0][i] = mat[1][i] = -1.0;
	}
	score = 0.0;
	now = 1;
	for (i=1; i<nbb; i++) {
		rowmax = -1.0;
		for (j=1; j<naa; j++)
		{ float dig, col, row, max;
		  int	cop, rop, top;
			rop = cop = top = 0;
			row = col = max = 0.0;
			if (j>na || i>nb) mat[now][j] = 0.0;
				     else mat[now][j] = sim[i][j];
			max = dig = mat[!now][j-1];
			if (colmax[j-1] > dig) {
				col = colmax[j-1];
				cop = i-maxcol[j-1]-1;
			} else {
				colmax[j-1] = dig;
				maxcol[j-1] = i-1;
			}
			if (rowmax > dig) {
				row = rowmax;
				rop = -(j-maxrow-1);
			} else {
				rowmax = dig;
				maxrow = j-1;
			}
                        if (row > max) { max = row; top = rop; }
                        if (col > max) { max = col; top = cop; }
			mat[now][j] += max;
			ptr[i][j] = top;
			if (mat[now][j] > score) {
				score = mat[now][j];
				maxi = i;
				maxj = j;
			}
		}
		now = !now;
	}
	*length = 0;
	if (score > 0.1) *length = trace(sim,ptr,aln,0,maxi,maxj);
	return score;
}

trace (s,p,a,n,i,j) float **s; int **p, **a, n, i, j;
{
	if (s[i][j] < 0.0) return n;
	n++;
	a[0][n] = j;
	a[1][n] = i;
	if (i<=1 || j<=1) return n;
	if (p[i][j] > 0) i -= p[i][j];
	if (p[i][j] < 0) j += p[i][j];
	return trace(s,p,a,n,--i,--j);
}

protin (prot,seq,id,m,z,flip)
Pdbentry_ *prot;
Seq *seq; Tri ***m; float z; int flip, id;
{
FILE	*pdb;
int	i, j, len;
Tri	**mat;
	len = copyca(prot->Chains,seq,flip,z);
        mat = (Tri**)malloc(sizeof(Tri*)*(len+2));
        for (i=0; i<=len+1; i++) {
                mat[i] = (Tri*)malloc(sizeof(Tri)*(len+2));
        }
	set_vect(seq->xyz,mat,len);
	*m = mat;
	return len;
}

set_vect (a,m,l) Vec *a; Tri **m; int l; {
int	i, j;
Mat	frame;
	for (i=1; i<=l; i++) {
		setframe(a[i-1],a[i],a[i+1],&frame);
		for (j=1; j<=l; j++) { Vec s, t;
			m[i][j].d = vdif(a[i],a[j]);
			vinit(&(m[i][j].v));
			if (i==j) continue;
			vsub(a[j],a[i],&s);
			VmulM(&frame,s,&(m[i][j].v));
		}
	}
}
/*
			vnorm(&(m[i][j].v));
*/

set_dist (a,m,l) Vec *a; Tri **m; int l; {
int	i, j;
	for (i=0; i<=l; i++) {
		for (j=0; j<=l; j++) {
			m[i][j].d = vdif(a[i],a[j]);
		}
	}
}

extend (res,i,j,k,new)
Vec	*res;
int	i, j, k, new;
{
	Vec	m, v;
	vave(res[j],res[k],&m);
	vsub(m,res[i],&v);
	vadd(m,v,&res[new]);
}
 
copyca (pdb,s,flip,z)
Chain_  *pdb;
Seq	*s;
int	flip;
float	z;
{	int	i, n;
	char	*seq;
	Vec	*xyz;
	float	*acc;
	n = pdb->Aano;
	seq = (char*)malloc(sizeof(char)*(n+3));
	acc = (float*)malloc(sizeof(float)*(n+3));
	xyz = (Vec*)malloc(sizeof(Vec)*(n+3));
	for (i=0; i<n; i++) {
		xyz[i+1].x = pdb->Atoms[i].X;
		xyz[i+1].y = pdb->Atoms[i].Y;
		xyz[i+1].z = pdb->Atoms[i].Z;
		acc[i+1] = pdb->Atoms[i].Bfact;
		seq[i+1] = pdb->Atoms[i].Aa;
		if (seq[i+1]<'A' || seq[i+1]>'Z') {
			printf("*NB* funny aa = %c\n", seq[i+1]);
			seq[i+1] = 'X';
		}
	}
	seq[0] = 'n';
        extend(xyz,3,2,1,0);    
        extend(xyz,n-2,n-1,n,n+1);
	seq[n+1] = 'c';
	seq[n+2] = 0;
	for (i=0; i<=n+1; i++) xyz[i].z *= z;
	if (flip) flipseq(xyz,seq,acc,n);
	s->res = seq;
	s->acc = acc;
	s->xyz = xyz;
	s->len = n;
	return n;
}

flipseq (xyz,seq,acc,n) Vec *xyz; char *seq; float *acc; int n;
{
int	i;
	for (i=0; i<=n/2; i++)
	{ Vec r; char c; float a;
	  int j = n+1-i;
		r = xyz[i]; xyz[i] = xyz[j]; xyz[j] = r;
		c = seq[i]; seq[i] = seq[j]; seq[j] = c;
		a = acc[i]; acc[i] = acc[j]; acc[j] = a;
	}
}
 
getca (res,pdb)
Vec    *res;
FILE	*pdb;
{	int	i = 1;
	char	line[225], junk[30];
        while(!feof(pdb)) {
		read_line(pdb,line);
		if (!strncmp(line,"TER",3)) break;
		if (strncmp(line,"ATOM",4)) continue;
		if (strncmp(line+13,"CA ",3)) continue;
		sscanf(line,"%30c %f%f%f",
                       	junk, &res[i].x, &res[i].y, &res[i].z);
		i++;
	}
	i--;
        extend(res,3,2,1,0);    
        extend(res,i-2,i-1,i,i+1);
	return i;
}
 
putpdb (res,out,len)
Vec    *res;
FILE    *out;
int     len;
{       int     i = 0, n = 0;
        for (i=1; i<=len; i++) {
                fprintf(out,"ATOM%7d  CA  GLY%6d     %7.3f %7.3f %7.3f   0.0   0.0\n",
                        i, i, res[i].x, res[i].y, res[i].z);
        }
        fprintf(out,"TER\n");
}

setframe (a, b, c, frame)
    Vec a, b, c;
    Mat *frame;
{
    int    i;
    Vec    x, y, z ;
	vsub(c,a,&x);
	vave(c,a,&c);
	vsub(c,b,&y);
	vprod(y,x,&z);
	vprod(z,x,&y);
	vnorm(&x);
	vnorm(&y);
	vnorm(&z);
	VtoM(x,y,z,frame);
}

normn (sigcut, data, m, n) float sigcut, **data; int m, n; {
float   d, fn, dmax,
        ave, var, sig;
int     i, j, k, mods;
        fn = (float)(m*n);
        ave = 0.0;
        for (i=1; i<=n; i++) for (j=1; j<=m; j++) ave += data[i][j];
        ave /= fn;
        var = 0.0;
        for (i=1; i<=n; i++) for (j=1; j<=m; j++) var += data[i][j]*data[i][j];
        var /= fn;
        sig = sqrt(var);
        mods = 0;
	dmax = 1.0;
        for (i=1; i<=n; i++) for (j=1; j<=m; j++) {
                data[i][j] /= sig;
                if (data[i][j] > sigcut) {
                        data[i][j] =  sigcut + 0.5*(data[i][j]-sigcut);
                        mods++;
                }
		if (data[i][j] > dmax) dmax = data[i][j];
        }
        for (i=1; i<=n; i++) for (j=1; j<=m; j++) data[i][j] /= dmax;
        return mods;
}

matin(file,mat)
        char    *file;
        int     mat[NACID][NACID];
{
        int     i, j, mat_const;
        char    acid[NACID], c;
        FILE    *mat_file;

        mat_file = fopen(file,"r");
        while( c = getc(mat_file), c != '\n' ) putchar(c); NL
        fscanf(mat_file,"%s\n",acid);
        printf("%s\n",acid);
        fscanf(mat_file,"%d\n",&mat_const);
        printf("matrix constant = %d\n",mat_const);
        for( i = 0; acid[i]; i++ ) 
        {       int     ai = acid[i]-'A';
                for( j = 0; acid[j]; j++ ) 
                {       int aj = acid[j]-'A';
                        fscanf(mat_file,"%d",&mat[ai][aj]);
                        mat[ai][aj] += mat_const;
                }
        }
}

super (seqa, seqb, sim, aln, len)
Seq	*seqa, *seqb;
float	**sim;
int	**aln, len;
{
double	*wt;
double	**va, **vb;
int	na = seqa->len;
int	nb = seqb->len;
int i, j;
	wt = (double*)alloca(sizeof(double)*len);
	va = (double**)alloca(sizeof(double*)*len);
	vb = (double**)alloca(sizeof(double*)*len);
	for (i=len; i>0; i--) 
	{ int	a = aln[0][i],
		b = aln[1][i];
	  char  ra = seqa->res[a],
		rb = seqb->res[b],
		aa1, aa2, ab1, ab2;
		aa1 = aa2 = ab1 = ab2 = ' ';
		if (seqa->acc[a]>0.0) aa1 = '*';
		if (seqb->acc[b]>0.0) ab1 = '*';
		if (seqa->acc[a]>0.5) aa2 = '*';
		if (seqb->acc[b]>0.5) ab2 = '*';
		printf("%c%c%c %4d %5.1f%4d %c%c%c\n",
			aa2,aa1,ra,a,sim[b][a],b,rb,ab1,ab2);
		va[i] = (double*)alloca(sizeof(double)*3);
	 	va[i][0] = seqa->xyz[a].x;
	 	va[i][1] = seqa->xyz[a].y;
	 	va[i][2] = seqa->xyz[a].z;
		vb[i] = (double*)alloca(sizeof(double)*3);
	 	vb[i][0] = seqb->xyz[a].x;
	 	vb[i][1] = seqb->xyz[a].y;
	 	vb[i][2] = seqb->xyz[a].z;
		wt[i] = sim[b][a];
	}
	supermac(wt,va,vb,len);
	
}
