/*
cc -O topol.c -o topol util/wt/util.o util/wt/geom.o util/wt/sort.o -lm -m32

*/
#include "util/wt/incl/util.h"
#include "util/wt/incl/geom.h"

#define NALLOC 1500
#define BUMP 0.5
#define CYCLES 500

Vec    res[NALLOC];
Vec    ret[NALLOC];
Vec    ray[NALLOC];

float	bond;
float	packs();

main(argc,argv)
int argc; char *argv[];
{
int	i, j, n, len, bumps;
char	line[225];
FILE	*pdb, *out;
int	beg, end, cycles, check;
float	r, rmax, bumper;
Vec    	v, cog;
	pdb = fopen(argv[1],"r");
	beg = end = 0;
	if (argc > 2) {
		sscanf(argv[2],"%d", &beg);
		sscanf(argv[3],"%d", &end);
		Pi(beg) Pi(end) NL
	}
	cycles = CYCLES;
	if (argc > 4) sscanf(argv[4],"%d", &cycles);
	check = 1;
	if (argc > 5) sscanf(argv[5],"%d", &check);
	bumper = BUMP;
	if (argc > 6) sscanf(argv[6],"%f", &bumper);
	Pi(check) Pi(cycles) Pr(bumper) NL
	len = getca(ray,pdb);
	if (!end) { beg = 1; end = len; }
	n = 0;
	for (i=beg; i<=end; i++) {
		n++;
		vcopy(ray[i],&res[n]);
		vcopy(ray[i],&ret[n]);
	}
	len = n;
	Pi(len) NL
	out = fopen("start.out","w");
	putpdb(res,out,len,1.0);
	fclose(out);
	out = fopen("soap.out","w");
	len = smooth(out,res,ret,len,cycles,check,bumper);
	Pi(len) NL
	fclose(out);
	out = fopen("fold.out","w");
	putpdb(res,out,len,1.0);
	if (len==2) exit(1);
        vinit(&cog);
        for (i=2; i<len; i++) vsum(res[i],&cog);
        vdiv(&cog,(float)(len-2));
        for (i=1; i<=len; i++) vsub(res[i],cog,ret+i);
	rmax = 0.0;
        for (i=2; i<len; i++) {
		r = vmod(ret[i]);
		if (r > rmax) rmax = r;
	}
	Pr(rmax) NL
	vcopy(ret[1],ret);
	vcopy(ret[len],ret+len+1);
	vsub(ret[1],ret[2],&v);
	vnorm(&v); vmul(&v,rmax*2.0); vadd(ret[2],v,ret+1);
	vnorm(&v); vmul(&v,rmax*3.0); vadd(ret[2],v,ret+0);
	vsub(ret[len],ret[len-1],&v);
	vnorm(&v); vmul(&v,rmax*2.0); vadd(ret[len-1],v,ret+len);
	vnorm(&v); vmul(&v,rmax*3.0); vadd(ret[len-1],v,ret+len+1);
	len += 2;
	rmax = 0.0;
        for (i=1; i<len; i++) {
        	for (j=i+1; j<=len; j++) {
			r = vdif(ret[i],ret[j]);
			if (r > rmax) rmax = r;
		}
	}
	Pr(rmax) NL
	for (i=1; i<=len; i++) vcopy(ret[i-1],&res[i]);
	for (i=1; i<=len; i++) vmul(&res[i],10.0/rmax);
	for (i=1; i<=len; i++) vcopy(res[i],&ret[i]);
	out = fopen("knot.out","w");
	putpdb(res+1,out,len-2,1.0);
}

smooth (out,res,ret,len,cycles,check,bumper)
FILE	*out;
Vec	*res, *ret;
int	len, cycles,check;
float	bumper;
{
int	i, j, n, last, bumps, bumpsum;
float	bumpmin, bumpj, bj;
float	sca=0.2;
	bj = bumper/(float)CYCLES;
	n = bumps = bumpsum = 0;
	bumpmin = bumper/1000.0;
	Pi(n) Pi(bumps) Pi(len) Pr(bumper) Pr(bj) NL
	last = len;
	bumpj = bumper;
	for (j=0; j<cycles; j++) { Vec mid;
		n++;
		bumps = move(res,ret,len,bumpj,check);
		bumpsum += bumps;
		while (1) { int las;
			las = len;
			len = shrink(ret,len,bumpj,check);
			if (len==las) break;
		}
		if (last==len) bumpj = bumpj - bj;
		if (bumpj<bumpmin) bumpj = bumpmin;
		Pi(n) Pr(bumpj) Pi(bumps) Pi(len) NL
		for (i=1; i<=len; i++) vcopy(ret[i],&res[i]);
		if (len==2) break;
		if (!(j-(j/2)*2)) putpdb(res,out,len,sca);
		last = len;
	}
	putpdb(res,out,len,sca);
	printf("%d cycles with %d bumps\n", n, bumpsum);
	return len;
}

move (a,b,len,bump,check)
Vec     *a, *b;
int     len;
float	bump;
int	check;
{
float	bbump2 = bump*bump*2.0;
int     i, bumps = 0;
	Pr(bump) NL
	for (i=2; i<len; i++) { Vec mid;
		vave(b[i-1],a[i+1],&mid);
		vave(a[i],mid,&b[i]);
		if (!check) continue;
		if (!stuck(a,b,i,len,bump)) continue;
		vcopy(a[i],&b[i]);
		bumps++;
	}
	return bumps;
}

bumps (a,b,c,d,p,q,bump)
Vec     a, b, c, d, p, q;
float   bump;
{
/* true if the move b -> c bumps or crosses p-q */
/* the new point c must lie on and within abd.  */
        if (clash(a,c,p,q,bump)) return 1;
        if (clash(d,c,p,q,bump)) return 1;
        if (line2tri(a,b,d,p,q)) {
                if (line2tri(a,b,c,p,q)) return 1;
                if (line2tri(d,b,c,p,q)) return 1;
        }
        return 0;
}

stuck (a,b,n,len,bump)
Vec     *a, *b;
int     n, len;
float   bump;
{
int     i, j;
        for (i=1; i<n-2; i++) {
		if (bumps(b[n-1],a[n],b[n],a[n+1],b[i],b[i+1],bump)) return 1;
	}
        for (i=n+2; i<len; i++) {
		if (bumps(b[n-1],a[n],b[n],a[n+1],a[i],a[i+1],bump)) return 1;
	}
        return 0;
}

shrink (a,len,bump,check)
Vec     *a;
int     len;
float	bump;
int	check;
{
int     i, j, remove, skip[10000];
        for (i=1; i<=len; i++) skip[i] = 0;
        for (i=2; i<len; i++) { Vec b, c;
		vsub(a[i-1],a[i],&b);
		vnorm(&b);
		vsub(a[i+1],a[i],&c);
		vnorm(&c);
		if (check >= 0) {
			if (vdot(b,c) > -0.99
			 && vdif(a[i-1],a[i+1]) > bump*2.0 ) continue;
		} else {
			 if (vdif(a[i-1],a[i+1]) > bump*2.0 ) continue;
		}
		remove = 1;
        	for (j=1; j<len; j++) {
			if (!check) continue;
			if (j>i-3 && j<i+2) continue;
			if (clash(a[i-1],a[i+1],a[j],a[j+1],bump)) {
				remove = 0;
				break;
			}
			if (line2tri(a[i-1],a[i],a[i+1],a[j],a[j+1])) {
				remove = 0;
				break;
			}
		}
		if (remove) skip[i] = 1;
		break;
	}
	j = 1;
        for (i=1; i<=len; i++) if (!skip[i]) { a[j++] = a[i]; }
	len = j-1;
	return len;
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
 
putpdb (res,out,len,sca)
Vec    *res;
FILE    *out;
int     len;
float	sca;
{       int     i = 0, n = 0;
        for (i=1; i<=len; i++) {
                	fprintf(out,"ATOM%7d  CA  GLY%6d     %7.3f %7.3f %7.3f   0.0   0.0\n",
                        i, i, sca*res[i].x, sca*res[i].y, sca*res[i].z);
        }
        fprintf(out,"TER\n");
}

clash (a,b,c,d,s) Vec a,b,c,d; float s; {
float   ss;
        if (line2line(a,b,c,d,s)) return 1;
        if (line2dot(a,b,c,s)) return 1;
        if (line2dot(a,b,d,s)) return 1;
        if (line2dot(c,d,a,s)) return 1;
        if (line2dot(c,d,b,s)) return 1;
        ss = s*s;
        if (vddif(a,c) < ss) return 1;
        if (vddif(a,d) < ss) return 1;
        if (vddif(b,c) < ss) return 1;
        if (vddif(b,d) < ss) return 1;
        return 0;
}
