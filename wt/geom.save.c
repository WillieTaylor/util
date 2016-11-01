#include "incl/util.h"
#include "incl/geom.h"

void vinit (	Vec *c )
{
	c->x = 0.0;
	c->y = 0.0;
	c->z = 0.0;
}

void vcopy (	Vec b, Vec *c )
{
	c->x = b.x;
	c->y = b.y;
	c->z = b.z;
}

void vunit (	Vec b, Vec *c )
{	float	d = DIST3v((b));
	c->x = b.x/d;
	c->y = b.y/d;
	c->z = b.z/d;
}

void vnorm (	Vec *c )
{	float	d = DIST3v((*c));
	c->x = c->x/d;
	c->y = c->y/d;
	c->z = c->z/d;
}

void vave ( Vec a, Vec b, Vec *c )
{
	c->x = 0.5 * (a.x + b.x);
	c->y = 0.5 * (a.y + b.y);
	c->z = 0.5 * (a.z + b.z);
}	 

void vsum ( Vec a, Vec *c )
{
	c->x += a.x;
	c->y += a.y;
	c->z += a.z;
}	 

void vadd ( Vec a, Vec b, Vec *c )
{
	c->x = a.x + b.x;
	c->y = a.y + b.y;
	c->z = a.z + b.z;
}	 

void vsub ( Vec a, Vec b, Vec *c)
{
	c->x = a.x - b.x;
	c->y = a.y - b.y;
	c->z = a.z - b.z;
}	 

void vmul ( Vec *c, float s )
{
	c->x = (c->x)*s;
	c->y = (c->y)*s;
	c->z = (c->z)*s;
}

void vdiv ( Vec *c, float s )
{
	c->x = (c->x)/s;
	c->y = (c->y)/s;
	c->z = (c->z)/s;
}

void vmul2 ( Vec a, float s, Vec *c )
{
	c->x = a.x*s;
	c->y = a.y*s;
	c->z = a.z*s;
}

void vdiv2 ( Vec a, float s, Vec *c )
{
	c->x = a.x/s;
	c->y = a.y/s;
	c->z = a.z/s;
}

void vrset ( Vec *c, float s )
{
	c->x = s*(2.0*drand48()-1.0);
	c->y = s*(2.0*drand48()-1.0);
	c->z = s*(2.0*drand48()-1.0);
}

void vradd ( Vec *c, float s )
{
	c->x = (c->x)+s*(2.0*drand48()-1.0);
	c->y = (c->y)+s*(2.0*drand48()-1.0);
	c->z = (c->z)+s*(2.0*drand48()-1.0);
}

float	vdif ( Vec a, Vec b )
{	Vec	c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	return sqrt(c.x*c.x + c.y*c.y + c.z*c.z);
}	 

float	vddif ( Vec a, Vec b )
{	Vec	c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;
	return (c.x*c.x + c.y*c.y + c.z*c.z);
}	 

float	vdad ( Vec a, Vec b )
{	Vec	c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;
	return sqrt(c.x*c.x + c.y*c.y + c.z*c.z);
}	 

float	vddad ( Vec a, Vec b )
{	Vec	c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;
	return (c.x*c.x + c.y*c.y + c.z*c.z);
}	 

float	vsqr ( Vec a )
{
	return (a.x*a.x + a.y*a.y + a.z*a.z);
}	 

float	vmod ( Vec a )
{
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}	 

float	vdot ( Vec a, Vec b )
{
	return (a.x*b.x + a.y*b.y + a.z*b.z);
}	 

void vprod ( Vec a, Vec b, Vec *c )
{
      c->x = a.y*b.z - b.y*a.z;
      c->y = a.z*b.x - b.z*a.x;
      c->z = a.x*b.y - b.x*a.y;
}

float	vtri ( Vec a, Vec b, Vec c )
{	Vec	d;
	vprod(a,b,&d);
	return vdot(c,d);
}	 

float	pdotp ( Vec a, Vec b, Vec c, Vec d )
{	Vec	x, y;
	vsub(b,a,&x);
	vsub(d,c,&y);
	return vdot(x,y);
}	 

float	pvol ( Vec a, Vec b, Vec c, Vec d )
{	Vec	x, y, z;
	vsub(b,a,&x);
	vsub(d,c,&y);
	vsub(c,b,&z);
	return vtri(x,y,z);
}	 

float	phand ( Vec a, Vec b, Vec c, Vec d )
{	Vec	x, y, z;
	vsub(b,a,&x);
	vsub(d,c,&y);
	vsub(c,b,&z);
	vnorm(&z);
	return vtri(x,y,z);
}	 

void identM ( Mat *M ) {
        M->A.x = 1.0; M->A.y = 0.0; M->A.z = 0.0;
        M->B.x = 0.0; M->B.y = 1.0; M->B.z = 0.0;
        M->C.x = 0.0; M->C.y = 0.0; M->C.z = 1.0;
}

void VtoM ( Vec a, Vec b, Vec c, Mat *M ) {
        M->A.x = a.x; M->A.y = a.y; M->A.z = a.z;
        M->B.x = b.x; M->B.y = b.y; M->B.z = b.z;
        M->C.x = c.x; M->C.y = c.y; M->C.z = c.z;
}

void MtoT ( Mat *M, Mat *T ) {
        T->A.x = M->A.x; T->A.y = M->B.x; T->A.z = M->C.x;
        T->B.x = M->A.y; T->B.y = M->B.y; T->B.z = M->C.y;
        T->C.x = M->A.z; T->C.y = M->B.z; T->C.z = M->C.z;
}

void Mprint ( Mat *m ) {
Mat M = *m;
        printf("\n");
	printf("%9.6f %9.6f %9.6f\n", M.A.x,M.A.y,M.A.z);
	printf("%9.6f %9.6f %9.6f\n", M.B.x,M.B.y,M.B.z);
	printf("%9.6f %9.6f %9.6f\n", M.C.x,M.C.y,M.C.z);
        printf("\n");
}

void MmulM ( Mat *p, Mat *q, Mat *R ) {
Mat P = *p, Q = *q;
Vec a = P.A, b = P.B, c = P.C;
Vec A = Q.A, B = Q.B, C = Q.C;

	R->A.x = a.x*A.x + b.x*A.y + c.x*A.z;
        R->B.x = a.x*B.x + b.x*B.y + c.x*B.z;
        R->C.x = a.x*C.x + b.x*C.y + c.x*C.z;

	R->A.y = a.y*A.x + b.y*A.y + c.y*A.z;
        R->B.y = a.y*B.x + b.y*B.y + c.y*B.z;
        R->C.y = a.y*C.x + b.y*C.y + c.y*C.z;

	R->A.z = a.z*A.x + b.z*A.y + c.z*A.z;
        R->B.z = a.z*B.x + b.z*B.y + c.z*B.z;
        R->C.z = a.z*C.x + b.z*C.y + c.z*C.z;
}

void Minv ( Mat *m, Mat *W, float d )
{
Mat M = *m;
Vec a = M.A, b = M.B, c = M.C;

	W->A.x =  (b.y*c.z - b.z*c.y)/d;
        W->A.y = -(a.y*c.z - a.z*c.y)/d;
	W->A.z =  (a.y*b.z - a.z*b.y)/d;

	W->B.x = -(b.x*c.z - b.z*c.x)/d;
        W->B.y =  (a.x*c.z - a.z*c.x)/d;
        W->B.z = -(a.x*b.z - a.z*b.x)/d;

	W->C.x =  (b.x*c.y - b.y*c.x)/d;
        W->C.y = -(a.x*c.y - a.y*c.x)/d;
        W->C.z =  (a.x*b.y - a.y*b.x)/d;
}

float Mdet ( Mat *m ) {
Mat M = *m;
Vec a = M.A, b = M.B, c = M.C;
	return a.x * (b.y*c.z - b.z*c.y)
             - b.x * (a.y*c.z - a.z*c.y)
             + c.x * (a.y*b.z - a.z*b.y);
}

void MmulV ( Mat *m, Vec d, Vec *e ) {
Mat M = *m;
Vec a = M.A, b = M.B, c = M.C;
        e->x = d.x*a.x + d.y*b.x + d.z*c.x;
        e->y = d.x*a.y + d.y*b.y + d.z*c.y;
        e->z = d.x*a.z + d.y*b.z + d.z*c.z;
}

void VmulM ( Mat *m, Vec d, Vec *e ) {
Mat M = *m;
Vec a = M.A, b = M.B, c = M.C;
	e->x = d.x*a.x + d.y*a.y + d.z*a.z;
	e->y = d.x*b.x + d.y*b.y + d.z*b.z;
	e->z = d.x*c.x + d.y*c.y + d.z*c.z;
}

void Mrot ( char axis, Mat *R, float r ) {
Vec	A = R->A, B = R->B, C = R->C;
Vec	a, b, c;
float	sr, cr;
	a.x = a.y = a.z = 0.0;
	b.x = b.y = b.z = 0.0;
	c.x = c.y = c.z = 0.0;
	cr = cos(r); sr = sin(r);
	if (axis=='X') {   a.x = 1.0;
		b.y =  cr; b.z = -sr;
		c.y =  sr; c.z =  cr;
	}
	if (axis=='Y') {   b.y = 1.0;
		a.x =  cr; a.z = -sr;
		c.x =  sr; c.z =  cr;
	}
	if (axis=='Z') {   c.z = 1.0;
		a.x =  cr; a.y = -sr;
		b.x =  sr; b.y =  cr;
	}

	R->A.x = a.x*A.x + b.x*A.y + c.x*A.z;
        R->B.x = a.x*B.x + b.x*B.y + c.x*B.z;
        R->C.x = a.x*C.x + b.x*C.y + c.x*C.z;

	R->A.y = a.y*A.x + b.y*A.y + c.y*A.z;
        R->B.y = a.y*B.x + b.y*B.y + c.y*B.z;
        R->C.y = a.y*C.x + b.y*C.y + c.y*C.z;

	R->A.z = a.z*A.x + b.z*A.y + c.z*A.z;
        R->B.z = a.z*B.x + b.z*B.y + c.z*B.z;
        R->C.z = a.z*C.x + b.z*C.y + c.z*C.z;
}

void Mframe ( Vec a, Vec b, Vec c, Mat *M ) {
// sets an orthogonal coordinate frame on b with y as the bisector of abc and x close to ac
Vec	x,y,z;
	x.x=c.x-a.x; x.y=c.y-a.y; x.z=c.z-a.z;
	y.x=b.x-0.5*(c.x+a.x); y.y=b.y-0.5*(c.y+a.y); y.z=b.z-0.5*(c.z+a.z);
	vprod(x,y,&z); vprod(y,z,&x);
	vnorm(&x); vnorm(&y); vnorm(&z);
	VtoM(x,y,z,M);
}

int line2tri ( Vec a, Vec b, Vec c, Vec d, Vec e ) {
/* TRUE if line segment d-e cuts triangle a-b-c */
Vec x, y, z;
Mat M[1], W[1];
float det;
	vsub(b,a,&x);
	vsub(c,a,&y);
	vsub(d,e,&z);
	VtoM(x,y,z,M);
	det = Mdet(M);
	if (fabs(det) < 0.00001) return 0;
	Minv(M,W,det);
	vsub(d,a,&d);
	MmulV(W,d,&e);
	if (e.x < 0.0) return 0;
	if (e.y < 0.0) return 0;
	if (e.z < 0.0) return 0;
	if (e.z > 1.0) return 0;
	if (e.x + e.y > 1.0) return 0;
	return 1;
}

int line2line ( Vec a, Vec b, Vec c, Vec d, float s ) {
/* true if line segments a-b and d-e are closer than s */
Vec x, y, z;
Mat M[1], W[1];
float det;
Vec e;
        vsub(b,a,&x);
        vsub(d,c,&y);
        vprod(x,y,&z);
        vnorm(&z);
        VtoM(x,y,z,M);
        det = Mdet(M);
        Minv(M,W,det);
        vsub(d,a,&d);
        MmulV(W,d,&e);
        if (e.x < 0.0) return 0;
        if (e.y < 0.0) return 0;
        if (e.x > 1.0) return 0;
        if (e.y > 1.0) return 0;
        if (e.z >  s ) return 0;
        if (e.z < -s ) return 0;
        return 1;
}

float dist2line ( Vec a, Vec b, Vec c, Vec d ) {
/* min dist between line segments a-b and d-e */
Vec x, y, z;
Mat M[1], W[1];
float det, dmin;
Vec e,f;
        vsub(b,a,&x);
        vsub(d,c,&y);
        vprod(x,y,&z);
        vnorm(&z);
        VtoM(x,y,z,M);
        det = Mdet(M);
        Minv(M,W,det);
        vsub(d,a,&f);
        MmulV(W,f,&e);
        if ((e.x>0.0 && e.x<1.0) && (e.y>0.0 && e.y<1.0)) { /* in both lines */
		return fabs(e.z);
	}
        if (e.x>0.0 && e.x<1.0) { /* just in a-b */
		return fmin(endOline(a,b,c),endOline(a,b,d));
	}
        if (e.y>0.0 && e.y<1.0) { /* just in c-d */
		return fmin(endOline(c,d,a),endOline(c,d,b));
	}
	return fmin(fmin(vdif(a,c),vdif(a,d)),fmin(vdif(b,c),vdif(b,d)));
}

float lineOline ( Vec a, Vec b, Vec c, Vec d, Vec *box ) {
/* return overlap length for line segments a-b and c-d */
Vec x, y, z;
Mat M[1], W[1];
float lap, det, aa, ga, gb, hc, hd, r, s, t[4];
float ab,cd,ac,bd;
Vec e,f,g,h, pox[4],aox[4];
int i, key[4], ley[4];
int parr=0, perp=0;
	vcopy(a,aox+0); vcopy(b,aox+1); vcopy(c,aox+2); vcopy(d,aox+3);
	vcopy(a,box+0); vcopy(b,box+1); vcopy(c,box+2); vcopy(d,box+3);
	vsub(b,a,&x); vnorm(&x);
	vsub(d,c,&y); vnorm(&y);
	vprod(x,y,&z); vnorm(&z);
	if (fabs(vdot(x,y)) > 0.9999) { int in = 0;
        // parallel or anti
                if (in==0) in = dot2pair(a,b,c,&g,&h);
                if (in==0) in = dot2pair(a,b,d,&g,&h);
                if (in==0) in = dot2pair(c,d,a,&h,&g);
                if (in==0) in = dot2pair(c,d,b,&h,&g);
                if (in==0) return 0.0;
		parr=1;
        } else {
		VtoM(x,y,z,M);
		det = Mdet(M);
		Minv(M,W,det);
		vsub(d,a,&e);
		MmulV(W,e,&f);
		vmul(&x,f.x);
		vadd(a,x,&g); /* g = top of mut.perp. line */
		vmul(&y,f.y);
		vsub(d,y,&h); /* h = bot of mut.perp. line */
	}
	ga = vdif(g,a); gb = vdif(g,b);
	hc = vdif(h,c); hd = vdif(h,d);
	vsub(a,g,&a); vsub(b,g,&b);
	vsub(c,h,&c); vsub(d,h,&d);
	if (vdot(a,c)==0.0) {
        // perpendicular
                if (fabs(ga-hc) > fabs(ga-hd)) hc = -hc;
		perp=1;
        } else {
                if (vdot(a,c)<0.0) hc = -hc;
        }
	if (vdot(a,b)<0.0) gb = -gb;
	if (vdot(a,d)<0.0) hd = -hd;
	// dists now relative to the contact normal gh: right+, left-
	t[0] = ga; t[1] = gb; t[2] = hc; t[3] = hd; 
	sort(0,t,0,key,4,1);
	for (i=0; i<4; i++) if (key[i]<2) ley[i] = 0; else ley[i] = 1;
	if (ley[0]==ley[1]) return 0.0;
	lap = fabs(t[key[1]]-t[key[2]]);
	if (box)
	{ int ke1 = key[1], ke2 = key[2],
	      le1 = ley[1], le2 = ley[2];
	  Vec tmp1, tmp2, tmp;
		if (ga > gb) vsub(box[0],box[1],&x);
			else vsub(box[1],box[0],&x);
		if (hc > hd) vsub(box[2],box[3],&y);
			else vsub(box[3],box[2],&y);
		vnorm(&x); vnorm(&y);
		vcopy(box[ke1],&tmp1);
		vcopy(box[ke2],&tmp2);
		vcopy(tmp1,box+0);
		if (le1==le2) { /* contained */
			vcopy(tmp2,box+1);
			if (le2==0) {
				 vcopy(y,&tmp); vmul(&tmp,t[ke2]); vadd(h,tmp,box+3);
			} else { vcopy(x,&tmp); vmul(&tmp,t[ke2]); vadd(g,tmp,box+3); }
		} else { /* staggered */
			vcopy(tmp2,box+3);
			if (le2==0) {
				 vcopy(y,&tmp); vmul(&tmp,t[ke2]); vadd(h,tmp,box+1);
			} else { vcopy(x,&tmp); vmul(&tmp,t[ke2]); vadd(g,tmp,box+1); }
		}
		if (le1==0) {
			 vcopy(y,&tmp); vmul(&tmp,t[ke1]); vadd(h,tmp,box+2);
		} else { vcopy(x,&tmp); vmul(&tmp,t[ke1]); vadd(g,tmp,box+2); }
	}
	/* restore ab cd line order in box (but box lines run parallel */
	if (ley[1]) {
		pox[0] = box[2]; pox[1] = box[3];
		pox[2] = box[0]; pox[3] = box[1];
		for (i=0; i<4; i++) box[i] = pox[i];
	}
        ab = vdif(box[0],box[1]);
        cd = vdif(box[2],box[3]);
        if (fabs(ab-cd) > 0.0001) { float gh;
                ac = vdif(box[0],box[2]);
                bd = vdif(box[1],box[3]);
                printf("*NB* box has diff edges: ab=%f, cd=%f, ac=%f, bd=%f\n", ab,cd,ac,bd);
                Pv(aox[0]) NL Pv(aox[1]) NL Pv(aox[2]) NL Pv(aox[3]) NL Pv(h) NL Pv(g) NL
                ab = vdif(aox[0],aox[1]); cd = vdif(aox[2],aox[3]);
                ac = vdif(aox[0],aox[2]); bd = vdif(aox[1],aox[3]);
		gh = vdif(g,h);
                printf("     end-end distances : ab=%f, cd=%f, ac=%f, bd=%f, gh=%f\n", ab,cd,ac,bd,gh);
                printf("     end-box distances : ga=%f, gb=%f, hc=%f, hd=%f\n", ga,gb,hc,hd);
                ga = vdif(aox[0],g); gb = vdif(aox[1],g);
                hc = vdif(aox[0],h); hd = vdif(aox[3],h);
                printf("     end-box distances : ga=%f, gb=%f, hc=%f, hd=%f\n", ga,gb,hc,hd);
		if(perp) printf("perpendicular\n");
		if(parr) printf("parallel\n");
        }
	return lap;
}

int dot2pair ( Vec a, Vec b, Vec c, Vec *g, Vec *h ) {
/* returns 1 if c lies over the line a-b and image of c on a-b in g plus copy of c in h */
Vec p, q;
float d;
	vsub(b,a,&p);
	vsub(c,a,&q);
	d = vdot(p,q)/vsqr(p);
	vmul(&p,d);
	vadd(a,p,g);
	vcopy(c,h);
	if ( d<0.0 || d>1.0 ) return 0;
	return 1;
}

float old_lineOline ( Vec a, Vec b, Vec c, Vec d, Vec *box ) {
/* return overlap length for line segments a-b and c-d */
Vec x, y, z;
Mat M[1], W[1];
float lap, det, aa, ga, gb, hc, hd, r, s, t[4];
Vec e,f,g,h, pox[4];
int i, key[4], ley[4];
	vcopy(a,box+0); vcopy(b,box+1); vcopy(c,box+2); vcopy(d,box+3);
	vsub(b,a,&x); vnorm(&x);
	vsub(d,c,&y); vnorm(&y);
	vprod(x,y,&z); vnorm(&z);
	if (vdot(x,y) < 0.0001) { 
		vcopy(x,&e); vmul(&e,0.0001);
		vsub(c,e,&c); vadd(d,e,&d);
		vsub(d,c,&y); vnorm(&y);
	 }
	VtoM(x,y,z,M);
	det = Mdet(M);
	Minv(M,W,det);
	vsub(d,a,&e);
	MmulV(W,e,&f);
	vmul(&x,f.x);
	vadd(a,x,&g); /* g = top of mut.perp. line */
	vmul(&y,f.y);
	vsub(d,y,&h); /* h = bot of mut.perp. line */
	ga = vdif(g,a); gb = vdif(g,b);
	hc = vdif(h,c); hd = vdif(h,d);
	vsub(a,g,&a); vsub(b,g,&b);
	vsub(c,h,&c); vsub(d,h,&d);
	aa = vsqr(a);
	r = vdot(a,c)/aa;
	s = vdot(a,d)/aa;
	if (vdot(a,b)<0.0) gb = -gb;
	if (r<0.0) hc = -hc;
	if (s<0.0) hd = -hd;
	t[0] = ga; t[1] = gb; t[2] = hc; t[3] = hd; 
	sort(0,t,0,key,4,1);
	for (i=0; i<4; i++) if (key[i]<2) ley[i] = 0; else ley[i] = 1;
	if (ley[0]==ley[1]) return 0.0;
	lap = fabs(t[key[1]]-t[key[2]]);
	if (box)
	{ int ke1 = key[1], ke2 = key[2],
	      le1 = ley[1], le2 = ley[2];
	  Vec tmp1, tmp2, tmp;
		if (ga > gb) vsub(box[0],box[1],&x);
			else vsub(box[1],box[0],&x);
		if (hc > hd) vsub(box[2],box[3],&y);
			else vsub(box[3],box[2],&y);
		vnorm(&x); vnorm(&y);
		vcopy(box[ke1],&tmp1);
		vcopy(box[ke2],&tmp2);
		vcopy(tmp1,box+0);
		if (le1==le2) { /* contained */
			vcopy(tmp2,box+1);
			if (le2==0) {
				 vcopy(y,&tmp); vmul(&tmp,t[ke2]); vadd(h,tmp,box+3);
			} else { vcopy(x,&tmp); vmul(&tmp,t[ke2]); vadd(g,tmp,box+3); }
		} else { /* staggered */
			vcopy(tmp2,box+3);
			if (le2==0) {
				 vcopy(y,&tmp); vmul(&tmp,t[ke2]); vadd(h,tmp,box+1);
			} else { vcopy(x,&tmp); vmul(&tmp,t[ke2]); vadd(g,tmp,box+1); }
		}
		if (le1==0) {
			 vcopy(y,&tmp); vmul(&tmp,t[ke1]); vadd(h,tmp,box+2);
		} else { vcopy(x,&tmp); vmul(&tmp,t[ke1]); vadd(g,tmp,box+2); }
	}
	/* restore ab cd line order in box (but box lines run parallel */
	if (ley[1]) {
		pox[0] = box[2]; pox[1] = box[3];
		pox[2] = box[0]; pox[3] = box[1];
		for (i=0; i<4; i++) box[i] = pox[i];
	}
	return lap;
}

int line2dot ( Vec a, Vec b, Vec c, float s ) {
/* TRUE if c lies over line segment a-b closer than s */
Vec p, q;
float d;
	vsub(b,a,&p);
	vsub(c,a,&q);
	d = vdot(p,q)/vsqr(p);
	if (d < 0.0 || d > 1.0) return 0;
	vmul(&p,d);
	vsub(q,p,&q);
	d = vsqr(q);
	if (d>s*s) return 0;
	return 1;
}

float endOline ( Vec a, Vec b, Vec c ) {
/* returns distance of c to an extended line a-b  */
Vec p, q;
float d;
	vsub(b,a,&p);
	vsub(c,a,&q);
	d = vdot(p,q)/vsqr(p);
	vmul(&p,d);
	vsub(q,p,&q);
	return vmod(q);
}

float dotOline ( Vec a, Vec b, Vec c, Vec *e ) {
/* returns distance of c to an extended line a-b and image of c on a-b in e */
Vec p, q;
float d;
	vsub(b,a,&p);
	vsub(c,a,&q);
	d = vdot(p,q)/vsqr(p);
	vmul(&p,d);
	vsub(q,p,&q);
	d = vsqr(q);
	vadd(a,p,&p);
	if (e) { e->x = p.x; e->y = p.y; e->z = p.z; }
	return sqrt(d);
}

void rotate ( Vec a, Vec b, Vec *c, float t ) {
/* rotate c about a-b by t (+t = clockwise viewed along a-->b)*/
Vec x,y,z, d;
Mat M[1],R[1];
	identM(R);
	R->B.y = cos(t); R->B.z = -sin(t);
	R->C.y = sin(t); R->C.z =  cos(t);
	vsub( b,a,&x);
	if (vsqr(x) < NOISE) return; // b lies on a
	vsub(*c,a,&y);
	if (vsqr(y) < NOISE) return; // c lies on a
	vprod(x,y,&z);
	if (vsqr(z) < NOISE) return; // c lies on a--b
	vprod(x,z,&y);
	vnorm(&x); vnorm(&y); vnorm(&z);
	VtoM(x,y,z,M);
	vsub(*c,a,&d); // a -> zero
	VmulM(M,d,&d); // a-b to X
	VmulM(R,d,&d); // rotate X
	MmulV(M,d,&d); // X to a-b
	vadd (a,d,c);  // back to a  
}

float angle ( Vec a, Vec b, Vec c ) {
/* returns the (unsigned) angle at b */
Vec x, y;
	vsub(a,b,&x); vnorm(&x);
	vsub(c,b,&y); vnorm(&y);
	return acos(vdot(x,y));
}

float torsion ( Vec a, Vec b, Vec c, Vec d ) {
/* returns the torsion angle down b-c */
Vec x, y, z, p, q, r;
float vol, cos, sin, tor;
	vsub(a,b,&x);
	vsub(b,c,&y);
	vsub(c,d,&z);
	vprod(x,y,&p);
	vnorm(&p);
	vprod(y,z,&q);
	vnorm(&q);
	cos = vdot(p,q);
	vprod(p,q,&r);
	vol = vtri(p,q,y);
	sin = vmod(r);
	if (vol < 0.0) sin = -sin;
	tor = angle1pi(sin,cos);
	return tor;
}

float	angle1pi ( float s, float c )
{
	if (s>=0.0 && c>=0.0) {
		if (s<0.5) return asin(s);
		      else return acos(c);
	}
	if (s>=0.0 && c<=0.0) {
		if (s<0.5) return PI - asin(s);
		      else return PI - acos(-c);
	}
	if (s<=0.0 && c<=0.0) {
		if (s>-.5) return asin(-s) - PI;
		      else return acos(-c) - PI;
	}
	if (s<=0.0 && c>=0.0) {
		if (s>-.5) return -asin(-s);
		      else return -acos(c);
	}
	printf("angle1pi(s,c) out of range: s = %f, c = %f\n", s,c);
	exit(1);
}

float	angle2pi ( float s, float c )
{
	if (s>=0.0 && c>=0.0) {
		if (s<0.5) return asin(s);
		      else return acos(c);
	}
	if (s>=0.0 && c<=0.0) {
		if (s<0.5) return PI - asin(s);
		      else return PI - acos(-c);
	}
	if (s<=0.0 && c<=0.0) {
		if (s>-.5) return PI + asin(-s);
		      else return PI + acos(-c);
	}
	if (s<=0.0 && c>=0.0) {
		if (s>-.5) return twoPI - asin(-s);
		      else return twoPI - acos(c);
	}
	printf("angle2pi(s,c) out of range: s = %f, c = %f\n", s,c);
	exit(1);
}

float angdif ( float a, float b ) {
// absolute difference of two angles (in +/- PI range)
        if (a<0.0 && b>0.0)
        { float d = b-a;
                if (d<PI) return d;
                     else return twoPI-d;
        }
        if (a>0.0 && b<0.0)
        { float d = a-b;
                if (d<PI) return d;
                     else return twoPI-d;
        }
        if (a>b) return a-b; else return b-a;
}

void separate ( Vec *a, Vec *b, float dist, float kick )
{
Vec     oldm, oldn, disp, mid;
float   gap, dif;
	if (kick < 0.0) {
		kick = -kick;	// -ve kick bypasses upper limit
	} else {
		if (kick > 0.3) kick = 0.3; // over 0.3 oscillates
	}
        gap = vdif(*a,*b);
        dif = gap-dist;
	if (dif < 0.0) dif = -dif;
	if (dif < dist*0.00001) return;
	if (dif > 1.0) dif = 1.0;
        vcopy(*a,&oldm); vcopy(*b,&oldn);
        vave(oldm,oldn,&mid);
        vsub(oldm,mid,&disp);
	vmul(&disp,kick*dif);
        if (dist<gap) {
                vsub(oldm,disp,a);
                vadd(oldn,disp,b);
        }
        if (dist>gap) {
                vadd(oldm,disp,a);
                vsub(oldn,disp,b);
        }
}
