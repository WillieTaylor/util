#include "code/util.def"
#include "code/geom.def"

main() {
Vec a, b, c, d;
Mat m;
	m.A.x = -0.6657;  m.A.y = -0.3078;  m.A.z = -0.6798;
  	m.B.x = -0.7459;  m.B.y =  0.2453;  m.B.z =  0.6193;
  	m.C.x = -0.0239;  m.C.y =  0.9193;  m.C.z = -0.3929;
	Mprint(&m);
	Pr(vmod(m.A));
	Pr(vmod(m.B));
	Pr(vmod(m.C));
	NL
	a.x = 1; a.y = 0; a.z = 0;
	b.x = 0; b.y = 1; b.z = 0;
	c.x = 0; c.y = 0; c.z = 1;
	Pv(a) NL Pv(b) NL Pv(c) NL
	NL
	MmulV(&m,a,&d); Pv(d) NL
	VmulM(&m,d,&d); Pv(d) NL
	NL
	MmulV(&m,b,&d); Pv(d) NL
	VmulM(&m,d,&d); Pv(d) NL
	NL
	MmulV(&m,c,&d); Pv(d) NL
	VmulM(&m,d,&d); Pv(d) NL
}
