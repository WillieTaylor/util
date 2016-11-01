#include "incl/util.h"
#include "incl/geom.h"
main() {
Vec	a, b, c, d, box[4];
float	over;
	vinit(&a); vinit(&b); vinit(&c); vinit(&d);
	b.x = 2.0; b.z = 2.0;
	c.x = -1.0; c.y = c.z = 0.9999 ;
	d.x = d.y = 1.0;
/*
 a.x = 44.498722; a.y = 17.620232; a.z =  9.540015;
 b.x = 32.492584; b.y = -2.832490; b.z = 10.433364;
 c.x = 34.181141; c.y = -3.054269; c.z =  3.592994;
 d.x = 41.066765; d.y = 12.441084; d.z = -0.621558;
*/
	Pv(a) NL Pv(b) NL Pv(c) NL Pv(d) NLL
	over = lineOline(a,b,c,d,box);
	Pr(over) NL
Pv(box[0]) NL Pv(box[1]) NL Pv(box[2]) NL Pv(box[3]) NL
	over = lineOline(b,a,c,d,box);
	Pr(over) NL
Pv(box[0]) NL Pv(box[1]) NL Pv(box[2]) NL Pv(box[3]) NL
	over = lineOline(a,b,d,c,box);
	Pr(over) NL
Pv(box[0]) NL Pv(box[1]) NL Pv(box[2]) NL Pv(box[3]) NL
	over = lineOline(b,a,d,c,box);
	Pr(over) NL
Pv(box[0]) NL Pv(box[1]) NL Pv(box[2]) NL Pv(box[3]) NL
}
