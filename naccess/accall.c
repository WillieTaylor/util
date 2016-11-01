/* accall.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__5 = 5;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__100 = 100;
static integer c__20000 = 20000;
static integer c__2000 = 2000;
static integer c__4 = 4;
static integer c__10 = 10;

/* accall_start */
/* Main program */ MAIN__()
{
    /* Initialized data */

    static real hyrad = (float)1.;
    static char rlab[3*3+1] = "RESHEMHOH";

    /* Format strings */
    static char fmt_102[] = "(\002 UNKNOWN residue type.............> \002,a\
10)";
    static char fmt_104[] = "(\002 NON-STANDARD atom.\002,a4,\002 in residue\
> \002,a10)";
    static char fmt_106[] = "(\002 ASSUMED vdw of \002,a4,\002 in \002,a10\
,\002 = \002,f5.2,\002 (same as \002,a3,\002)\002)";
    static char fmt_108[] = "(\002 GUESSED vdw of \002,a4,\002 in \002,a10\
,\002 = \002,f5.2)";
    static char fmt_110[] = "(\002 CHAINS   \002,i5,/,\002 RESIDUES \002,i5,\
/,\002 ATOMS    \002,i5)";
    static char fmt_120[] = "(\002REM  File of summed (Sum) and % (per.\
)\002,\002 accessibilities for \002,a)";
    static char fmt_125[] = "(\002REM RES _ NUM      All atoms   Non P sid\
e  \002,\002 Polar Side   Total Side   Main Chain\002)";
    static char fmt_130[] = "(\002REM                ABS   REL    ABS   RE\
L\002,\002    ABS   REL    ABS   REL    ABS   REL\002)";
    static char fmt_150[] = "(a3,1x,a10,1x,5(f7.2,f6.1))";
    static char fmt_160[] = "(\002END  Absolute sums over accessible surfa\
ce \002,/,\002TOTAL\002,8x,5(f8.1,5x))";

    /* System generated locals */
    address a__1[2];
    integer i__1[2], i__2;
    real r__1, r__2;
    char ch__1[260];
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_cmp();
    /* Subroutine */ int s_copy();
    integer i_indx(), s_wsfe(), do_fio(), e_wsfe();
    /* Subroutine */ int s_stop(), s_cat();
    integer f_open(), s_rsfi(), e_rsfi(), f_clos();

    /* Local variables */
    static real accs[20000];
    static char card[256];
    static integer flen;
    static logical falt;
    static real rads[20000];
    static char atom[4];
    static integer slen;
    static char last[10];
    static integer vlen, nats, resindex[20000], atomtype[20000];
    static char c__[256*256];
    static integer i__, j, k, l[256], n;
    static real bfact[20000];
    static char label[30*20000], fname[256], sname[256];
    static logical hetas, conta;
    static char vname[256];
    extern integer fopen_();
    static real occup[20000];
    extern integer parse_();
    static real probe;
    static logical dorsa;
    extern /* Subroutine */ int vanin_(), gatom_();
    extern doublereal readfloat_();
    static logical fullo;
    static integer atype;
    static logical resok;
    extern /* Subroutine */ int ratom_();
    static logical hydro;
    extern /* Subroutine */ int solva_();
    static logical start;
    extern integer what_atom__();
    extern /* Subroutine */ int tolow_();
    static real tsums[5];
    extern /* Subroutine */ int which3_();
    static char aacids[3*100], ch_nam__[1*100];
    static logical ok;
    static integer nacids;
    static char anames[4*100*100];
    static real vradii[10000]	/* was [100][100] */;
    static char resnam[10*2000], firsta[1];
    static real zslice;
    static integer num_chains__, numats[100];
    extern /* Subroutine */ int summer_(), vguess_();
    extern integer readstring_();
    static logical aok;
    static integer len;
    static char alt[1], res[3];
    static real vdw;
    static integer rty[2000], num_res__;
    static real xyz[60000]	/* was [20000][3] */;
    static logical wwaters;
    static real ressums[20000]	/* was [2000][5][2] */;

    /* Fortran I/O blocks */
    static cilist io___26 = { 0, 6, 0, "(2a)", 0 };
    static cilist io___32 = { 0, 4, 0, "(a)", 0 };
    static cilist io___33 = { 0, 4, 0, "(2a)", 0 };
    static cilist io___34 = { 0, 4, 0, "(a,f6.2)", 0 };
    static cilist io___35 = { 0, 4, 0, "(a,f6.3)", 0 };
    static cilist io___36 = { 0, 4, 0, "(2a)", 0 };
    static cilist io___37 = { 0, 4, 0, "(a)", 0 };
    static cilist io___38 = { 0, 4, 0, "(a)", 0 };
    static cilist io___39 = { 0, 4, 0, "(a)", 0 };
    static cilist io___40 = { 0, 4, 0, "(a)", 0 };
    static cilist io___41 = { 0, 4, 0, "(a)", 0 };
    static cilist io___42 = { 0, 4, 0, "(a)", 0 };
    static cilist io___43 = { 0, 4, 0, "(a,i3,a)", 0 };
    static cilist io___58 = { 0, 4, 0, fmt_102, 0 };
    static cilist io___63 = { 0, 4, 0, fmt_104, 0 };
    static cilist io___64 = { 0, 4, 0, fmt_106, 0 };
    static cilist io___65 = { 0, 4, 0, fmt_108, 0 };
    static icilist io___68 = { 0, card+54, 0, "(f)", 6, 1 };
    static icilist io___70 = { 0, card+60, 0, "(f)", 6, 1 };
    static icilist io___73 = { 0, card, 0, "(30x,3f8.3)", 256, 1 };
    static cilist io___75 = { 0, 4, 0, "(a)", 0 };
    static cilist io___76 = { 0, 4, 0, fmt_110, 0 };
    static cilist io___78 = { 0, 2, 0, "(a30,3f8.3,f6.2,f5.1,f7.3,1x,f5.2)", 
	    0 };
    static cilist io___79 = { 0, 2, 0, "(a30,3f8.3,f8.3,1x,f5.2)", 0 };
    static cilist io___80 = { 0, 4, 0, "(a)", 0 };
    static cilist io___83 = { 0, 4, 0, "(a)", 0 };
    static cilist io___84 = { 0, 3, 0, fmt_120, 0 };
    static cilist io___85 = { 0, 3, 0, fmt_125, 0 };
    static cilist io___86 = { 0, 3, 0, fmt_130, 0 };
    static cilist io___87 = { 0, 3, 0, fmt_150, 0 };
    static cilist io___88 = { 0, 3, 0, fmt_160, 0 };



/*     --  AIM:- */
/*     --  Input a Brookhaven entry file and output an */
/*     --  PDB format file after filtering/cleaning, including Van der */
/*     --  Waal radii, contained in an external file "vdw.radii". */

/*     --  INPUT:- */
/*     --  PDB format file, van der Waal radii file */

/*     --  OPTIONS:- */
/*     --  Inclusion of non-standard amino acids, het-atoms, waters */
/*     --  etc. Flagging of missing residues, chain-breaks, */
/*     --  non-standard atoms names, missing atoms. */

/*     --  AUTHOR: S. Hubbard 3/92. EMBL. */

/*     --  PARAMETERS: maxa = max atoms per residue */
/*     --              maxr = max number of "standard" residue types */
/*     --              maxs = max number of atoms in PDB file */


/*     -- functions */


/*     -- variables */

/*     data card /'accall.input'/ */

/*  f77 accall.f -o accall */

/*     -- defaults */

    hetas = FALSE_;
    hydro = FALSE_;
    wwaters = FALSE_;
    fullo = FALSE_;
    dorsa = TRUE_;
    conta = FALSE_;

/*     -- Get USER directives */

    while(readstring_(&c__5, card, &len, (ftnlen)256) >= 0) {
	n = parse_(card, &len, " ", c__, l, (ftnlen)256, (ftnlen)1, (ftnlen)
		256);
	tolow_(c__, l, (ftnlen)256);
	if (s_cmp(c__, "pdbf", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(fname, c__ + 256, (ftnlen)256, (ftnlen)256);
	    flen = l[1];
	} else if (s_cmp(c__, "vdwf", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(vname, c__ + 256, (ftnlen)256, (ftnlen)256);
	    vlen = l[1];
	} else if (s_cmp(c__, "stdf", (ftnlen)4, (ftnlen)4) == 0) {
	    s_copy(sname, c__ + 256, (ftnlen)256, (ftnlen)256);
	    slen = l[1];
	} else if (s_cmp(c__, "prob", (ftnlen)4, (ftnlen)4) == 0) {
	    probe = readfloat_(c__ + 256, &l[1], (ftnlen)256);
	} else if (s_cmp(c__, "zsli", (ftnlen)4, (ftnlen)4) == 0) {
	    zslice = readfloat_(c__ + 256, &l[1], (ftnlen)256);
	} else if (s_cmp(c__, "heta", (ftnlen)4, (ftnlen)4) == 0) {
	    hetas = TRUE_;
	} else if (s_cmp(c__, "hydr", (ftnlen)4, (ftnlen)4) == 0) {
	    hydro = TRUE_;
	} else if (s_cmp(c__, "wate", (ftnlen)4, (ftnlen)4) == 0) {
	    wwaters = TRUE_;
	} else if (s_cmp(c__, "full", (ftnlen)4, (ftnlen)4) == 0) {
	    fullo = TRUE_;
	} else if (s_cmp(c__, "asao", (ftnlen)4, (ftnlen)4) == 0) {
	    dorsa = FALSE_;
	} else if (s_cmp(c__, "cont", (ftnlen)4, (ftnlen)4) == 0) {
	    conta = TRUE_;
	}
    }

/*     --   open files */

    i__ = i_indx(fname, ".", (ftnlen)256, (ftnlen)1) - 1;
    k = 1;
    ok = FALSE_;
    for (j = i__; j >= 1; --j) {
	if (*(unsigned char *)&fname[j - 1] == '/' && ! ok) {
	    k = j + 1;
	    ok = TRUE_;
	}
    }
    if (fopen_(&c__1, fname, &flen, "old", (ftnlen)256, (ftnlen)3) == 0) {
	s_wsfe(&io___26);
	do_fio(&c__1, " PDB FILE INPUT ", (ftnlen)16);
	do_fio(&c__1, fname, flen);
	e_wsfe();
	s_stop("ERROR: unable to open PDB file", (ftnlen)30);
    }
    o__1.oerr = 0;
    o__1.ounit = 2;
    o__1.ofnmlen = i__ - (k - 1) + 4;
/* Writing concatenation */
    i__1[0] = i__ - (k - 1), a__1[0] = fname + (k - 1);
    i__1[1] = 4, a__1[1] = ".asa";
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)260);
    o__1.ofnm = ch__1;
    o__1.orl = 0;
    o__1.osta = "unknown";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
    if (dorsa) {
	o__1.oerr = 0;
	o__1.ounit = 3;
	o__1.ofnmlen = i__ - (k - 1) + 4;
/* Writing concatenation */
	i__1[0] = i__ - (k - 1), a__1[0] = fname + (k - 1);
	i__1[1] = 4, a__1[1] = ".rsa";
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)260);
	o__1.ofnm = ch__1;
	o__1.orl = 0;
	o__1.osta = "unknown";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
    }
    o__1.oerr = 0;
    o__1.ounit = 4;
    o__1.ofnmlen = i__ - (k - 1) + 4;
/* Writing concatenation */
    i__1[0] = i__ - (k - 1), a__1[0] = fname + (k - 1);
    i__1[1] = 4, a__1[1] = ".log";
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)260);
    o__1.ofnm = ch__1;
    o__1.orl = 0;
    o__1.osta = "unknown";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);

/*     -- Read in VDW radii for all residues/atoms */

    vanin_(vname, &vlen, &nacids, aacids, anames, numats, vradii, &c__100, &
	    c__100, (ftnlen)256, (ftnlen)3, (ftnlen)4);
    s_wsfe(&io___32);
    do_fio(&c__1, " ACCALL - Accessibility calculations", (ftnlen)36);
    e_wsfe();
    s_wsfe(&io___33);
    do_fio(&c__1, " PDB FILE INPUT ", (ftnlen)16);
    do_fio(&c__1, fname, flen);
    e_wsfe();
    s_wsfe(&io___34);
    do_fio(&c__1, " PROBE SIZE     ", (ftnlen)16);
    do_fio(&c__1, (char *)&probe, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___35);
    do_fio(&c__1, " Z-SLICE WIDTH  ", (ftnlen)16);
    do_fio(&c__1, (char *)&zslice, (ftnlen)sizeof(real));
    e_wsfe();
    s_wsfe(&io___36);
    do_fio(&c__1, " VDW RADII FILE ", (ftnlen)16);
    do_fio(&c__1, vname, vlen);
    e_wsfe();
    if (hetas) {
	s_wsfe(&io___37);
	do_fio(&c__1, " INCL HETATOMS", (ftnlen)14);
	e_wsfe();
    } else {
	s_wsfe(&io___38);
	do_fio(&c__1, " EXCL HETATOMS", (ftnlen)14);
	e_wsfe();
    }
    if (hydro) {
	s_wsfe(&io___39);
	do_fio(&c__1, " INCL HYDROGENS", (ftnlen)15);
	e_wsfe();
    } else {
	s_wsfe(&io___40);
	do_fio(&c__1, " EXCL HYDROGENS", (ftnlen)15);
	e_wsfe();
    }
    if (wwaters) {
	s_wsfe(&io___41);
	do_fio(&c__1, " INCL WATERS", (ftnlen)12);
	e_wsfe();
    } else {
	s_wsfe(&io___42);
	do_fio(&c__1, " EXCL WATERS", (ftnlen)12);
	e_wsfe();
    }
    s_wsfe(&io___43);
    do_fio(&c__1, " READVDW ", (ftnlen)9);
    do_fio(&c__1, (char *)&nacids, (ftnlen)sizeof(integer));
    do_fio(&c__1, " residues input", (ftnlen)15);
    e_wsfe();

/* --  Initialise variables/logicals */

    falt = FALSE_;
    start = TRUE_;
    s_copy(last, "          ", (ftnlen)10, (ftnlen)10);

/* -- Read data & decode */

    nats = 0;
    while(readstring_(&c__1, card, &len, (ftnlen)256) >= 0) {
	atype = 0;
	if (s_cmp(card, "ATOM", (ftnlen)4, (ftnlen)4) == 0) {
	    atype = 1;
	}
	if (s_cmp(card, "HETATM", (ftnlen)6, (ftnlen)6) == 0) {
	    atype = 2;
	}
	if (s_cmp(card + 17, "HOH", (ftnlen)3, (ftnlen)3) == 0) {
	    atype = 3;
	}
	if (atype == 1 || atype == 2 && hetas || atype == 3 && wwaters) {

/*     -- Ignore Alternate positions, other than blanks or 1st */
/*     -- encountered */

	    *(unsigned char *)alt = *(unsigned char *)&card[16];
	    if (*(unsigned char *)alt != ' ') {
		if (! falt) {
		    *(unsigned char *)firsta = *(unsigned char *)alt;
		    falt = TRUE_;
		}
		if (*(unsigned char *)alt != *(unsigned char *)firsta) {
		    goto L5;
		}
	    }

/*     -- Ignore hydrogens & deuteriums (unless flagged) */

	    if (*(unsigned char *)&card[13] == 'H' || *(unsigned char *)&card[
		    13] == 'D') {
		if (! hydro) {
		    goto L5;
		}
		vdw = hyrad;
		++nats;
		goto L6;
	    }

/*     -- Next atom */

	    ++nats;

/*     -- First residue ? */

	    if (start) {
		start = FALSE_;
		*(unsigned char *)&ch_nam__[0] = *(unsigned char *)&card[21];
		num_chains__ = 1;
	    }

/*     -- New residue ? */

	    if (s_cmp(last, card + 17, (ftnlen)10, (ftnlen)10) != 0) {
		s_copy(last, card + 17, (ftnlen)10, (ftnlen)10);
		++num_res__;
		s_copy(res, card + 17, (ftnlen)3, (ftnlen)3);
		s_copy(resnam + (num_res__ - 1) * 10, last, (ftnlen)10, (
			ftnlen)10);
		which3_(res, aacids, &i__, &nacids, &c__100, &resok, (ftnlen)
			3, (ftnlen)3);
		if (! resok) {
		    s_wsfe(&io___58);
		    do_fio(&c__1, card + 17, (ftnlen)10);
		    e_wsfe();
		}
		rty[num_res__ - 1] = atype;
	    }

/*     -- Get atom type */

	    s_copy(atom, card + 12, (ftnlen)4, (ftnlen)4);
	    atomtype[nats - 1] = what_atom__(atom, (ftnlen)4);

/*     -- Assign radius to atom */
/*     -- Special case(s) */

	    if (s_cmp(atom, "OXT", (ftnlen)4, (ftnlen)3) == 0) {
		vdw = (float)1.4;
		goto L6;
	    }

/*     -- known residue type */

	    if (resok) {
		vdw = (float)0.;
		ratom_(atom, anames, &i__, numats, &c__100, &c__100, &aok, &j,
			 (ftnlen)4, (ftnlen)4);
	    }

/*     -- Not OK, then try atoms in all residue types */

	    if (! resok || ! aok) {
		gatom_(atom, anames, &nacids, numats, &c__100, &c__100, &aok, 
			&i__, &j, (ftnlen)4, (ftnlen)4);
		s_wsfe(&io___63);
		do_fio(&c__1, atom, (ftnlen)4);
		do_fio(&c__1, card + 17, (ftnlen)10);
		e_wsfe();
		if (aok) {
		    s_wsfe(&io___64);
		    do_fio(&c__1, atom, (ftnlen)4);
		    do_fio(&c__1, card + 17, (ftnlen)10);
		    do_fio(&c__1, (char *)&vradii[i__ + j * 100 - 101], (
			    ftnlen)sizeof(real));
		    do_fio(&c__1, aacids + (i__ - 1) * 3, (ftnlen)3);
		    e_wsfe();
		}
	    }

/*     -- Still not OK, make a guess */

	    if (! aok) {
		vguess_(atom, &vdw, (ftnlen)4);
		s_wsfe(&io___65);
		do_fio(&c__1, atom, (ftnlen)4);
		do_fio(&c__1, card + 17, (ftnlen)10);
		do_fio(&c__1, (char *)&vdw, (ftnlen)sizeof(real));
		e_wsfe();
	    } else {
		vdw = vradii[i__ + j * 100 - 101];
	    }

/*     -- Store data */

L6:
	    rads[nats - 1] = vdw;
	    s_copy(label + (nats - 1) * 30, card, (ftnlen)30, (ftnlen)30);
	    if (fullo) {
		s_rsfi(&io___68);
		do_fio(&c__1, (char *)&occup[nats - 1], (ftnlen)sizeof(real));
		e_rsfi();
		s_rsfi(&io___70);
		do_fio(&c__1, (char *)&bfact[nats - 1], (ftnlen)sizeof(real));
		e_rsfi();
	    }
	    resindex[nats - 1] = num_res__;
	    s_rsfi(&io___73);
	    for (k = 1; k <= 3; ++k) {
		do_fio(&c__1, (char *)&xyz[nats + k * 20000 - 20001], (ftnlen)
			sizeof(real));
	    }
	    e_rsfi();
	}
L5:
	;
    }

/*     -- output */

    cl__1.cerr = 0;
    cl__1.cunit = 1;
    cl__1.csta = 0;
    f_clos(&cl__1);
    s_wsfe(&io___75);
    do_fio(&c__1, " ADDED VDW RADII", (ftnlen)16);
    e_wsfe();
    s_wsfe(&io___76);
    do_fio(&c__1, (char *)&num_chains__, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&num_res__, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nats, (ftnlen)sizeof(integer));
    e_wsfe();

/*     -- calculate atomic accessibilities */

    solva_(&nats, xyz, rads, accs, &probe, &zslice, &c__20000);
    if (conta) {
	i__2 = nats;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing 2nd power */
	    r__1 = rads[i__ - 1];
/* Computing 2nd power */
	    r__2 = rads[i__ - 1] + probe;
	    accs[i__ - 1] = accs[i__ - 1] * (r__1 * r__1) / (r__2 * r__2);
	}
    }
    i__2 = nats;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (fullo) {
	    s_wsfe(&io___78);
	    do_fio(&c__1, label + (i__ - 1) * 30, (ftnlen)30);
	    for (j = 1; j <= 3; ++j) {
		do_fio(&c__1, (char *)&xyz[i__ + j * 20000 - 20001], (ftnlen)
			sizeof(real));
	    }
	    do_fio(&c__1, (char *)&occup[i__ - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&bfact[i__ - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&accs[i__ - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&rads[i__ - 1], (ftnlen)sizeof(real));
	    e_wsfe();
	} else {
	    s_wsfe(&io___79);
	    do_fio(&c__1, label + (i__ - 1) * 30, (ftnlen)30);
	    for (j = 1; j <= 3; ++j) {
		do_fio(&c__1, (char *)&xyz[i__ + j * 20000 - 20001], (ftnlen)
			sizeof(real));
	    }
	    do_fio(&c__1, (char *)&accs[i__ - 1], (ftnlen)sizeof(real));
	    do_fio(&c__1, (char *)&rads[i__ - 1], (ftnlen)sizeof(real));
	    e_wsfe();
	}
    }
    s_wsfe(&io___80);
    do_fio(&c__1, " CALCULATED ATOMIC ACCESSIBILITES", (ftnlen)33);
    e_wsfe();
    if (dorsa) {
	summer_(sname, &slen, &nats, accs, &c__20000, &c__2000, atomtype, 
		resindex, resnam, ressums, tsums, (ftnlen)256, (ftnlen)10);
	s_wsfe(&io___83);
	do_fio(&c__1, " SUMMED ACCESSIBILITIES OVER RESIDUES", (ftnlen)37);
	e_wsfe();
	s_wsfe(&io___84);
	e_wsfe();
	s_wsfe(&io___85);
	e_wsfe();
	s_wsfe(&io___86);
	e_wsfe();
	i__2 = resindex[nats - 1];
	for (i__ = 1; i__ <= i__2; ++i__) {
	    s_wsfe(&io___87);
	    do_fio(&c__1, rlab + (rty[i__ - 1] - 1) * 3, (ftnlen)3);
	    do_fio(&c__1, resnam + (i__ - 1) * 10, (ftnlen)10);
	    for (j = 1; j <= 5; ++j) {
		do_fio(&c__1, (char *)&ressums[i__ + (j + 5) * 2000 - 12001], 
			(ftnlen)sizeof(real));
		do_fio(&c__1, (char *)&ressums[i__ + (j + 10) * 2000 - 12001],
			 (ftnlen)sizeof(real));
	    }
	    e_wsfe();
	}
	s_wsfe(&io___88);
	for (i__ = 1; i__ <= 5; ++i__) {
	    do_fio(&c__1, (char *)&tsums[i__ - 1], (ftnlen)sizeof(real));
	}
	e_wsfe();
    }
/* ---------------------------------------------------------- */
/* ---------------------------------------------------------- */
} /* MAIN__ */

/* Subroutine */ int vguess_(atom, vdw, atom_len)
char *atom;
real *vdw;
ftnlen atom_len;
{
    /* Builtin functions */
    integer s_cmp();

    *vdw = (float)1.8;

/* -- Make a guess then ! */

    if (*(unsigned char *)&atom[1] == 'C') {
	*vdw = (float)1.8;
    }
    if (*(unsigned char *)&atom[1] == 'N') {
	*vdw = (float)1.6;
    }
    if (*(unsigned char *)&atom[1] == 'S') {
	*vdw = (float)1.85;
    }
    if (*(unsigned char *)&atom[1] == 'O') {
	*vdw = (float)1.4;
    }
    if (*(unsigned char *)&atom[1] == 'P') {
	*vdw = (float)1.9;
    }
    if (s_cmp(atom, "CA", (ftnlen)2, (ftnlen)2) == 0) {
	*vdw = (float)2.07;
    }
    if (s_cmp(atom, "FE", (ftnlen)2, (ftnlen)2) == 0) {
	*vdw = (float)1.47;
    }
    if (s_cmp(atom, "CU", (ftnlen)2, (ftnlen)2) == 0) {
	*vdw = (float)1.78;
    }
    if (s_cmp(atom, "ZN", (ftnlen)2, (ftnlen)2) == 0) {
	*vdw = (float)1.39;
    }
    if (s_cmp(atom, "MG", (ftnlen)2, (ftnlen)2) == 0) {
	*vdw = (float)1.73;
    }
    return 0;
} /* vguess_ */

/* Subroutine */ int gatom_(atom, anames, nres, nats, maxr, maxa, ok, ir, ia, 
	atom_len, anames_len)
char *atom, *anames;
integer *nres, *nats, *maxr, *maxa;
logical *ok;
integer *ir, *ia;
ftnlen atom_len;
ftnlen anames_len;
{
    /* System generated locals */
    integer anames_dim1, anames_offset;

    /* Local variables */
    extern /* Subroutine */ int ratom_();

    /* Parameter adjustments */
    --nats;
    anames_dim1 = *maxr;
    anames_offset = 1 + anames_dim1 * 1;
    anames -= anames_offset * 4;

    /* Function Body */
    *ok = FALSE_;
    *ir = 0;
    while(*ir < *nres && ! (*ok)) {
	++(*ir);
	ratom_(atom, anames + (anames_offset << 2), ir, &nats[1], maxr, maxa, 
		ok, ia, (ftnlen)4, (ftnlen)4);
    }
    return 0;
} /* gatom_ */

/* Subroutine */ int ratom_(atom, anames, ires, nats, maxr, maxa, ok, find, 
	atom_len, anames_len)
char *atom, *anames;
integer *ires, *nats, *maxr, *maxa;
logical *ok;
integer *find;
ftnlen atom_len;
ftnlen anames_len;
{
    /* System generated locals */
    integer anames_dim1, anames_offset, i__1;

    /* Builtin functions */
    integer s_cmp();

    /* Local variables */
    static integer i__;


/* --  Checks to see if a "standard" atom name has been read in. */
/* --  Standard atom names are read from the file "vdw.radii", for */
/* --  the defined residue types therein. OK=.true. if found, and */
/* --  residue type = ires, and find = integer identifier of atom. */

    /* Parameter adjustments */
    --nats;
    anames_dim1 = *maxr;
    anames_offset = 1 + anames_dim1 * 1;
    anames -= anames_offset * 4;

    /* Function Body */
    *ok = FALSE_;
    *find = 0;
    if (*ires == 0 || *ires > *maxr) {
	return 0;
    }
    i__1 = nats[*ires];
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (s_cmp(atom, anames + (*ires + i__ * anames_dim1 << 2), (ftnlen)4, 
		(ftnlen)4) == 0) {
	    *find = i__;
	    *ok = TRUE_;
	    return 0;
	}
    }
    return 0;
} /* ratom_ */

/* Subroutine */ int which3_(res, acids, ires, nacids, maxr, ok, res_len, 
	acids_len)
char *res, *acids;
integer *ires, *nacids, *maxr;
logical *ok;
ftnlen res_len;
ftnlen acids_len;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_cmp();

    /* Local variables */
    static integer i__;


/* -- Search array "acids" for existence of residue "res". */
/* -- OK = .true. if found, and index of res returned in "ires". */

    /* Parameter adjustments */
    acids -= 3;

    /* Function Body */
    *ires = *nacids + 1;
    i__1 = *nacids;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (s_cmp(res, acids + i__ * 3, (ftnlen)3, (ftnlen)3) == 0) {
	    *ires = i__;
	    *ok = TRUE_;
	    return 0;
	}
    }
    *ok = FALSE_;
    return 0;
} /* which3_ */

integer fopen_(iochan, filename, len, fstat, filename_len, fstat_len)
integer *iochan;
char *filename;
integer *len;
char *fstat;
ftnlen filename_len;
ftnlen fstat_len;
{
    /* System generated locals */
    integer ret_val, i__1;
    olist o__1;

    /* Builtin functions */
    integer f_open();

    o__1.oerr = 1;
    o__1.ounit = *iochan;
    o__1.ofnmlen = *len;
    o__1.ofnm = filename;
    o__1.orl = 0;
    o__1.osta = fstat;
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L100;
    }
    ret_val = 1;
    return ret_val;
L100:
    ret_val = 0;
    return ret_val;
} /* fopen_ */

integer readstring_(file, card, len, card_len)
integer *file;
char *card;
integer *len;
ftnlen card_len;
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop();
    integer s_rsfe(), do_fio(), e_rsfe();
    /* Subroutine */ int s_copy();

    /* Fortran I/O blocks */
    static cilist io___91 = { 1, 0, 1, "(q,a)", 0 };


    if (*file > 200) {
	s_stop("ERROR: file number too large", (ftnlen)28);
    }
    io___91.ciunit = *file;
    i__1 = s_rsfe(&io___91);
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = do_fio(&c__1, (char *)&(*len), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = do_fio(&c__1, card, (*len));
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = *len;
    s_copy(card + i__1, "\000", *len + 1 - i__1, (ftnlen)1);
    ret_val = *len;
    return ret_val;
L100:
    ret_val = -1;
    return ret_val;
} /* readstring_ */

integer parse_(card, length, separator, chars, clen, card_len, separator_len, 
	chars_len)
char *card;
integer *length;
char *separator, *chars;
integer *clen;
ftnlen card_len;
ftnlen separator_len;
ftnlen chars_len;
{
    /* System generated locals */
    integer ret_val, i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy();

    /* Local variables */
    static integer i__;
    static logical search;
    static integer charpos[512]	/* was [256][2] */;

    /* Parameter adjustments */
    --clen;
    chars -= chars_len;

    /* Function Body */
    ret_val = 0;
    i__ = 0;
    search = FALSE_;
    while(i__ < *length) {
	++i__;
	if (! search) {
	    if (*(unsigned char *)&card[i__ - 1] != *(unsigned char *)
		    separator) {
		++ret_val;
		charpos[ret_val - 1] = i__;
		search = TRUE_;
	    }
	} else {
	    if (*(unsigned char *)&card[i__ - 1] == *(unsigned char *)
		    separator) {
		charpos[ret_val + 255] = i__ - 1;
		search = FALSE_;
	    }
	}
    }
    if (search) {
	charpos[ret_val + 255] = *length;
    }
    i__1 = ret_val;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = charpos[i__ - 1] - 1;
	s_copy(chars + i__ * chars_len, card + i__2, chars_len, charpos[i__ + 
		255] - i__2);
	clen[i__] = charpos[i__ + 255] - charpos[i__ - 1] + 1;
    }
    return ret_val;
} /* parse_ */

doublereal readfloat_(card, len, card_len)
char *card;
integer *len;
ftnlen card_len;
{
    /* System generated locals */
    integer i__1;
    real ret_val;
    icilist ici__1;

    /* Builtin functions */
    integer s_rsli(), do_lio(), e_rsli();

    /* Local variables */
    static real value;

    ici__1.icierr = 1;
    ici__1.iciend = 0;
    ici__1.icirnum = 1;
    ici__1.icirlen = *len;
    ici__1.iciunit = card;
    ici__1.icifmt = 0;
    i__1 = s_rsli(&ici__1);
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = do_lio(&c__4, &c__1, (char *)&value, (ftnlen)sizeof(real));
    if (i__1 != 0) {
	goto L100;
    }
    i__1 = e_rsli();
    if (i__1 != 0) {
	goto L100;
    }
    ret_val = value;
    return ret_val;
L100:
    ret_val = (float)-999.9;
    return ret_val;
} /* readfloat_ */

/* Subroutine */ int tolow_(text, len, text_len)
char *text;
integer *len;
ftnlen text_len;
{
    /* Initialized data */

    static char alph[1*26*2+1] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnop\
qrstuvwxyz";

    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j;
    static logical ok;

    i__1 = *len;
    for (i__ = 1; i__ <= i__1; ++i__) {
	j = 0;
	ok = TRUE_;
	while(j < 26 && ok) {
	    ++j;
	    if (*(unsigned char *)&alph[j - 1] == *(unsigned char *)&text[i__ 
		    - 1]) {
		ok = FALSE_;
		*(unsigned char *)&text[i__ - 1] = *(unsigned char *)&alph[j 
			+ 25];
	    }
	}
    }
    return 0;
} /* tolow_ */

/* Subroutine */ int solva_(nats, xyz, rads, accs, probe, zslice, maxs)
integer *nats;
real *xyz, *rads, *accs, *probe, *zslice;
integer *maxs;
{
    /* Initialized data */

    static real xmin = (float)9999.;
    static real ymin = (float)9999.;
    static real zmin = (float)9999.;
    static real xmax = (float)-9999.;
    static real ymax = (float)-9999.;
    static real zmax = (float)-9999.;

    /* System generated locals */
    integer xyz_dim1, xyz_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Builtin functions */
    double acos();
    /* Subroutine */ int s_stop();
    double sqrt(), atan2();
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static real area, arcf[2000], beta;
    static integer cube[20000];
    static real arci[2000];
    static integer itab[10000], karc, idim, mkji, natm[1500000]	/* was [150][
	    10000] */;
    static real rmax;
    static integer inov[2000];
    static real zres, rrsq, b, d__[2000];
    static integer i__, j, k, l, m, n;
    static real t, alpha, parea;
    static integer jidim;
    static real radsq[20000], rsecn, rsecr, zgrid, trig_test__, rsec2n, 
	    rsec2r;
    static integer in, io;
    static real pi, tf;
    static integer ir;
    static real dx[2000], dy[2000];
    static integer nm;
    static real ti, rr;
    static integer kjidim;
    static real tt, xr, yr, zr, arcsum;
    extern /* Subroutine */ int sortag_();
    static real rad[20000];
    static integer tag[2000], kji, ict;
    static real dsq[2000];
    static integer nzp;
    static real pix2, rrx2;

    /* Fortran I/O blocks */
    static cilist io___163 = { 0, 4, 0, "(a)", 0 };


/* ** */
/* ************************************************************* */
/* **  SOLVA - LEE & RICHARDS TYPE ACCESSIBLITY CALCULATIONS ** */
/* ************************************************************* */
/* ** */
/* ** Calculate accessible surface area for a group of atoms. */
/* ** The accessible area for a given atom is calculated by the */
/* ** formula: */
/* **     (arcsum) x (atom radius+probe radius) x (deltaz) */
/* ** Numerical integration is carried out over z. in each z- */
/* ** section, the arcsum for a given atom is the arclength of */
/* ** the circle (intersection of the atom sphere with the z- */
/* ** section) that is not interior to any other atom circles */
/* ** in the same z-section. */
/* ** */
/* ************************************************************* */
/* ** */
/* **  error parameter  - this gives accuracy of calculation */
/* **                   - suitable values are 0.01 (high */
/* **                     accuracy) to 0.1 (low accuracy) */
/* **                   - in detail the z sections are spaced */
/* **                     at about error*diameter of atom */
/* ** */
/* **  probe size       - radius of probe in angstroms */
/* **                   - suitable value for water = 1.4 */
/* ** */
/* ************************************************************* */
/* ============================================================= */

/* PARAMETERS: */
/* ncube = maximum number of cubes allowed for placing of atoms */
/* nac   = maximum number of atoms per cube */
/* nint  = maximum number of sphere intersections */
/* smax  = maximum number of atoms */


/*     the following are dimensioned to the max no of atoms (maxs) */


/*     the following are dimensioned to the max no of intersections */
/*     of neighbouring spheres (nint) */

    /* Parameter adjustments */
    --accs;
    --rads;
    xyz_dim1 = *maxs;
    xyz_offset = 1 + xyz_dim1 * 1;
    xyz -= xyz_offset;

    /* Function Body */

/*     initialise variables, constants */

    ict = 2000;
    pi = acos((float)-1.);
    pix2 = pi * (float)2.;

/*     -- Radius of an atom sphere = atom radius + probe radius */
/*     -- Find maxima and minima */

    rmax = (float)0.;
    karc = ict;
    i__1 = *nats;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rad[i__ - 1] = rads[i__] + *probe;
/* Computing 2nd power */
	r__1 = rad[i__ - 1];
	radsq[i__ - 1] = r__1 * r__1;
	if (rad[i__ - 1] > rmax) {
	    rmax = rad[i__ - 1];
	}
	if (xmin > xyz[i__ + xyz_dim1]) {
	    xmin = xyz[i__ + xyz_dim1];
	}
	if (ymin > xyz[i__ + (xyz_dim1 << 1)]) {
	    ymin = xyz[i__ + (xyz_dim1 << 1)];
	}
	if (zmin > xyz[i__ + xyz_dim1 * 3]) {
	    zmin = xyz[i__ + xyz_dim1 * 3];
	}
	if (xmax < xyz[i__ + xyz_dim1]) {
	    xmax = xyz[i__ + xyz_dim1];
	}
	if (ymax < xyz[i__ + (xyz_dim1 << 1)]) {
	    ymax = xyz[i__ + (xyz_dim1 << 1)];
	}
	if (zmax < xyz[i__ + xyz_dim1 * 3]) {
	    zmax = xyz[i__ + xyz_dim1 * 3];
	}
    }

/*     rmax = max diameter */

    rmax *= (float)2.;

/*     -- Cubicals containing the atoms are setup. */
/*     -- The dimension of an edge equals the largest atom sphere radius */
/*     -- The cubes have a single index */
/*     -- Minimum of 3 by 3 cubic grid */
/*     -- EXIT if max cubes exceeded */

    idim = (xmax - xmin) / rmax + (float)1.;
    if (idim < 3) {
	idim = 3;
    }
    jidim = (ymax - ymin) / rmax + (float)1.;
    if (jidim < 3) {
	jidim = 3;
    }
    jidim = idim * jidim;
    kjidim = (zmax - zmin) / rmax + (float)1.;
    if (kjidim < 3) {
	kjidim = 3;
    }
    kjidim = jidim * kjidim;
    if (kjidim > 10000) {
	s_stop("SOLVA_ERROR: max cubes exceeded", (ftnlen)31);
    }

/*     -- Prepare upto ncube cubes each containing upto nac atoms. The cube index */
/*     -- is kji. The atom index for each cube is in itab */

    for (l = 1; l <= 10000; ++l) {
	itab[l - 1] = 0;
    }
    i__1 = *nats;
    for (l = 1; l <= i__1; ++l) {
	i__ = (xyz[l + xyz_dim1] - xmin) / rmax + (float)1.;
	j = (xyz[l + (xyz_dim1 << 1)] - ymin) / rmax;
	k = (xyz[l + xyz_dim1 * 3] - zmin) / rmax;
	kji = k * jidim + j * idim + i__;
	n = itab[kji - 1] + 1;
	if (n > 150) {
	    s_stop("SOLVA_ERROR: max atoms per cube exceeded", (ftnlen)40);
	}
	itab[kji - 1] = n;
	natm[n + kji * 150 - 151] = l;
	cube[l - 1] = kji;
    }

/*     -- Process each atom in turn */

    nzp = (float)1. / *zslice + (float).5;
    i__1 = *nats;
    for (ir = 1; ir <= i__1; ++ir) {
	kji = cube[ir - 1];
	io = 0;
	area = (float)0.;
	xr = xyz[ir + xyz_dim1];
	yr = xyz[ir + (xyz_dim1 << 1)];
	zr = xyz[ir + xyz_dim1 * 3];
	rr = rad[ir - 1];
	rrx2 = rr * (float)2.;
	rrsq = radsq[ir - 1];

/*     -- Find the 'mkji' cubes neighboring the kji cube */

	for (k = -1; k <= 1; ++k) {
	    for (j = -1; j <= 1; ++j) {
		for (i__ = -1; i__ <= 1; ++i__) {
		    mkji = kji + k * jidim + j * idim + i__;
		    if (mkji >= 1) {
			if (mkji > kjidim) {
			    goto L14;
			}
			nm = itab[mkji - 1];
			if (nm >= 1) {

/*     -- record the atoms in inov that neighbor atom ir */

			    i__2 = nm;
			    for (m = 1; m <= i__2; ++m) {
				in = natm[m + mkji * 150 - 151];
				if (in != ir) {
				    ++io;
				    if (io > ict) {
					s_stop("SOLVA_ERROR: intrsctns > max",
						 (ftnlen)28);
				    }
				    dx[io - 1] = xr - xyz[in + xyz_dim1];
				    dy[io - 1] = yr - xyz[in + (xyz_dim1 << 1)
					    ];
/* Computing 2nd power */
				    r__1 = dx[io - 1];
/* Computing 2nd power */
				    r__2 = dy[io - 1];
				    dsq[io - 1] = r__1 * r__1 + r__2 * r__2;
				    d__[io - 1] = sqrt(dsq[io - 1]);
				    inov[io - 1] = in;
				}
			    }
			}
		    }
		}
	    }
	}
L14:
	if (io >= 1) {

/*     z resolution determined */

	    zres = rrx2 / nzp;
	    zgrid = xyz[ir + xyz_dim1 * 3] - rr - zres / (float)2.;
	} else {
	    area = pix2 * rrx2;
	    goto L18;
	}

/*     section atom spheres perpendicular to the z axis */

	i__2 = nzp;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    zgrid += zres;

/*     find the radius of the circle of intersection of */
/*     the ir sphere on the current z-plane */

/* Computing 2nd power */
	    r__1 = zgrid - zr;
	    rsec2r = rrsq - r__1 * r__1;
	    rsecr = sqrt(rsec2r);
	    i__3 = karc;
	    for (k = 1; k <= i__3; ++k) {
		arci[k - 1] = (float)0.;
	    }
	    karc = 0;
	    i__3 = io;
	    for (j = 1; j <= i__3; ++j) {
		in = inov[j - 1];

/*     find radius of circle locus */

/* Computing 2nd power */
		r__1 = zgrid - xyz[in + xyz_dim1 * 3];
		rsec2n = radsq[in - 1] - r__1 * r__1;
		if (rsec2n <= (float)0.) {
		    goto L10;
		}
		rsecn = sqrt(rsec2n);

/*     find intersections of n.circles with ir circles in section */

		if (d__[j - 1] >= rsecr + rsecn) {
		    goto L10;
		}

/*     do the circles intersect, or is one circle completely inside the other? */

		b = rsecr - rsecn;
		if (d__[j - 1] > dabs(b)) {
		    goto L20;
		}
		if (b <= (float)0.) {
		    goto L9;
		}
		goto L10;

/*     if the circles intersect, find the points of intersection */

L20:
		++karc;
		if (karc >= ict) {
		    s_stop("SOLVA_ERROR: max intersections exceeded2", (
			    ftnlen)40);
		}

/*     Initial and final arc endpoints are found for the ir circle intersected */
/*     by a neighboring circle contained in the same plane. The initial endpoint */
/*     of the enclosed arc is stored in arci, and the final arc in arcf */
/*     law of cosines */

		trig_test__ = (dsq[j - 1] + rsec2r - rsec2n) / (d__[j - 1] * (
			float)2. * rsecr);
		if (trig_test__ >= (float)1.) {
		    trig_test__ = (float).99999;
		}
		if (trig_test__ <= (float)-1.) {
		    trig_test__ = (float)-.99999;
		}
		alpha = acos(trig_test__);

/*     alpha is the angle between a line containing a point of intersection and */
/*     the reference circle center and the line containing both circle centers */

		beta = atan2(dy[j - 1], dx[j - 1]) + pi;

/*     beta is the angle between the line containing both circle centers and the x-axis */

		ti = beta - alpha;
		tf = beta + alpha;
		if (ti < (float)0.) {
		    ti += pix2;
		}
		if (tf > pix2) {
		    tf -= pix2;
		}
		arci[karc - 1] = ti;
		if (tf >= ti) {
		    goto L3;
		}

/*     if the arc crosses zero, then it is broken into two segments. */
/*     the first ends at pix2 and the second begins at zero */

		arcf[karc - 1] = pix2;
		++karc;
L3:
		arcf[karc - 1] = tf;
L10:
		;
	    }

/*     find the accessible surface area for the sphere ir on this section */

	    if (karc != 0) {
		goto L19;
	    }
	    arcsum = pix2;
	    goto L25;

/*     The arc endpoints are sorted on the value of the initial arc endpoint */

L19:
	    sortag_(arci, &karc, tag);

/* *************************************** */
/*     calculate the accessible area */
/* *************************************** */

	    arcsum = arci[0];
	    t = arcf[tag[0] - 1];
	    if (karc == 1) {
		goto L11;
	    }
	    i__3 = karc;
	    for (k = 2; k <= i__3; ++k) {
		if (t < arci[k - 1]) {
		    arcsum = arcsum + arci[k - 1] - t;
		}
		tt = arcf[tag[k - 1] - 1];
		if (tt > t) {
		    t = tt;
		}
	    }
L11:
	    arcsum = arcsum + pix2 - t;

/*     The area/radius is equal to the accessible arc length x the section thickness. */

L25:
	    parea = arcsum * zres;

/*     Add the accessible area for this atom in this section to the area for this */
/*     atom for all the section encountered thus far */

	    area += parea;
L9:
	    ;
	}

/*     scale area to vdw shell */

L18:
	b = area * rr;
	accs[ir] = b;
/* ------------------------------------------------------------------ */
/* The following line converts from accessible to contact surface */
/*         c=(b*(rad(ir)-probe)**2)/(rad(ir)**2) */
/* ------------------------------------------------------------------ */
    }
    s_wsfe(&io___163);
    do_fio(&c__1, " SOLVA: PROGRAM ENDS CORRECTLY", (ftnlen)30);
    e_wsfe();
    return 0;
} /* solva_ */

/* Subroutine */ int sortag_(a, n, tag)
real *a;
integer *n, *tag;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, l, m;
    static real t;
    static integer ij, il[16], tg, iu[16];
    static real tt;

    /* Parameter adjustments */
    --tag;
    --a;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	tag[i__] = i__;
    }
    m = 1;
    i__ = 1;
    j = *n;
L5:
    if (i__ >= j) {
	goto L70;
    }
L10:
    k = i__;
    ij = (j + i__) / 2;
    t = a[ij];
    if (a[i__] <= t) {
	goto L20;
    }
    a[ij] = a[i__];
    a[i__] = t;
    t = a[ij];
    tg = tag[ij];
    tag[ij] = tag[i__];
    tag[i__] = tg;
L20:
    l = j;
    if (a[j] >= t) {
	goto L40;
    }
    a[ij] = a[j];
    a[j] = t;
    t = a[ij];
    tg = tag[ij];
    tag[ij] = tag[j];
    tag[j] = tg;
    if (a[i__] <= t) {
	goto L40;
    }
    a[ij] = a[i__];
    a[i__] = t;
    t = a[ij];
    tg = tag[ij];
    tag[ij] = tag[i__];
    tag[i__] = tg;
    goto L40;
L30:
    a[l] = a[k];
    a[k] = tt;
    tg = tag[l];
    tag[l] = tag[k];
    tag[k] = tg;
L40:
    --l;
    if (a[l] > t) {
	goto L40;
    }
    tt = a[l];
L50:
    ++k;
    if (a[k] < t) {
	goto L50;
    }
    if (k <= l) {
	goto L30;
    }
    if (l - i__ <= j - k) {
	goto L60;
    }
    il[m - 1] = i__;
    iu[m - 1] = l;
    i__ = k;
    ++m;
    goto L80;
L60:
    il[m - 1] = k;
    iu[m - 1] = j;
    j = l;
    ++m;
    goto L80;
L70:
    --m;
    if (m == 0) {
	return 0;
    }
    i__ = il[m - 1];
    j = iu[m - 1];
L80:
    if (j - i__ >= 1) {
	goto L10;
    }
    if (i__ == 1) {
	goto L5;
    }
    --i__;
L90:
    ++i__;
    if (i__ == j) {
	goto L70;
    }
    t = a[i__ + 1];
    if (a[i__] <= t) {
	goto L90;
    }
    tg = tag[i__ + 1];
    k = i__;
L100:
    a[k + 1] = a[k];
    tag[k + 1] = tag[k];
    --k;
    if (t < a[k]) {
	goto L100;
    }
    a[k + 1] = t;
    tag[k + 1] = tg;
    goto L90;
} /* sortag_ */

/* Subroutine */ int summer_(sname, slen, nats, accs, m1, m2, atomtype, 
	resindex, resnam, ressums, tsums, sname_len, resnam_len)
char *sname;
integer *slen, *nats;
real *accs;
integer *m1, *m2, *atomtype, *resindex;
char *resnam;
real *ressums, *tsums;
ftnlen sname_len;
ftnlen resnam_len;
{
    /* System generated locals */
    integer ressums_dim1, ressums_offset, i__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe(), s_cmp();
    /* Subroutine */ int s_copy();
    integer s_rsfi(), e_rsfi(), f_clos();

    /* Local variables */
    static char line[256];
    static integer ires, maxr, i__, j;
    static char acids[3*100];
    extern integer fopen_();
    static logical stand;
    static real standarea[500]	/* was [100][5] */;
    extern /* Subroutine */ int which3_();
    static logical ok;
    static integer ir, nacids, rindex[2000];
    extern integer readstring_();
    static integer len;
    static char res[3];

    /* Fortran I/O blocks */
    static cilist io___176 = { 0, 3, 0, "(4a)", 0 };
    static icilist io___181 = { 0, line+16, 0, "(f7.2)", 7, 1 };
    static icilist io___183 = { 0, line+29, 0, "(f7.2)", 7, 1 };
    static icilist io___184 = { 0, line+42, 0, "(f7.2)", 7, 1 };
    static icilist io___185 = { 0, line+55, 0, "(f7.2)", 7, 1 };
    static icilist io___186 = { 0, line+68, 0, "(f7.2)", 7, 1 };
    static cilist io___187 = { 0, 4, 0, "(a,i3,a)", 0 };
    static cilist io___188 = { 0, 4, 0, "(a)", 0 };


/*     -- 	program to sum atomic accessibilities by residue. */
/*     --	copes with atom and hetatom records */
/*     --	produces relative accessibilities for the 20 common aminos */
/*     --	ouput written to .rsa file (channel 4) */

/*     if "standard.data" exists in, read them in. */

    /* Parameter adjustments */
    --resindex;
    --atomtype;
    --accs;
    ressums_dim1 = *m2;
    ressums_offset = 1 + ressums_dim1 * 6;
    ressums -= ressums_offset;
    resnam -= 10;
    --tsums;

    /* Function Body */
    stand = FALSE_;
    if (fopen_(&c__1, sname, slen, "old", (ftnlen)256, (ftnlen)3) != 0) {
	s_wsfe(&io___176);
	do_fio(&c__1, "REM  Relative accessibilites read from", (ftnlen)38);
	do_fio(&c__1, " external file \"", (ftnlen)16);
	do_fio(&c__1, sname, (*slen));
	do_fio(&c__1, "\"", (ftnlen)1);
	e_wsfe();
	stand = TRUE_;
	i__ = 0;
	while(readstring_(&c__1, line, &len, (ftnlen)256) >= 0 && (real) i__ <
		 (float)100.) {
	    if (s_cmp(line, "ATOM", (ftnlen)4, (ftnlen)4) == 0) {
		++i__;
		s_copy(acids + (i__ - 1) * 3, line + 12, (ftnlen)3, (ftnlen)3)
			;
		s_rsfi(&io___181);
		do_fio(&c__1, (char *)&standarea[i__ - 1], (ftnlen)sizeof(
			real));
		e_rsfi();
		s_rsfi(&io___183);
		do_fio(&c__1, (char *)&standarea[i__ + 99], (ftnlen)sizeof(
			real));
		e_rsfi();
		s_rsfi(&io___184);
		do_fio(&c__1, (char *)&standarea[i__ + 199], (ftnlen)sizeof(
			real));
		e_rsfi();
		s_rsfi(&io___185);
		do_fio(&c__1, (char *)&standarea[i__ + 299], (ftnlen)sizeof(
			real));
		e_rsfi();
		s_rsfi(&io___186);
		do_fio(&c__1, (char *)&standarea[i__ + 399], (ftnlen)sizeof(
			real));
		e_rsfi();
	    }
	}
	s_wsfe(&io___187);
	do_fio(&c__1, " RELATIVE (STANDARD) ACCESSIBILITIES READFOR ", (
		ftnlen)45);
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
	do_fio(&c__1, " AMINO ACIDS", (ftnlen)12);
	e_wsfe();
	cl__1.cerr = 0;
	cl__1.cunit = 1;
	cl__1.csta = 0;
	f_clos(&cl__1);
    } else {
	s_wsfe(&io___188);
	do_fio(&c__1, " NO STANDARD VALUES INPUT", (ftnlen)25);
	e_wsfe();
    }
    nacids = i__;
    i__1 = resindex[*nats];
    for (i__ = 1; i__ <= i__1; ++i__) {
	rindex[i__ - 1] = 0;
	if (stand) {
	    s_copy(res, resnam + i__ * 10, (ftnlen)3, (ftnlen)3);
	    maxr = (float)100.;
	    which3_(res, acids, &ires, &nacids, &maxr, &ok, (ftnlen)3, (
		    ftnlen)3);
	    if (ok) {
		rindex[i__ - 1] = ires;
	    }
	}
    }

/*     -- sum the values */

    i__1 = *nats;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ir = resindex[i__];
	tsums[1] += accs[i__];
	ressums[ir + ressums_dim1 * 6] += accs[i__];
	if (rindex[ir - 1] != 0) {
	    if (atomtype[i__] == 3) {
		ressums[ir + ressums_dim1 * 10] += accs[i__];
		tsums[5] += accs[i__];
	    } else if (atomtype[i__] != 0) {
		ressums[ir + ressums_dim1 * 9] += accs[i__];
		tsums[4] += accs[i__];
		if (atomtype[i__] == 1) {
		    ressums[ir + ressums_dim1 * 7] += accs[i__];
		    tsums[2] += accs[i__];
		} else if (atomtype[i__] == 2) {
		    ressums[ir + (ressums_dim1 << 3)] += accs[i__];
		    tsums[3] += accs[i__];
		}
	    }
	}
    }

/*     -- calculate realtive accessibilities */

    i__1 = resindex[*nats];
    for (i__ = 1; i__ <= i__1; ++i__) {
	ires = rindex[i__ - 1];
	if (stand && ires != 0) {
	    for (j = 1; j <= 5; ++j) {
		if (standarea[ires + j * 100 - 101] > (float)0.) {
		    ressums[i__ + (j + 10) * ressums_dim1] = ressums[i__ + (j 
			    + 5) * ressums_dim1] * (float)100. / standarea[
			    ires + j * 100 - 101];
		} else {
		    ressums[i__ + (j + 10) * ressums_dim1] = (float)0.;
		}
	    }
	} else {
	    for (j = 1; j <= 5; ++j) {
		ressums[i__ + (j + 10) * ressums_dim1] = (float)-99.9;
	    }
	}
    }

/*     -- COMPLETED */

    return 0;
} /* summer_ */


/*     --what atom is it ? ie non-polar/polar, main/sidechain. */

integer what_atom__(atom, atom_len)
char *atom;
ftnlen atom_len;
{
    /* Initialized data */

    static char phobs[4*20+1] = " CA  CB  CD  CD1 CD2 CG  CG1 CG2 CE  CE1 CE\
2 CE3 CH2 CZ  CZ2 CZ3 SD  SG ********";
    static char phils[4*20+1] = " AD1 AD2 AE1 AE2 ND1 ND2 NE  NE1 NE2 NH1 NH\
2 NZ  OD1 OD2 OE3 OE2 OE1 OG  OG1 OH ";
    static char mc[4*4+1] = " N   C   O   OXT";

    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer s_cmp();

    /* Local variables */
    static integer i__;


/* 	atoms classed as non-polar in sidechains */


/*     ATOMS CLASSED AS POLAR IN SIDECHAINS */


/*       MAIN CHAIN ATOMS */

    ret_val = 3;
    for (i__ = 1; i__ <= 4; ++i__) {
	if (s_cmp(atom, mc + (i__ - 1 << 2), (ftnlen)4, (ftnlen)4) == 0) {
	    return ret_val;
	}
    }
    ret_val = 1;
    for (i__ = 1; i__ <= 20; ++i__) {
	if (s_cmp(atom, phobs + (i__ - 1 << 2), (ftnlen)4, (ftnlen)4) == 0) {
	    return ret_val;
	}
    }
    ret_val = 2;
    for (i__ = 1; i__ <= 20; ++i__) {
	if (s_cmp(atom, phils + (i__ - 1 << 2), (ftnlen)4, (ftnlen)4) == 0) {
	    return ret_val;
	}
    }
    ret_val = 0;
    if (*(unsigned char *)&atom[1] == 'O' || *(unsigned char *)&atom[1] == 
	    'N') {
	ret_val = 2;
    } else {
	ret_val = 3;
    }
    return ret_val;
} /* what_atom__ */

/* Subroutine */ int vanin_(vname, vlen, nacids, aacids, anames, numats, 
	vradii, maxr, maxa, vname_len, aacids_len, anames_len)
char *vname;
integer *vlen, *nacids;
char *aacids, *anames;
integer *numats;
real *vradii;
integer *maxr, *maxa;
ftnlen vname_len;
ftnlen aacids_len;
ftnlen anames_len;
{
    /* System generated locals */
    integer vradii_dim1, vradii_offset, anames_dim1, anames_offset;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop();
    integer s_cmp();
    /* Subroutine */ int s_copy();
    integer f_clos();

    /* Local variables */
    static char card[256], c__[256*256];
    static integer l[256], n, nlabs;
    extern integer fopen_(), parse_();
    extern doublereal readfloat_();
    extern integer readstring_();
    static integer len;


/* -- Read in van der Waal radii from external file "vfile" */
/* -- nacids = number of residues read in */
/* -- aacids = *3 character array containing amino acid names */
/* -- anames = *4 character array containing atom names for each residue */
/* -- numats = array containing number of atoms for each residue */
/* -- radii  = vdw radii for each atom, indexed identically to anames */

    /* Parameter adjustments */
    --numats;
    aacids -= 3;
    vradii_dim1 = *maxr;
    vradii_offset = 1 + vradii_dim1 * 1;
    vradii -= vradii_offset;
    anames_dim1 = *maxr;
    anames_offset = 1 + anames_dim1 * 1;
    anames -= anames_offset * 4;

    /* Function Body */
    if (fopen_(&c__10, vname, vlen, "old", (ftnlen)256, (ftnlen)3) == 0) {
	s_stop("ERROR: unable to open \"vdw.radii\"", (ftnlen)33);
    }
    *nacids = 0;
    nlabs = 0;
    while(readstring_(&c__10, card, &len, (ftnlen)256) >= 0) {
	n = parse_(card, &len, " ", c__, l, (ftnlen)256, (ftnlen)1, (ftnlen)
		256);
	if (s_cmp(c__, "RESIDUE", (ftnlen)256, (ftnlen)7) == 0) {
	    ++(*nacids);
	    if (*nacids > *maxr) {
		s_stop("ERROR: increase max_r", (ftnlen)21);
	    }
	    s_copy(aacids + *nacids * 3, c__ + 512, (ftnlen)3, (ftnlen)3);
	    numats[*nacids] = 0;
	}
	if (s_cmp(c__, "ATOM", (ftnlen)256, (ftnlen)4) == 0) {
	    ++numats[*nacids];
	    if (numats[*nacids] > *maxa) {
		s_stop("ERROR: increase max_a", (ftnlen)21);
	    }
	    s_copy(anames + (*nacids + numats[*nacids] * anames_dim1 << 2), 
		    card + 5, (ftnlen)4, (ftnlen)4);
	    vradii[*nacids + numats[*nacids] * vradii_dim1] = readfloat_(c__ 
		    + (n - 1 << 8), &l[n - 1], (ftnlen)256);
	}
    }
    cl__1.cerr = 0;
    cl__1.cunit = 10;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
} /* vanin_ */

/* Main program alias */ int accall_ () { MAIN__ (); }
