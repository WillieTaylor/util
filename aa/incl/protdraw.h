#ifndef __PROTDRAW_H__
#define __PROTDRAW_H__

/* ==== HEADER protdraw.h ==== */

/* Draws a colour-coded PDB chain which can
 be interactively rotated/moved on the screen by the mouse. */

/* PORTABILITY WARNING: Highly Silicon Graphics-specific code! */
/* Link in the SGI shared graphics library (-lgl_s) into executable */

/* ANSI C (otherwise...), Silicon Graphics Iris Indigo, IRIX 5.2 */

/* 26. July 1994, Andris Aszodi */

/* ---- DEFINITIONS ---- */

/* 32 bit colour constants, full brightness. The alpha byte is set to
 FF in all cases (no transparency). The byte order is ABGR. */
#define RGB_WHITE 0xFFFFFFFF
#define RGB_YELLOW 0xFF00FFFF
#define RGB_MAGENTA 0xFFFF00FF
#define RGB_CYAN 0xFFFFFF00
#define RGB_RED 0xFF0000FF
#define RGB_ORANGE 0xFF0080FF
#define RGB_GREEN 0xFF00FF00
#define RGB_BLUE 0xFFFF0000
#define RGB_BLACK 0xFF000000

/* default window position */
#define CXORIG 100
#define CYORIG 10

/* symbolic constants for manipulating the static variables in
 * "pdbdraw.c" (menu options etc.) See set_pdrvars()
 */
#define PDR_OPTNO 11	    /* no. of Pdropt_ constants */
#define PDR_ON 1	    /* switches */
#define PDR_OFF 0
typedef enum {PDR_DEPTH, PDR_CALPHA, PDR_SIDECH, 
	    PDR_HBOND, PDR_VIEW, PDR_LCD, 
	    PDR_INV, PDR_COL, PDR_LABEL} Pdropt_ ;
typedef enum {VIEW_ORTHO=1, VIEW_WIDE, VIEW_NORM,
	    VIEW_TELE} Viewopt_ ;
typedef enum {COLOUR_LEVITT=1, COLOUR_BINPHOB, COLOUR_OCCU,
	    COLOUR_BFACT, COLOUR_OCBF, COLOUR_CHAIN, COLOUR_POS, COLOUR_USER} Colopt_ ;

/* ---- TYPES ---- */

/* Point_ : an atom's description on the screen */
typedef struct
{
    float X[3];	    /* X,Y,Z coordinates */
    unsigned long Col;	/* RGB colour info */
    float Occu, Bfact;	/* occupancy and B-factor */
    Str4_ Id;	/* atom ID, in PDB style (see "pdbprot.h") */
    int Atno;	/* PDB atom number */
}
Point_ ;

/* Res_: a residue description */
typedef struct
{
    Point_ *First;  /* points to the first atom in it */
    Point_ *Pca, *Pcb;	/* points to the C-alpha and C-beta */
    int Atomno;	    /* total no. of atoms: backbone+sidechain */
    char Aa, Chid, Rid;	    /* amino acid code, chain ID, ins code */
    int Resno;	    /* PDB residue number */
}
Res_ ;

/* Pair_ : record for H-bonds and disulfide bridges */
typedef enum {HBOND, DISULF} Pairtype_;
typedef struct
{
    Pairtype_ Ptype;	/* H-bond or disulfide */
    const Point_ *P1, *P2;	/* pointers to the two endpoints */
}
Pair_ ;

/* Subunit_ : holds draw information for a whole chain,
 * complete with residue sequence, intrachain H-bond
 * and S-S bonds (similar to Chain_ in "pdbprot.h") 
 */
typedef struct
{
    Point_ *Points;  /* list of atoms */
    int Pno;	    /* no. of points */
    Res_ *Res;	/* list of residues */
    int Rno;	    /* no. of residues */
    Pair_ *Hpairs;  /* intrachain H-bonds array */
    int Hpno;	    /* no. of H-bonds */
    Pair_ *Spairs;  /* intrachain disulfides array */
    int Spno;	    /* no. of disulfides */
}
Subunit_ ;

/* Drawentry_ : draw information for a molecule that may
 * contain several subunits and interchain H-bonds and
 * disulfides. Only polypeptides are supported as yet.
 */
typedef struct
{
    Subunit_ *Subs;	/* array of chains */
    int Subno;	    /* no. of chains */
    Pair_ *Hpairs;  /* interchain H-bonds */
    int Hpno;	    /* no. of H-bonds */
    Pair_ *Spairs;  /* interchain disulfides */
    int Spno;	    /* no. of disulfides */
    float Rot[3][3];	/* overall rotation matrix */
    int Atot;	    /* total no. of atoms */
    int Status[PDR_OPTNO];	/* draw status vector */
    float Lcd;	/* largest coordinate */
    double Lowval, Upval;	/* colour scale limits */
    const unsigned long *Col;	/* pointer to colour code array */
}
Drawentry_ ;

/* ---- GLOBAL VARS ---- */

extern unsigned long Phob2col[26];
extern unsigned long Levittcol[26];

/* ---- PROTOTYPES ---- */

Drawentry_ *prot_draw(const Pdbentry_ *Prot, 
	const unsigned long *Col, float Bval, float Yval);
void destruct_draw(Drawentry_ *Draw);
long init_draw(Drawentry_ *Draw, char *Windowtitle, 
	long Xorig, long Yorig);
void display_draw(long Gid, Drawentry_ *Draw, Device Dev);
void draw_chains(long Gid, Drawentry_ *Draw);
void set_pdrvars(Drawentry_ *Draw, Pdropt_ Pdropt, int Value);

/* ==== END OF HEADER protdraw.h ==== */

#endif
