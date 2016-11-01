#ifndef __PDB_H__
#define __PDB_H__

/* ==== HEADER pdb.h ==== */

/* PDB-handling function collection.  */

/* ANSI C cc SGI Iris Indigo IRIX 4.0.5  30. Nov. 1993. Andris */

/* ---- INCLUDE FILES ---- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h> /* SEEK_SET on Suns */

/* ---- GLOBAL CONSTANTS ---- */

#define MAXCHAINSIZE 12000  /* max. no. of atoms per chain */
#define MAXCHAINNO 128       /* max. no. of chains per PDB entry */

#define FALSE 0             /* logical values */
#define TRUE 1

#define ALLATOMS 0     /* read all atoms */
#define CALPHA 1       /* read C-alpha backbone only */
#define RELAXED 0      /* read non-standard AA-s ('X') */
#define STRICT 1       /* read the 20 standard AA-s + B,Z only */

/* A global array for converting 1-letter AA codes to 3-letter ones.
 * Use the construct Aacode[Onelettercode-'A']
 */
extern char *Aacode[26];
	    
/* ---- GLOBAL PDB TYPES ---- */

typedef char Str3_[4];  /* 3-letter words... */

typedef struct    /* entry for an atom */
{
  int Atno;	/* atom serial number as in the PDB entry */
  Str3_ Id;  /* atom type like CA or OD2: 3 chars max. + \0 */
  char Aa;     /* 1-letter AA code or ' ' for OXT */
  int Resno;  /* residue no. within a chain */
  char Rid;   /* insertion code like in "27A" */
  float X,Y,Z; /* the coordinates */
  float Occu, Bfact; /* occupancy and temperature factor */
} Atom_ ;

typedef struct	    /* entry for a H-bond */
{
    int At;	    /* the atom */
    int Don, Acc;   /* atom nos of donor and acceptor to At */
}
Hbond_ ;

typedef struct	    /* entry for a secondary structure element */
{
    char Helix;	    /* 1 if helical, 0 if sheet */
    int No;	    /* serial number of structure */
    Str3_ Id;	    /* identifier string */
    int Strandno;   /* for sheets only */
    int Beg, End;   /* residue nos (1..N) */
    char Chid;	    /* chain ID, for begin only */
    char Begaa, Endaa;	/* amino acid 1-letter codes */
    char Begrid, Endrid;    /* insertion codes */
    int Type;	    /* 1..7 for helices, -1..1 for sheets */
	    /* these entries are for sheets only */
    Str3_ Thisat, Otherat;  /* atom names for beta registration */
    char Thisaa, Otheraa;   /* amino acid codes for beta regs */
    int This, Other;	/* registration AA positions */
    char Thisrid, Otherid, Otherchid;	/* ins codes and chain ID for the OTHER */
}
Secstr_ ;

typedef struct   /* entry for a chain */
{
  Atom_ *Chhd; /* array of atoms in this chain */
  Hbond_ *Hbonds;   /* array of H-bonds */
  int Hbno;	/* length of Hbonds[] */
  Secstr_ *Secs;    /* secondary structure info */
  int Secno;	/* length of Secs[] */
  int Atomno;   /* no. of atoms in chain */
  int Aano;     /* no. of amino acids in chain */
  char Chid;    /* chain ID (for multichains) */
  char Type;    /* chain type: 'P'-rotein, 'A'-lpha or 'X' */
} Chain_ ;

/* ---- PROTOTYPES ---- */

#ifdef __cplusplus
extern "C" {
#endif

int read_pdb(FILE *Pdb, Chain_ Chains[],int Ca,int Strict);
char aa_code31(const char *Aa3);
char *pdb_seq(Chain_ *Chain);
float atom_dist(Atom_ *At1, Atom_ *At2);
void write_pdb(FILE *Pdb, Chain_ Chains[], int Chainno,
	const char Pdbcode[5], const char *Header, const char *Compound,
	char *Remarks[], int Remno);
void remove_chain(Chain_ *Chain);

#ifdef __cplusplus
}
#endif

/* ==== END OF HEADER pdb.h ==== */

#endif

