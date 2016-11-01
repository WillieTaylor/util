/* ==== FUNCTIONS dsspca.c ==== */

/* Little program to get secondary structure info and C-alpha coordinates
   from DSSP files. H-bonds and SS-bridges are not read. */

/* ANSI C, IRIX 5.3, 16. May 1995. Andris Aszodi */

/* ---- STANDARD HEADERS ---- */

#include <string.h>
#include <math.h>

/* ---- MODULE HEADER ---- */

#include "dsspca.h"

/* ---- FUNCTIONS ---- */

/* dssp_ca: reads the (opened) Dssp stream. This is a text file produced
   by DSSP, Version Oct. 1988 (ref: W. Kabsch, C. Sander, Biopolymers
   22:2577-2637 (1983)). Returns Entry, an array of size Nres, containing
   chain ID, residue no, 1-letter AA code, Kabsch/Sander secondary
   structure code and solvent accessibility for each residue. 
   Returns the size of Entry (i.e. no. of AAs) or 0 on error.
*/
#define LINELEN 132

int dssp_ca(FILE *Dssp, Dssprec_ **Entry, int *Chainno)
{
    char Line[LINELEN],Dummy[18];  /* read-in buffers */
    Dssprec_ *Cur;            /* temp pointer */
    int i,Entryj,Nres;
    char c;
    
    /* init reading, skip lines */
    Nres=*Chainno=0;
    *Entry=NULL;

    rewind(Dssp);
    /* read total no. of residues, no. of chains */
    do
    {
	if(NULL==fgets(Line,LINELEN,Dssp)) return(0);
	if (NULL!=strstr(Line, "TOTAL NUMBER OF RESIDUES"))
	    break;
    } while(1); 
    if (2>sscanf(Line,"%d %d",&Nres,Chainno)) return(0);
    
#if 0
    /* skip total accessible surf */
    if (NULL==fgets(Line,LINELEN,Dssp)) return(0);
    
    /* skip total no. of hydrogen bonds */
    if (NULL==fgets(Line,LINELEN,Dssp)) return(0);
#endif
    
    /* skip all info until the last text line is encountered */
    do 
    {
	if(NULL==fgets(Line,LINELEN,Dssp)) return(0);
	if(NULL!=strstr(Line, "#  RESIDUE AA STRUCTURE"))
	    break;
    }
    while (1);
        
    /* adjust storage, allocate */
    Nres+=*Chainno-1;
    *Entry= (Dssprec_ *) calloc(Nres,sizeof(Dssprec_));

    while (NULL!=fgets(Line,LINELEN,Dssp)) 
    {
	sscanf(Line,"%d%1s",&Entryj,&c);  /* reads pos and 1st non-wspace */
	Entryj--;   /* indexes res position now */
        Cur=*Entry+Entryj;  /* current pos in array */
	if (c=='!')
	{                   /* chain break */
	    Cur->Resno=0; Cur->Res=c;
	    Cur->Chain=Cur->Secstruct=Cur->Turns3=Cur->Turns4=Cur->Turns5=' ';
	    Cur->Bend=Cur->Chir=Cur->Bridge1=Cur->Bridge2=Cur->Sheet=' ';
	    Cur->Access=0;
            continue;              /* nothing else to be done */
	}
	else 
	{                   /* normal entry */
            if (*Chainno>1)
	    {
                   /* multiple chains, chain ID char read, too */
	    	sscanf(Line,"%*d%d%1s%1s",
                       &(Cur->Resno),&(Cur->Chain),&(Cur->Res));
	    }
	    else
	    {      /* there's no chain char if Chainno==1 */
	    	sscanf(Line,"%*d%d%1s",&(Cur->Resno),&(Cur->Res));
		Cur->Chain=' ';
	    }
	    
	    /* transfer the funny char descriptors */
	    Cur->Secstruct=Line[16];
	    Cur->Turns3=Line[18]; Cur->Turns4=Line[19]; Cur->Turns5=Line[20];
	    Cur->Bend=Line[21]; Cur->Chir=Line[22];
	    Cur->Bridge1=Line[23]; Cur->Bridge2=Line[24];
	    sscanf(Line+25, "%d%d", &(Cur->Beta1), &(Cur->Beta2));
	    Cur->Sheet=Line[33];
	    
	    /* read accessibility and C-alpha coordinates */
            sscanf(Line+34,
              "%d%*d,%*f%*d,%*f%*d,%*f%*d,%*f%*f%*f%*f%*f%*f%f%f%f",
              &(Cur->Access),&(Cur->Xca),&(Cur->Yca),&(Cur->Zca));
	    
            /* Rename half-cystines. These are coded as 'a'..'z' */
            if (Cur->Res>='a' && Cur->Res<='z')
		Cur->Res='C';
	}
    }        /* while */
    return(Nres);
    
}     /* END of dssp_ca */

#undef LINELEN

/* dssp_cadist(): returns the CA:CA distance between two DSSP records
 * pointed to by Dp1 and Dp2. Checks if they're NULL.
 */
double dssp_cadist(const Dssprec_ *Dp1, const Dssprec_ *Dp2)
{
    double Dist, D;
    
    if (Dp1==NULL || Dp2==NULL)
    {
	fprintf(stderr, "\n? dssp_cadist(%x, %x): NULL argument\n", Dp1, Dp2);
	return(0.0);
    }
    
    D=Dp1->Xca-Dp2->Xca; Dist=D*D;
    D=Dp1->Yca-Dp2->Yca; Dist+=D*D;
    D=Dp1->Zca-Dp2->Zca; Dist+=D*D;
    return(sqrt(Dist));
}
/* END of dssp_cadist() */

/* ==== END OF FUNCTIONS dsspca.c ==== */    
