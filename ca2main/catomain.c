/*************************************************************************

   Program:    
   File:       
   
   Version:    
   Date:       
   Function:   
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1996
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 419 3890
   EMail:      INTERNET: martin@biochem.ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdio.h>

/************************************************************************/
/* Defines and macros
*/
#define TRUE 1
#define FALSE 0
#define MAXBUFF 160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
int ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);
extern void atomg_(char *);
extern void caord_(char *);

/************************************************************************/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   
   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      atomg_(infile);
      caord_(outfile);
   }
   else
   {
      Usage();
   }
}

/************************************************************************/
/*>int ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
   Returns: int                 Success?

   Parse the command line
   
   14.11.96 Original    By: ACRM
*/
int ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = ' ';
   infile[1] = outfile[1] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}

/************************************************************************/
void Usage(void)
{
   fprintf(stderr,"\ncatomain (c) Willie Taylor, modified by Andrew C.R. \
Martin, UCL\n");
   fprintf(stderr,"\nUsage: catomain [in.pdb [out.pdb]]\n");
   fprintf(stderr,"\ncatomain adds backbone atoms to a CA PDB file using \
a modified version\n");
   fprintf(stderr,"of the method of Levitt (1977).\n\n");
}

