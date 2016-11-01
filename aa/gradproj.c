/* ==== FUNCTIONS gradproj.c ==== */

/* Gradual Projection algorithm. */

/* ANSI C, IRIX 5.3, 27. Nov. 1995. */

/* ---- HEADER ---- */

#include "incl/gradproj.h"

/* ---- DEFINITIONS ---- */

/* Take care of the silly SunOS4.1 non-ANSI math library */
#ifdef __SUN_OS__
#define sqrtf sqrt
#define expf exp
#endif

/* ---- PROTOTYPES ---- */

static void dist_metric(Trimat_ Dist, Trimat_ Metric, double Cdist2[], int Rno);
static void metric_dist(Trimat_ Metric, Trimat_ Dist, int Rno);
static float trieq_bal(Trimat_ Metric, int Rno, float *Diagshf);

static int metric_project(Trimat_ Metric, int Rno, Sqmat_ Xyz,
	float Diagshf, float Equ, int Oldim, float *Projqual);
static int centre_dist(Trimat_ Dist, int Rno, double Cdist2[]);
static void new_distmat(Sqmat_ Xyz, int Rno, int Dim, Trimat_ Newdist);

/* ==== FUNCTIONS ==== */

/* ---- GRADUAL PROJECTION ---- */

/* grad_proj(): performs a gradual projection iteration down to
 * 3D. The initial UN-squared distances are in Dist, the corresponding
 * strictness values in Strict (both are triangular matrices with
 * size Size x Size). In each iteration, the distances in Dist
 * are applied to the actual distance matrix with appropriate
 * strictnesses. The final 3D embedding is obtained in Xyz
 * (should be an Rno x Rno square matrix).
 * Return value: the number of iterations done.
 */
int grad_proj(const Trimat_ Dist, const Trimat_ Strict, int Rno, Sqmat_ Xyz)
{
    static const int MAX_TRINEQ=10;
    
    Trimat_ D, M;
    double *Cdist2;
    register unsigned int i, j, Dim, Oldim, Itno;
    int Cviol, Projno;
    float Diagshf, Projqual, Tricorr;
    
    /* set up storage */
    D=alloc_trimat(Rno); M=alloc_trimat(Rno);
    Cdist2=(double *) calloc(Rno, sizeof(double));

    /* init the actual distance matrix D */
    for (i=0; i<Rno; i++) {
	for (j=0; j<=i; j++) {
	    Dist[i][j] *= Dist[i][j];
	    D[i][j]=Dist[i][j];
	}
    }
    
    /* projection cycle */
    
    for (Dim=Oldim=Rno, Projno=1; Dim>3; Projno++,Oldim=Dim)
    {
	/* triangle inequality balancing */
	Itno=0; Tricorr=0.0;
	do
	{
	    if (Tricorr>0.0)	/* get back D */
		metric_dist(M,D,Rno);
	
	    /* regularisation: iron out metric matrix triangle ineq's */
	    Cviol=centre_dist(D,Rno,Cdist2);
	    dist_metric(D,M,Cdist2,Rno);
	    Tricorr=trieq_bal(M,Rno,&Diagshf);
	}
	while (++Itno<MAX_TRINEQ && (Tricorr>0.0 || Cviol));
	dist_metric(D,M,Cdist2,Rno);   /* final */
	
	/* embed the metric matrix: keep all pos. eigenvalues */
	Dim=metric_project(M, Rno, Xyz, Diagshf, 0.99, Oldim, &Projqual);
	printf("Projection dimension = %3d  (quality = %f)\n", Dim, Projqual);
	if (Dim==3) break;
	
	/* generate new dist matrix, blend in desired values */
	new_distmat(Xyz, Rno, Dim, D);
	for (i=0; i<Rno; i++)
	    for (j=0; j<=i; j++)
		D[i][j]=(1.0-Strict[i][j])*D[i][j]+Strict[i][j]*Dist[i][j];
    }

    /* free up things */
    free_matrix(D, Rno); free_matrix(M, Rno);
    free(Cdist2);
    
    return(Projno);
}
/* END of grad_proj() */

/* ---- TRIANGLE INEQUALITY BALANCING ---- */

/* dist_metric: given a dist matrix Dist (with squared entries) and a
 vector of squared distances from the centre (Cdist2), then
 from these the metric matrix which is the matrix of scalar products
 of the points in a coord.system centred on the centre of gravity
 is made here. The metric matrix is calc'd using the Cosine Rule.
*/
static void dist_metric(Trimat_ Dist, Trimat_ Metric, double Cdist2[], int Rno)
{
    register int i,j;

    Metric[0][0]=Cdist2[0];
    for (i=1; i<Rno; i++)
    {
	Metric[i][i]=Cdist2[i];
	for (j=0; j<i; j++)
	    Metric[i][j]=(Cdist2[i]+Cdist2[j]-Dist[i][j])/2.0;
    }	
}
/* END of dist_metric */

/* metric_dist: from the metric matrix Metric, a corresponding distance
 matrix Dist is calculated. (The inverse of dist_metric(..) above.)
*/
static void metric_dist(Trimat_ Metric, Trimat_ Dist, int Rno)
{
    register int i,j;

    for (i=0; i<Rno; i++)
    {
	Dist[i][i]=0.0;
	for (j=0; j<i; j++)
	    Dist[i][j]=Metric[i][i]+Metric[j][j]
		-2.0*Metric[i][j];
    }
}
/* END of metric_dist */

/* trieq_bal: balances the triangle inequalities in Metric.
 If its diagonal has non-negative values, then it is "shifted": 
 twice the abs value of the worst non-negative diagonal element
 is added to all diag elements. The shift value is returned in
 *Diagshf (for the projection). Then every off-diag entry is
 checked for triangle inequality violations: being a scalar product,
 cos(phi)=Metric[i][j]/sqrt(Metric[i][i]*Metric[j][j]), which should
 fall in the range -1.0 ... +1.0. If not, then Metric[i][j] is scaled
 a bit. If a particular fabs(cos(phi)) was >1.0, then its value-1 is
 added to the violation score (ideally,  this should be 0.0).
 The total violation score,  normalised by the metric size, is
 returned.
*/
static float trieq_bal(Trimat_ Metric, int Rno, float *Diagshf)
{
    const double ADJFACTOR=0.95;
    register int i,j;
    int Diagcorr;
    register double Shift,Sqroots;
    register float Violsco;

    /* get minimal metric diag (or 0.0 if all non-neg) */
    Shift=0.0; Diagcorr=0;
    for (i=0; i<Rno; i++)
	if (Metric[i][i]<Shift)
	{
	    Shift=Metric[i][i];
	    Diagcorr++;	/* correction needed */
	}

    /* shift centredists if necessary */
    if (Diagcorr)
    {
	Shift*=2.0;
	for (i=0; i<Rno; i++)
	    Metric[i][i]-=Shift;
    }

    /* generate, check and adjust off-diagonals */
    Violsco=0.0;
    for (i=1; i<Rno; i++)
	for (j=0; j<i; j++)	/* Cosine rule */
	{
	    Sqroots=sqrt(Metric[i][i]*Metric[j][j]);
	    if (Metric[i][j]<-Sqroots)
	    {
		Violsco+=(-Metric[i][j]/Sqroots);
		Metric[i][j]=-ADJFACTOR*Sqroots;
	    }
	    else if (Metric[i][j]>Sqroots)
	    {
		Violsco+=(Metric[i][j]/Sqroots);
		Metric[i][j]=ADJFACTOR*Sqroots;
	    }
	}
    *Diagshf=(float)Shift;
    return(Violsco/(Rno*(Rno-1)/2));
}
/* END of trieq_bal */

/* ---- PROJECTION ---- */

/* metric_project: the Crippen/Havel distance matrix projector. Takes
 a metric matrix Metric (lower triangle) of size Rno*Rno.If it is diagonal
 then the vectors to the Rno points form a linearly independent base
 and therefore can live in a Rno-dimensional space only.
 The points can be forced into a lower-dimensional subspace,
 however, by diagonalising the metric matrix
 and then calculating the coordinates using the largest eigenvalues
 only. These coord's are stored in Xyz (a square matrix).
 If there was a shift in trieq_bal(..), then that
 value (Diagshf) is subtracted from all eigenvalues.
 The new subspace cannot have more than Oldim dimensions;
 also, it cannot be less than 3. If all eigenvalues are non-negative, 
 then an Equ fraction of them is kept. 
 Return value: the no. of dimensions of the subspace into which the
 points were projected. The sum of the first abs eigenvalues accounting
 for the new subspace divided by the sum of all abs eigenvalues are
 returned in Projqual: the closer to 1.0,  the better.
*/
static int metric_project(Trimat_ Metric, int Rno, Sqmat_ Xyz,
	float Diagshf, float Equ, int Oldim, float *Projqual)
{
    Sqmat_ Evec; 
    double *Eval, *Sqeval, *Sumeval;
    float Qu;
    double Abseval, Abssum, Newsum;
    register int i,j,Dim;
    
    /* build arrays */
    /* eigenvalues: real, and square roots */
    Eval=(double *) calloc(Rno,sizeof(double));
    Sqeval=(double *) calloc(Rno,sizeof(double));
    Sumeval=(double *) calloc(Rno,sizeof(double));
    /* eigenvectors: Rno*Rno */
    Evec=alloc_sqmat(Rno);
    
    /* get eigenvalues and eigenvectors of the metric matrix:
     in general, use the fast and precise algorithm of
     Housholder tridiag+QL transform: Numerical Recipes;
     Eigenvectors are in the rows of Evec, sorted */
    eigen_ql(Metric,Rno,Eval,Evec); /* Housholder + QL */

    /* subtract the Diagshf value from the eigenvalues */
    for (i=0; i<Oldim; i++)
	Eval[i]-=Diagshf;

    /* find the new dimension */
    if (Eval[Oldim-1]<=0.0)
    {
	/* go up to first positive */
	for (Dim=Oldim-1; Dim>0 && Eval[Dim-1]<=0.0; Dim--);
	if (Dim<1) Dim=1;
    }
    else
    {
	/* make up sums of eigenvalues */
	Sumeval[0]=Eval[0];
	for (i=1; i<Oldim; i++)
	    Sumeval[i]=Eval[i]+Sumeval[i-1];

	/* determine no. of dimensions for Equ-th of total */
	Qu=Sumeval[Oldim-1]*Equ;
	for (Dim=1; Dim<=Oldim && Sumeval[Dim-1]<Qu; Dim++);
    }

    if (Dim>=Oldim) Dim=Oldim-1;	/* force shrink */
    if (Dim<3 && Eval[2]>=0.0) Dim=3;	/* no planars */

    /* make up Cartesian coordinates in a Dim-dimensional subspace */
    for (j=0; j<Dim; j++)
	Sqeval[j]=sqrt(Eval[j]);	/* sq.root of largest evals */
    for (i=0; i<Rno; i++)
	for (j=0; j<Dim; j++)
	    Xyz[i][j]=Sqeval[j]*Evec[j][i];	/* row eigenvectors */

    /* calculate projection quality */
    Abssum=Newsum=0.0;
    for (i=0; i<Dim; i++)
    {
	Abseval=fabs(Eval[i]);
	Abssum+=Abseval; Newsum+=Abseval;
    }
    for (i=Dim; i<Rno; i++)
	Abssum+=fabs(Eval[i]);
	
    /* return new dimension and projection quality */
    *Projqual=(float)(Newsum/Abssum);

    /* cleanup */
    free(Eval); free(Sqeval);
    free(Sumeval);
    free_matrix(Evec, Rno);

    return(Dim);
}
/* END of metric_project */

/* centre_dist: given the distance matrix of Rno points (lower triangle
 is sufficient), the distances of the points from their common centre
 of gravity is returned in Cdist2[]. Based on Lagrange's Theorem. 
 Dist should contain the square of distances. Returns the number of
 Cdist2^-s which happen to be <0.0
*/
static int centre_dist(Trimat_ Dist, int Rno, double Cdist2[])
{
    register double Trisum,Isum;
    register int i,j,k;
    int Cderr;

    /* get lower triangle squared sum */
    Trisum=0.0;
    for (j=0; j<Rno; j++)
	for (k=0; k<j; k++)
	    Trisum+=Dist[j][k];
    Trisum/=(Rno*Rno);

    /* get squared distance sums for the i-th point */
    Cderr=0;
    for (i=0; i<Rno; i++)
    {
	Isum=0.0;
	for (j=0; j<i; j++)
	    Isum+=Dist[i][j];
	for (j=i+1; j<Rno; j++)
	    Isum+=Dist[j][i];
	Cdist2[i]=Isum/Rno-Trisum;
	if (Cdist2[i]<0.0)
	{
	    fprintf(stderr,"! centre_dist: Cdist2[%d]=%e\n",i,Cdist2[i]);
	    Cderr++;
	}
    }
    return(Cderr);
}
/* END of centre_dist */

/* new_distmat: reconstructs the dist matrix (Newdist) from the
 projected coordinates in the matrix Xyz (Rno rows, for each point,
 and Dim columns, for coordinates of points). Newdist is a
 lower triangle distance matrix (Rno*Rno), squared.
*/
static void new_distmat(Sqmat_ Xyz, int Rno, int Dim, Trimat_ Newdist)
{
    register double Dx,Nd;
    register int i,j,k;

    for (i=0; i<Rno; i++)
    {
	Newdist[i][i]=0.0;
	for (j=0; j<i; j++)
	{
	    Nd=0.0;
	    for (k=0; k<Dim; k++)
	    {
		Dx=Xyz[i][k]-Xyz[j][k];
		Nd+=(Dx*Dx);
	    }
	    Newdist[i][j]=Nd;
	}
    }
}
/* END of new_distmat */

#ifdef __SUN_OS__
#undef sqrtf
#undef expf
#endif

/* ==== END OF FUNCTIONS gradproj.c ==== */
