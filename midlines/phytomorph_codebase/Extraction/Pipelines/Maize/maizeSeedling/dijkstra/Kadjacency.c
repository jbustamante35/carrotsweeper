/*

 
 Returns the Sparse adjacency matrix of the K neighbours of xi.

 Inputs:
 -------

   X      - Data (d x N) in double precision
   K      - Number of neighbourg

   
 Outputs:
 -------

   A      - Sparse adjacency matrix (N x N) with KN non-zeros elements

Example
-------


s                 = 1;
e                 = 3;
d                 = 2;
N                 = 1000;
K                 = 5;
X                 = rand(d , N);
A                 = Kadjacency(X , K);
[path , pathcost] = dijkstra(A , s , e);
gplot(A , X');hold on,plot(X(1 , :) , X(2 , :) , 'k+',X(1 , path) , X(2 , path) , 'r' , 'markersize' , 5 , 'linewidth' , 2), hold off


To compile
----------

mex Kadjacency.c

mex -g  Kadjacency.c

mex Kadjacency.c -largeArrayDims


 Author : Sébastien PARIS : sebastien.paris@lsis.org
 -------  Date : 11/01/2007

 Reference ""


*/

#include <math.h>
#include <limits.h>
#include "mex.h"

#define MAX_INF 1.7977e+307

/*-------------------------------------------------------------------------------------*/

void qsindex( double * , int * , int , int  ); 

/*-------------------------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray** plhs ,  int nrhs, const mxArray** prhs)
{
    double *X;
	int K;
	double *A;
	int d , N;
	int i , j, l, k , nzmax;
	double *dist;
	mwIndex *jc , *ir; 
	int id , jd , max_inx , ind;
	double temp, adist, max_dist , val;

	if (nrhs != 2)
	{
		mxErrMsgTxt("Too few arguments");
	}
    
	/* Input 1 */
    X        = mxGetPr(prhs[0]);
    d        = mxGetM(prhs[0]);
    N        = mxGetN(prhs[0]);
	
	/* Input 2 */
	
    K        = (int) mxGetScalar(prhs[1]);
    if (K < 0 || K > N)
	{
        mxErrMsgTxt("K must be < N");
	}
    dist      = (double *)malloc(N*sizeof(double));
	nzmax     = K*N;

    plhs[0]   = mxCreateSparse(N , N , nzmax , mxREAL);
    A         = mxGetPr(plhs[0]);  /* (KN x 1) vector */
    ir        = mxGetIr(plhs[0]);  /* (KN x 1) index vector */
	jc        = mxGetJc(plhs[0]);  /* (N+1 x 1) index vector */
	
	/* Main Call */
		
	k         = 0;
	for( j = 0 ; j < N ; j++ ) 
	{
		jd            = j*d;
		jc[j]         = k;
		for( i = 0 ; i < N ; i++ ) 
		{
			id        = i*d;
			adist     = 0.0;
			for( l = 0 ; l < d ; l++ ) 
			{
				temp     = (X[l + id] - X[l + jd]);
				adist   += (temp*temp); 
			}
			if (adist != 0.0)
			{
				dist[i]   = adist;
			}
			else
			{
				dist[i]   = MAX_INF;
			}
		}
		for( l = 0 ; l < K ; l++) 
		{
			max_dist = MAX_INF;
			max_inx  = 0;
			for( i = 0 ; i < N ; i++ ) 
			{
				if( max_dist > dist[i] ) 
				{
					max_inx  = i;
					max_dist = dist[i];
				}
			}
			ir[k]           = max_inx;
			A[k]            = max_dist;
			dist[ max_inx ] = MAX_INF;
			k++;
		}
	}
	jc[N] = k;
	free(dist);
}
/* ------------------------------------------------------------------------------------------------------------------- */
void qsindex (double  *a, int *index , int lo, int hi)
{
/*  lo is the lower index, hi is the upper index
  of the region of array a that is to be sorted
*/
    int i=lo, j=hi , ind;
    double x=a[(lo+hi)/2] , h;

    /*  partition */
    do
    {    
        while (a[i]<x) i++; 
        while (a[j]>x) j--;
        if (i<=j)
        {
            h        = a[i]; 
			a[i]     = a[j]; 
			a[j]     = h;
			ind      = index[i];
			index[i] = index[j];
			index[j] = ind;
            i++; 
			j--;
        }
    }
	while (i<=j);

    /*  recursion */
    if (lo<j) qsindex(a , index , lo , j);
    if (i<hi) qsindex(a , index , i , hi);
}
/* ------------------------------------------------------------------------------------------------------------------- */

