/*
 
 Returns the Sparse adjacency matrix of neighbours of xi within the ball radius R, such {x_j}, j<>i, d(xi,xj)<R�

 Inputs:
 -------

   X        Data (d x N) in double precision
   R        Radius

 
 Outputs:
 -------

  A        Sparse adjacency matrix (N x N) with KN non-zeros elements


Example
-------

s                 = 1;
e                 = 30;
d                 = 2;
N                 = 1000;
R                 = 2/sqrt(N);
X                 = rand(d , N);
A                 = Radjacency(X , R);
[path , pathcost] = dijkstra(A , s , e);
gplot(A , X');hold on,plot(X(1 , :) , X(2 , :) , 'k+',X(1 , path) , X(2 , path) , 'r' , 'markersize' , 5 , 'linewidth' , 2), hold off


To compile
----------

mex Radjacency.c

mex Radjacency.c -largeArrayDims

mex -g -output Radjacency.dll Radjacency.c

mex -f mexopts_intel10.bat -output Radjacency.dll Radjacency.c


 Author : S�bastien PARIS : sebastien.paris@lsis.org
 -------  Date : 11/01/2007

 Reference ""


*/

#include <math.h>
#include <limits.h>
#include "mex.h"

#define MAX_INF 1.7977e+307

/*-------------------------------------------------------------------------------------*/

void mexFunction(int nlhs, mxArray** plhs ,  int nrhs, const mxArray** prhs)
{
    double *X;
	double R;
	double *A;
	int d , N;
	int i , j, l, k , oldnzmax , nzmax ;	
	mwIndex *jc , *ir;
	int id , jd ;
	double temp, adist , percent_sparse;
		
	if (nrhs != 2)
	{    
        /*
		mxErrMsgTxt("Too few arguments");	
         */
        mexErrMsgTxt("Too few arguments");	
	}
    
	/* Input 1 */
	
    X        = mxGetPr(prhs[0]);
    d        = mxGetM(prhs[0]);
    N        = mxGetN(prhs[0]);
	
	/* Input 2 */
	
    R              = mxGetScalar(prhs[1]);
    if (R < 0.0)
	{
        /*
        mxErrMsgTxt("R must be > 0");
         */
        mexErrMsgTxt("R must be > 0");
	}
	R              *= R;
		
    percent_sparse = 0.2;
    nzmax          = (int)ceil(N*N*percent_sparse);
	
    plhs[0]        = mxCreateSparse(N , N , nzmax , mxREAL);
    A              = mxGetPr(plhs[0]);  /* (nzmax x 1) vector */
    ir             = mxGetIr(plhs[0]);  /* (nzmax x 1) index vector */
	jc             = mxGetJc(plhs[0]);  /* (N+1 x 1) index vector */

	/* ----- Main Call ------*/

	k              = 0;
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
			if ((adist != 0.0) && (adist < R))
			{
				if (k >= nzmax)
				{
					oldnzmax        = nzmax;
					percent_sparse += 0.1;
					nzmax           = (int)ceil(N*N*percent_sparse);
					if (oldnzmax == nzmax)
					{
						nzmax++;
					}
					
					mxSetNzmax(plhs[0], nzmax); 
					mxSetPr(plhs[0], mxRealloc(A , nzmax*sizeof(double)));
					mxSetIr(plhs[0], mxRealloc(ir , nzmax*sizeof(mwIndex)));
					A                = mxGetPr(plhs[0]);
					ir               = mxGetIr(plhs[0]);
				}
				A[k]  = adist;
				ir[k] = i;
				k++;
			}
		}
	}
	jc[N] = k;
}
/*-------------------------------------------------------------------------------------*/

