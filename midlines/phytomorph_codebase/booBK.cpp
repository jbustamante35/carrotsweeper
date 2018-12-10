#include <string.h> /* needed for memcpy() */
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
/*
 * doubleelement.c - example found in API guide
 *
 * constructs a 2-by-2 matrix with unsigned 16-bit integers, doubles
 * each element, and returns the matrix
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2007 The MathWorks, Inc.
 */


#define NDIMS 2
#define TOTAL_ELEMENTS 4

/* the computational subroutine */
void dbl_elem(unsigned short *x)
{
  unsigned short scalar=2;
  int i,j;

  for(i=0;i<2;i++) {
    for(j=0;j<2;j++) {
      *(x+i*2+j) = scalar * *(x+i*2+j);
    }
  }
}

static void
analyze_uint8(const mxArray *array_ptr)
{
  unsigned char *pr, *pi; 
  mwSize total_num_of_elements, index; 
  
  pr = (unsigned char *)mxGetData(array_ptr);
  //pi = (unsigned char *)mxGetImagData(array_ptr);
  total_num_of_elements = mxGetNumberOfElements(array_ptr);
  
  for (index=0; index<total_num_of_elements; index++)  {
    mexPrintf("\t");
    //display_subscript(array_ptr, index);
    if (mxIsComplex(array_ptr)) {
      mexPrintf(" = %u + %ui\n", *pr, *pi++); 
    } else {
      mexPrintf(" = %u\n", *pr++);
    }
  } 
}

static void
myBITOP(const mxArray *array_ptrX,const mxArray *array_ptrY)
{
  unsigned char *prX, *prY,A; 
  mwSize total_num_of_elements, index; 
  
  prX = (unsigned char *)mxGetData(array_ptrX);
  prY = (unsigned char *)mxGetData(array_ptrY);
  
  //pi = (unsigned char *)mxGetImagData(array_ptr);
  total_num_of_elements = mxGetNumberOfElements(array_ptrX);
  
  for (index=0; index<total_num_of_elements; index++)  {
    mexPrintf("\n");
    A = *prX++ & *prY++;
    mexPrintf(" = %u\n", A);
  } 
}


static void
myBITOPP(const mxArray *array_ptrX,const mxArray *array_ptrY)
{
  unsigned char *prX, *prY,A; 
  mwSize total_num_of_elements, index; 
  size_t T;
  mxArray *array_ptrZ;
  const mwSize *dims;
  mwSize ndim;
  mwIndex *dynamicData;
  double  *temp;
  mexPrintf("loop3");
  total_num_of_elements  = mxGetNumberOfElements(array_ptrX);
  dynamicData = (mwIndex *)mxCalloc(total_num_of_elements, sizeof(UINT8_T));
  /* Create a local array and load data */
  
    //for ( index = 0; index < ELEMENTS; index++ ) {
    //    dynamicData[index] = data[index];
    //}
  mexPrintf("loop2");
  ndim = mxGetNumberOfDimensions(array_ptrX);
  dims = mxGetDimensions(array_ptrX);
  mexPrintf("loop1");
  prX = (unsigned char *)mxGetData(array_ptrX);
  prY = (unsigned char *)mxGetData(array_ptrY);
  
 
  //array_ptrZ = mxCreateNumericArray(ndim,dims,mxUINT8_CLASS,mxREAL);
  
  
  
 
  
  //temp = mxGetPr(array_ptrX);
  mexPrintf("loop");
  for (index=0; index<total_num_of_elements; index++)  {
    mexPrintf("\n");
    A = *prX++ & *prY++;
    mexPrintf("loop0");
    //mxSetData(array_ptrZ, A);
    //dynamicData[index] = A;
    //array_ptrZ[index] = A;
    mexPrintf(" = %u\n", A[0]);
  } 
  
  //mxDestroyArray(array_ptrZ);
}




/* the gataway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  const mwSize dims[]={2,2};
  
  const mwSize *INdimsX;
  const mwSize *INdimsY;
  unsigned char *start_of_pr;
  unsigned short data[]={1,2,3,4};
  size_t bytes_to_copy;
  int dimxX, dimyX,dimxY, dimyY, numdims;
  
  double *inX;
  
  inX = mxGetPr(prhs[0]);
  
  
  // figure out dimensions
  INdimsX = mxGetDimensions(prhs[0]);
  dimyX = (int)INdimsX[0]; 
  dimxX = (int)INdimsX[1];
  mexPrintf("\n\tX size is %i x %i ",dimyX,dimxX);
  
  // figure out dimensions
  INdimsY = mxGetDimensions(prhs[1]);
  dimyY = (int)INdimsY[0]; 
  dimxY = (int)INdimsY[1];
  mexPrintf("\n\tY size is %i x %i ",dimyY,dimxY);
  
  analyze_uint8(prhs[0]);
  
  myBITOP(prhs[0],prhs[1]);
  
  (void) nlhs; (void) nrhs;  /* unused parameters */

  /* call the computational subroutine */
  dbl_elem(data);

  /* create a 2-by-2 array of unsigned 16-bit integers */
  plhs[0] = mxCreateNumericArray(1,dims,mxUINT8_CLASS,mxREAL);

  /* populate the real part of the created array */
  start_of_pr = (unsigned char *)mxGetData(plhs[0]);
  bytes_to_copy = TOTAL_ELEMENTS * mxGetElementSize(plhs[0]);
  memcpy(start_of_pr,data,bytes_to_copy);
}
