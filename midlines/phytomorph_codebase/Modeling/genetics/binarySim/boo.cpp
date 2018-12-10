#include <string.h> /* needed for memcpy() */
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <math.h>
#include <matrix.h>
#include <mex.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdint.h>
#define min(a,b) ((a)<(b)?(a):(b))




#define VERBOSE 0
/*
 * doubleelement.c - example found in API guide
 *
 * constructs a 2-by-2 matrix with unsigned 16-bit integers, doubles
 * each element, and returns the matrix
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2007 The MathWorks, Inc.
 */

typedef struct {
    UINT8_T *inX;
    UINT8_T *inY;
    
    unsigned long rowSTART;
    unsigned long rowEND;
    unsigned long JUMPx;
    unsigned long JUMPy;
    unsigned long columnN;
    unsigned long JUMPp;
    unsigned long progN;
    unsigned long threadNUM;
    bool gatherOutput;
    UINT32_T *outZ;
    UINT8_T *out;
} pt_info_t;




/* thread function */
void *doThreadSum(void *tinfo){

    /* extract info from the struct */
    pt_info_t pinfo = *((pt_info_t *)tinfo);
    
    UINT8_T *inX = pinfo.inX;
    UINT8_T *inY = pinfo.inY;
    unsigned long rowSTART = pinfo.rowSTART;
    unsigned long rowEND = pinfo.rowEND;
    unsigned long JUMPx = pinfo.JUMPx;
    unsigned long JUMPy = pinfo.JUMPy;
    unsigned long columnN = pinfo.columnN;
    unsigned long JUMPp = pinfo.JUMPp;
    unsigned long progN = pinfo.progN;
    
    unsigned long outPTR;
    unsigned long threadNUM = pinfo.threadNUM;
    UINT8_T *out = pinfo.out;
    UINT32_T *outZ = pinfo.outZ;
    bool gatherOutput = pinfo.gatherOutput;
    unsigned long xPTR,yPTR;
    
    
    unsigned long row, col, A,prog;
   
    uint32_t acc;
    for (prog=0;prog<progN;prog++) {
        if (VERBOSE==1) {
            mexPrintf("*************************\n");
            mexPrintf("PROG:%d\n",prog);
            mexPrintf("*************************\n");
        }
        for (row=rowSTART;row<rowEND;row++) {
            xPTR = row;
            yPTR = prog;
            if (VERBOSE==1) { mexPrintf("-----------------------\n"); }
            outPTR = JUMPp*prog + row;
            acc = 0;
            for (col=0;col<columnN;col++) {
                A = inX[xPTR] & inY[yPTR];
                if (VERBOSE==1) { mexPrintf("COLUMN:%d:%d:%d->%d\n",col,inX[xPTR],inY[yPTR], A);}
                xPTR = xPTR + JUMPx;
                yPTR = yPTR + JUMPy;
                //out[outPTR] = A;

                //colXptr = colXptr + INdimsX[0];
                //colYptr = colYptr + INdimsY[0];

                if (gatherOutput) {
                    acc = acc + A;
                }
                
                if (A > 0) {
                    out[outPTR] = 1;
                    if (!gatherOutput) {
                        break;
                    }
                }
            }
            if (gatherOutput) {
                //mexPrintf("%d",acc);
                outZ[outPTR] = acc;
            }
            if (VERBOSE==1) {mexPrintf("-----------------------\n");}
            //Aptr = Aptr + 1;
            /* increment row counter */
            //rowXptr = rowXptr + 1;
          }
    }
    
    //mexPrintf("Hello World--%d--[%d].\n",loopVALUE,in[inP]);
    //loopVALUE = 1;
    //printf("start:stop:pnum:%llu:%llu:%llu\n",str,stp,pnum);
    //for (unsigned long i = rowSTART;i < rowEND;i++){
    //    mexPrintf("Hello World--%d--[%d].\n",i,inX[xPTR]);
    //    xPTR = xPTR + JUMPx;
    //}
    
    //mexPrintf("DONE:%d\n",threadNUM);
    /* exit the thread */
    if (VERBOSE==0) {
        pthread_exit(NULL); 
    }
}


#define NDIMS 2
#define TOTAL_ELEMENTS 4


UINT8_T *
myBITOP(const mxArray *array_ptrX,const mxArray *array_ptrY)
{
    
  UINT8_T *prX, *prY; 
  mwSize total_num_of_elements, index; 
  UINT8_T *dynamicData;
  
  total_num_of_elements  = mxGetNumberOfElements(array_ptrY);
  dynamicData = (UINT8_T *)mxCalloc(total_num_of_elements, sizeof(UINT8_T));
  /* Create a local array and load data */
  prX = (UINT8_T *)mxGetData(array_ptrX);
  prY = (UINT8_T *)mxGetData(array_ptrY);
  
 
  for (index=0; index<total_num_of_elements; index++)  {
    dynamicData[index] = *prX++ & *prY++;
  } 
  return dynamicData;
}




/* the gataway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  UINT8_T *dynamicData, *progY;
  const mwSize *INdimsX;
  mwSize A;
  const mwSize *INdimsY;
  int numdims;
  UINT8_T *rowXptr, *colXptr, *rowYptr, *colYptr;
  mwSize Aptr;
  double *inX;
  pthread_t *threads;
  UINT32_T *outZ;
  mwSize total_num_of_elements, index, row, col;
          
  unsigned long rowSTART,rowEND, prog,block;
  mwSize ROW_BLOCK_SIZE,BLOCK_NUMBER,ROW_BLOCK_REM;
  
  
  
  
  //pt_info_t *pinfo;
  
  /* figure out dimensions */
  INdimsX = mxGetDimensions(prhs[0]); 
  /* mexPrintf("\n\tX size is %i x %i ",INdimsX[0],INdimsX[1]); */
  /* figure out dimensions */
  INdimsY = mxGetDimensions(prhs[1]);
  /* mexPrintf("\n\tY size is %i x %i ",INdimsY[0],INdimsY[1]); */
  
  
  (void) nrhs;  /* unused parameters */
  
  /* create a 2-by-2 array of unsigned 16-bit integers */
  plhs[0] = mxCreateNumericArray(0,0,mxUINT8_CLASS,mxREAL);

  rowXptr = (UINT8_T *)mxGetData(prhs[0]);
  rowYptr = (UINT8_T *)mxGetData(prhs[1]);
  progY = (UINT8_T *)mxCalloc(INdimsX[0]*INdimsY[0], sizeof(UINT8_T));
  
  if (nlhs == 2) {
    outZ = (UINT32_T *)mxCalloc(INdimsX[0]*INdimsY[0], sizeof(UINT32_T));
    plhs[1] = mxCreateNumericArray(0,0,mxUINT32_CLASS,mxREAL);
  } else {
    outZ = NULL;
  }
  
  
  BLOCK_NUMBER = min(20,INdimsX[0]);
  ROW_BLOCK_SIZE = floor(INdimsX[0] / BLOCK_NUMBER);
  ROW_BLOCK_REM = remainder(INdimsX[0],BLOCK_NUMBER);
  
  //if (ROW_BLOCK_REM != 0) {
  //    BLOCK_NUMBER = BLOCK_NUMBER + 1;
  //}
  if (VERBOSE==1) {
    mexPrintf("ROWS:%d-->blockN:%d x blockS:%d + blockR:%d\n",INdimsX[0],BLOCK_NUMBER,ROW_BLOCK_SIZE,ROW_BLOCK_REM);
  }
  
  
  /* ************************************ */
  
  threads = new pthread_t[(int)BLOCK_NUMBER];
  
  pt_info_t pinfo[(int)BLOCK_NUMBER];
  
  
  Aptr = 0;
  //for (prog=0;prog<INdimsY[0];prog++) {
      
       for (block=0;block<BLOCK_NUMBER;block++) {
            if (VERBOSE==1) { mexPrintf("BLOCK START...\n"); }
            rowSTART = block*ROW_BLOCK_SIZE;
            rowEND = rowSTART + ROW_BLOCK_SIZE;
            
            
            if (block == (BLOCK_NUMBER-1)) {
               rowEND = rowEND + ROW_BLOCK_REM;
            }
            
           
            pinfo[block].out = progY;
            pinfo[block].inX = rowXptr;
            pinfo[block].inY = rowYptr;
            pinfo[block].rowSTART = rowSTART;
            pinfo[block].rowEND = rowEND;
            pinfo[block].columnN = (unsigned long) INdimsX[1];
            pinfo[block].JUMPx = (unsigned long) INdimsX[0];
            pinfo[block].JUMPy = (unsigned long) INdimsY[0];
            pinfo[block].JUMPp = (unsigned long) INdimsX[0];
            pinfo[block].progN = INdimsY[0];
            pinfo[block].gatherOutput = (nlhs == 2);
            pinfo[block].outZ = outZ;
            pinfo[block].threadNUM = block;

            if (VERBOSE==1) { mexPrintf("--->%d\n",rowXptr[block]);}

            if (VERBOSE==1) {
                doThreadSum(&pinfo[block]);
            }
            
            
            if (VERBOSE==1) {
                mexPrintf("******\n");
                mexPrintf("BLOCK END...\n");
            }
            
            if (VERBOSE==0) {
                pthread_create(&threads[block], NULL, doThreadSum, (void *)(&pinfo[block]));
            }
            //usleep(1000000);
       }
       
       if (VERBOSE==0) {
           for (block=0;block<BLOCK_NUMBER;block++) {
                pthread_join(threads[block], NULL);
           }
       }
  //}
  
  /* ************************************ */
  
  
  
//   Aptr = 0;
//   for (prog=0;prog<INdimsY[0];prog++) {
//       
//       
//       
//       for (row=0;row<INdimsX[0];row++) {
//         colYptr = rowYptr;
//         colXptr = rowXptr+row;
//         for (col=0;col<INdimsX[1];col++) {
//             A = *colXptr & *colYptr;
//             /* mexPrintf("%d:%d->%d\n",*colXptr,*colYptr, A); */
//             
//             colXptr = colXptr + INdimsX[0];
//             colYptr = colYptr + INdimsY[0];
//             
//             if (A > 0) {
//                 progY[Aptr] = 1;
//                 break;
//             }
//         }
//         Aptr = Aptr + 1;
//         /* increment row counter */
//         //rowXptr = rowXptr + 1;
//       }
//       
//       /* increment program counter */
//       rowYptr = rowYptr + 1;
//   }
  
  /* ************************************ */
  
  /* dynamicData = myBITOP(prhs[0],prhs[1]); */

  /* Point mxArray to dynamicData */
  /* mxSetData(plhs[1],dynamicData); */
  /* mxSetM(plhs[0], INdimsX[0]); */
  /* mxSetN(plhs[0], INdimsX[1]); */
  
  mxSetData(plhs[0],progY);
  mxSetM(plhs[0], INdimsX[0]);
  mxSetN(plhs[0], INdimsY[0]);
  
  if (nlhs==2) {
    mxSetData(plhs[1],outZ);
    mxSetM(plhs[1], INdimsX[0]);
    mxSetN(plhs[1], INdimsY[0]);
  }
  
  
}
