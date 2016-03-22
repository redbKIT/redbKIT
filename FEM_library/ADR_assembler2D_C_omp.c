#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "blas.h"

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    /* Check for proper number of arguments. */
    if(nrhs!=19) {
        mexErrMsgTxt("19 inputs are required.");
    } else if(nlhs>6) {
        mexErrMsgTxt("Too many output arguments.");
    }

       
    int noe     = mxGetN(prhs[3]);
    int* nln_ptr = (int*)mxGetData(prhs[4]);/*mxGetPr(prhs[1]);*/
    int nln     = mxGetM(prhs[3])-1;
    int nln2    = nln*nln;
    
    /**/
    plhs[0] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL); 
    plhs[2] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(nln*noe,1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(nln*noe,1, mxREAL);
       
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    double* myMcoef    = mxGetPr(plhs[3]);
    double* myRrows    = mxGetPr(plhs[4]);
    double* myRcoef    = mxGetPr(plhs[5]);
    
    
    /*Get pointer to output
    double myArows[nln2*noe];
    double myAcols[nln2*noe];
    double myAcoef[nln2*noe];
    double myMcoef[nln2*noe];
    double myRrows[nln*noe];
    double myRcoef[nln*noe];
    */
    
    int OP[4] = {1, 1, 1, 1};
    double C_t[2];
    double C_d[2][2];
    /*Get pointer to input data*/
    /* int* TC_d = (int*)mxGetData(prhs[1]);mxGetPr(prhs[1]);*/
    /* int* TC_t = (int*)mxGetData(prhs[2]); mxGetPr(prhs[2]);*/
    double* TC_d   = mxGetPr(prhs[1]);
    double* TC_t   = mxGetPr(prhs[2]);

    int k,l;
    for (k = 0; k < 2; k = k + 1 )
    {
        for (l = 0; l < 2; l = l + 1 )
        {
            C_d[k][l] = 0;
        }
        C_t[k] = 0;
    }
    
    if (TC_d[0]==10 && TC_d[1]==10)
    {
        for (l = 0; l < 2; l = l + 1 )
        {
            C_d[l][l] = 1;
        }
    }
    else
    {
        C_d[(int)(TC_d[0])][(int)(TC_d[1])] = 1;
    }
    
    if (TC_t[0]==10)
    {
        for (l = 0; l < 2; l = l + 1 )
        {
            C_t[l] = 1;
        }
    }
    else
    {
        C_t[(int)(TC_t[0])] = 1;
    }
    
    
    /* Local mass matrix (computed only once) with quadrature nodes */
    double LocalMass[nln][nln];
    int q;
    int NumQuadPoints     = mxGetN(prhs[10]);
    
    double* mu   = mxGetPr(prhs[5]);
    double* bx   = mxGetPr(prhs[6]);
    double* by   = mxGetPr(prhs[7]);
    double* si   = mxGetPr(prhs[8]);
    double* f    = mxGetPr(prhs[9]);
    double* w   = mxGetPr(prhs[10]);
    double* dcdx = mxGetPr(prhs[11]);
    double* dcdy = mxGetPr(prhs[12]);
    double* dedx = mxGetPr(prhs[13]);
    double* dedy = mxGetPr(prhs[14]);
    double* phi = mxGetPr(prhs[15]);
    double* dcsiphi = mxGetPr(prhs[16]);
    double* detaphi = mxGetPr(prhs[17]);
    double* detjac = mxGetPr(prhs[18]);

    for (k = 0; k < nln; k = k + 1 )
    {
        for (l = 0; l < nln; l = l + 1 )
        {
            double tmp = 0;
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                tmp = tmp + phi[k+q*nln] * phi[l+q*nln] * w[q];
            }
            LocalMass[k][l] = tmp;
        }
    }

    
    double gradx[nln][NumQuadPoints];
    double grady[nln][NumQuadPoints];
    double* elements   = mxGetPr(prhs[3]);


    /* Assembly: loop over the elements */
    int ie;
    int a, b;
    int iii = 0;
    int ii = 0;
    double aloc = 0;
    double floc = 0;
    
    #pragma omp parallel for shared(dcdx, dedx, dcdy, dedy, mu, bx, by, si, f, detjac, elements, myRrows, myRcoef,myAcols, myArows, myAcoef, myMcoef) private(gradx,grady,ie,ii,iii,a,b,k,l,q,aloc,floc) firstprivate(phi, dcsiphi, detaphi, w, nln2, nln, OP, C_t, C_d, LocalMass)
    
    /*private(phi, gradx, grady, dcsiphi, detaphi, nln, NumQuadPoints, nln2, LocalMass, OP, C_t, C_d, w, k, l, a, b, ii, iii, aloc, floc, q)*/
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        for (k = 0; k < nln; k = k + 1 )
        {
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                gradx[k][q] = dcdx[ie]*dcsiphi[k+q*nln] + dedx[ie]*detaphi[k+q*nln];
                grady[k][q] = dcdy[ie]*dcsiphi[k+q*nln] + dedy[ie]*detaphi[k+q*nln];
            }
        }

        iii = 0;
        ii = 0;
        
        for (a = 0; a < nln; a = a + 1 )
        {
            
            for (b = 0; b < nln; b = b + 1 )
            {

                aloc = 0;
                
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {
                    
                    aloc = aloc + (
                            OP[0] * (   C_t[0] * bx[ie+q*noe] * gradx[b][q] * phi[a+q*nln]
                            +  C_t[1] * by[ie+q*noe] * grady[b][q] * phi[a+q*nln] )
                            + OP[1] * (    C_d[0][0] * mu[ie+q*noe] * gradx[b][q] * gradx[a][q]
                                        +  C_d[0][1] * mu[ie+q*noe] * gradx[b][q] * grady[a][q]
                                        +  C_d[1][0] * mu[ie+q*noe] * grady[b][q] * gradx[a][q]
                                        +  C_d[1][1] * mu[ie+q*noe] * grady[b][q] * grady[a][q] )
                            + OP[2] * (  si[ie+q*noe] * phi[b+q*nln] * phi[a+q*nln] )
                            ) * w[q];             
                }
 
                myArows[ie*nln2+iii] = elements[a+ie*(nln+1)];
                myAcols[ie*nln2+iii] = elements[b+ie*(nln+1)];
                myAcoef[ie*nln2+iii] = aloc*detjac[ie];
                myMcoef[ie*nln2+iii] = LocalMass[a][b]*detjac[ie];
                
                iii = iii + 1;

            }
            
            
            floc = 0;
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                floc = floc + ( OP[3] * phi[a+q*nln] * f[ie+q*noe] ) * w[q];
            }
            myRrows[ie*nln+ii] = elements[a+ie*(nln+1)];
            myRcoef[ie*nln+ii] = floc*detjac[ie];
            /*
             printf("\n ie*nln+ii = %d, myRrows = %f, myRrows2 = %f, from thread %d", ie*nln+ii, myRrows[ie*nln+ii], elements[a+ie*(nln+1)], omp_get_thread_num());
            */
            ii = ii + 1;
        }
        
    }
    
    
    /*
    int jj;
    for (jj = 0; jj < noe*nln; jj = jj + 1 )
    {
         printf("\n jj = %d, Rr = %f from thread %d", jj, myRrows[jj], omp_get_thread_num());
    }
    */
    /*create space for output
    plhs[0] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL); 
    plhs[2] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nln2*noe,1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(nln*noe,1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(nln*noe,1, mxREAL);
    
    double* Arows    = mxGetPr(plhs[0]);
    double* Acols    = mxGetPr(plhs[1]);
    double* Acoef    = mxGetPr(plhs[2]);
    double* Mcoef    = mxGetPr(plhs[3]);
    double* Rrows    = mxGetPr(plhs[4]);
    double* Rcoef    = mxGetPr(plhs[5]);
    Arows = &myArows;
    Acols = &myAcols;
    Acoef = &myAcoef;
    Mcoef = &myMcoef;
    Rrows = &myRrows;
    Rcoef = &myRcoef;
     */
    
            
}

