#include "mex.h"
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "blas.h"
#include <string.h>

void mexFunction(int nlhs,mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    
    /* Check for proper number of arguments. */
    if(nrhs!=26) {
        mexErrMsgTxt("26 inputs are required.");
    } else if(nlhs>6) {
        mexErrMsgTxt("Too many output arguments.");
    }

       
    int noe     = mxGetN(prhs[3]);
    double* nln_ptr = mxGetPr(prhs[4]);
    int nln     = (int)(nln_ptr[0]); 
    int numRowsElements  = mxGetM(prhs[3]); 
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
    
    /* copy the string data from prhs[0] into a C string input_ buf.    */
    char *OP_string = mxArrayToString(prhs[0]);
    int OP[4] = {0, 0, 0, 0};
    if (strcmp(OP_string, "diffusion"))
    {
        OP[0] = 1;
    }
    
    if (strcmp(OP_string, "transport"))
    {
        OP[1] = 1;
    }
    
    if (strcmp(OP_string, "reaction"))
    {
        OP[2] = 1;
    }
    
    if (strcmp(OP_string, "source"))
    {
        OP[3] = 1;
    }
    
    if (strcmp(OP_string, "all"))
    {
        OP[0] = 1;
        OP[1] = 1;
        OP[2] = 1;
        OP[3] = 1;
    }
    
    double C_t[3];
    double C_d[3][3];
    
    double* TC_d   = mxGetPr(prhs[1]);
    double* TC_t   = mxGetPr(prhs[2]);

    int k,l;
    for (k = 0; k < 3; k = k + 1 )
    {
        for (l = 0; l < 3; l = l + 1 )
        {
            C_d[k][l] = 0;
        }
        C_t[k] = 0;
    }
    
    if (TC_d[0]==10 && TC_d[1]==10)
    {
        for (l = 0; l < 3; l = l + 1 )
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
        for (l = 0; l < 3; l = l + 1 )
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
    int NumQuadPoints     = mxGetN(prhs[11]);
    
    double* mu   = mxGetPr(prhs[5]);
    double* bx   = mxGetPr(prhs[6]);
    double* by   = mxGetPr(prhs[7]);
    double* bz   = mxGetPr(prhs[8]);
    double* si   = mxGetPr(prhs[9]);
    double* f    = mxGetPr(prhs[10]);
    double* w   = mxGetPr(prhs[11]);
    double* dcdx = mxGetPr(prhs[12]);
    double* dcdy = mxGetPr(prhs[13]);
    double* dcdz = mxGetPr(prhs[14]);
    double* dedx = mxGetPr(prhs[15]);
    double* dedy = mxGetPr(prhs[16]);
    double* dedz = mxGetPr(prhs[17]);
    double* dtdx = mxGetPr(prhs[18]);
    double* dtdy = mxGetPr(prhs[19]);
    double* dtdz = mxGetPr(prhs[20]);
    double* phi = mxGetPr(prhs[21]);
    double* dcsiphi = mxGetPr(prhs[22]);
    double* detaphi = mxGetPr(prhs[23]);
    double* dtauphi = mxGetPr(prhs[24]);
    double* detjac = mxGetPr(prhs[25]);

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
    double gradz[nln][NumQuadPoints];
    double* elements   = mxGetPr(prhs[3]);

    /* Assembly: loop over the elements */
    int ie;
    int a, b;
    int iii = 0;
    int ii = 0;
    double aloc = 0;
    double floc = 0;
    
    #pragma omp parallel for shared(dcdx,dedx,dtdx,dcdy,dedy,dtdy,dcdz,dedz,dtdz,mu,bx,by,bz,si,f,detjac,elements,myRrows, myRcoef,myAcols, myArows, myAcoef, myMcoef) private(gradx,grady,gradz,ie,ii,iii,a,b,k,l,q,aloc,floc) firstprivate(phi,dcsiphi,detaphi,dtauphi,w,numRowsElements, nln2, nln, OP, C_t, C_d, LocalMass)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        for (k = 0; k < nln; k = k + 1 )
        {
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                gradx[k][q] = dcdx[ie]*dcsiphi[k+q*nln] + dedx[ie]*detaphi[k+q*nln] + dtdx[ie]*dtauphi[k+q*nln];
                grady[k][q] = dcdy[ie]*dcsiphi[k+q*nln] + dedy[ie]*detaphi[k+q*nln] + dtdy[ie]*dtauphi[k+q*nln];
                gradz[k][q] = dcdz[ie]*dcsiphi[k+q*nln] + dedz[ie]*detaphi[k+q*nln] + dtdz[ie]*dtauphi[k+q*nln];
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
                                     +  C_t[1] * by[ie+q*noe] * grady[b][q] * phi[a+q*nln]
                                     +  C_t[2] * bz[ie+q*noe] * grady[b][q] * phi[a+q*nln]
                                    )
                            + OP[1] * (    C_d[0][0] * mu[ie+q*noe] * gradx[b][q] * gradx[a][q]
                                        +  C_d[0][1] * mu[ie+q*noe] * gradx[b][q] * grady[a][q]
                                        +  C_d[0][2] * mu[ie+q*noe] * gradx[b][q] * gradz[a][q]
                                        +  C_d[1][0] * mu[ie+q*noe] * grady[b][q] * gradx[a][q]
                                        +  C_d[1][1] * mu[ie+q*noe] * grady[b][q] * grady[a][q] 
                                        +  C_d[1][2] * mu[ie+q*noe] * grady[b][q] * gradz[a][q]
                                        +  C_d[2][0] * mu[ie+q*noe] * gradz[b][q] * gradx[a][q]
                                        +  C_d[2][1] * mu[ie+q*noe] * gradz[b][q] * grady[a][q] 
                                        +  C_d[2][2] * mu[ie+q*noe] * gradz[b][q] * gradz[a][q]
                                      )
                            + OP[2] * (  si[ie+q*noe] * phi[b+q*nln] * phi[a+q*nln] )
                            ) * w[q];             
                }
 
                myArows[ie*nln2+iii] = elements[a+ie*numRowsElements];
                myAcols[ie*nln2+iii] = elements[b+ie*numRowsElements];
                myAcoef[ie*nln2+iii] = aloc*detjac[ie];
                myMcoef[ie*nln2+iii] = LocalMass[a][b]*detjac[ie];
                
                iii = iii + 1;

            }
            
            
            floc = 0;
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                floc = floc + ( OP[3] * phi[a+q*nln] * f[ie+q*noe] ) * w[q];
            }
            myRrows[ie*nln+ii] = elements[a+ie*numRowsElements];
            myRcoef[ie*nln+ii] = floc*detjac[ie];
    
            ii = ii + 1;
        }
        
    }
            
}

