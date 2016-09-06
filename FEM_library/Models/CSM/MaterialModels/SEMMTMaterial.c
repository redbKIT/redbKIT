/*   This file is part of redbKIT.
 *   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
 *   Author: Federico Negri <federico.negri@epfl.ch>
 */

#include "SEMMTMaterial.h"

/*************************************************************************/
void SEMMTMaterial_forces(mxArray* plhs[], const mxArray* prhs[])
{
    
    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;
    
    plhs[0] = mxCreateDoubleMatrix(nln*noe*dim,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln*noe*dim,1, mxREAL);
    
    double* myRrows    = mxGetPr(plhs[0]);
    double* myRcoef    = mxGetPr(plhs[1]);
    
    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    int NumNodes          = (int)(mxGetM(prhs[3]) / dim);
    
    double* U_h   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phi = mxGetPr(prhs[9]);
    double* gradrefphi = mxGetPr(prhs[10]);
    
    double gradphi[dim][nln][NumQuadPoints];
    double* elements  = mxGetPr(prhs[4]);
    
    double GradV[dim][dim];
    double GradU[dim][dim];
    double GradUh[dim][dim][NumQuadPoints];
    
    double Id[dim][dim];
    int d1,d2;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            Id[d1][d2] = 0;
            if (d1==d2)
            {
                Id[d1][d2] = 1;
            }
        }
    }
    
    double F[dim][dim];
    double EPS[dim][dim];
    double dP[dim][dim];
    double P_Uh[dim][dim];
    
    double* material_param = mxGetPr(prhs[2]);
    double Young = material_param[0];
    double Poisson = material_param[1];
    double Stiffening_power = material_param[2];
    double mu = Young / (2 + 2 * Poisson);
    double lambda =  Young * Poisson /( (1+Poisson) * (1-2*Poisson) );
    
    /* Assembly: loop over the elements */
    int ie;
    
#pragma omp parallel for shared(invjac,detjac,elements,myRrows,myRcoef,U_h) private(gradphi,F,EPS,dP,P_Uh,GradV,GradU,GradUh,ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,Id,mu,lambda)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        for (k = 0; k < nln; k = k + 1 )
        {
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphi[d1][k][q] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphi[d1][k][q] = gradphi[d1][k][q] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                    }
                }
            }
        }
        
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                {
                    GradUh[d1][d2][q] = 0;
                    for (k = 0; k < nln; k = k + 1 )
                    {
                        int e_k;
                        e_k = (int)(elements[ie*numRowsElements + k] + d1*NumNodes - 1);
                        GradUh[d1][d2][q] = GradUh[d1][d2][q] + U_h[e_k] * gradphi[d2][k][q];
                    }
                }
            }
        }
        
        int iii = 0;
        int ii = 0;
        int a, b, i_c, j_c;
        
        /* loop over test functions --> a */
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < dim; i_c = i_c + 1 )
            {
                /* set gradV to zero*/
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        GradV[d1][d2] = 0;
                    }
                }
                                
                double rloc = 0;
                for (q = 0; q < NumQuadPoints; q = q + 1 )
                {
                    
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        GradV[i_c][d2] = gradphi[d2][a][q];
                    }
                    
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                        {
                            F[d1][d2] = Id[d1][d2] + GradUh[d1][d2][q];
                        }
                    }
                    
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                        {
                            EPS[d1][d2] = 0.5 * ( F[d1][d2] + F[d2][d1] ) - Id[d1][d2];
                        }
                    }
                    
                    double trace = Trace(dim, EPS);
                    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                    {
                        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                        {
                            P_Uh[d1][d2] = 2 * mu * EPS[d1][d2] + lambda * trace * Id[d1][d2];
                        }
                    }
                    rloc  = rloc + Mdot( dim, GradV, P_Uh) * w[q];
                }
                                            
                myRrows[ie*nln*dim+ii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                myRcoef[ie*nln*dim+ii] = rloc * detjac[ie] * pow( detjac[0] / detjac[ie], Stiffening_power );
                ii = ii + 1;
            }
        }
    }
        
}


/*************************************************************************/
void SEMMTMaterial_jacobianFast3D(mxArray* plhs[], const mxArray* prhs[])
{
    
    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;
    
    plhs[0] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    
    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    int NumNodes          = (int)(mxGetM(prhs[3]) / dim);
    
    double* U_h   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phi = mxGetPr(prhs[9]);
    double* gradrefphi = mxGetPr(prhs[10]);
    
    double* elements  = mxGetPr(prhs[4]);
    
    int d1,d2;
    
    double* material_param = mxGetPr(prhs[2]);
    double Young = material_param[0];
    double Poisson = material_param[1];
    double Stiffening_power = material_param[2];

    double mu = Young / (2 + 2 * Poisson);
    double lambda =  Young * Poisson /( (1+Poisson) * (1-2*Poisson) );
    
    double detjac_ref = detjac[0];
    
    /* Assembly: loop over the elements */
    int ie;
    
#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef,U_h) private(ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,mu,lambda,detjac_ref,Stiffening_power)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        double gradphi[NumQuadPoints][dim][nln];
        
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* Compute Gradient of Basis functions*/
            for (k = 0; k < nln; k = k + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphi[q][d1][k] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphi[q][d1][k] = gradphi[q][d1][k] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                    }
                }
            }
        }
        
        int iii = 0;
        int a, b, i_c, j_c;
        
        double aloc[nln][dim][nln][dim];
        /* loop over test functions --> a */
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < 3; i_c = i_c + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < 3; j_c = j_c + 1 )
                    {
                        aloc[a][i_c][b][j_c]  = 0.0;
                    }
                }
            }
        }
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* loop over test functions --> a */
            for (a = 0; a < nln; a = a + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    
                    aloc[a][0][b][0] += ( gradphi[q][0][a]*(lambda*gradphi[q][0][b] + 2.0*mu*gradphi[q][0][b]) + mu*gradphi[q][1][b]*gradphi[q][1][a] + mu*gradphi[q][2][b]*gradphi[q][2][a] ) * w[q];
                    
                    aloc[a][0][b][1] += ( lambda*gradphi[q][1][b]*gradphi[q][0][a] + mu*gradphi[q][0][b]*gradphi[q][1][a] ) * w[q];
                    
                    aloc[a][0][b][2] += ( lambda*gradphi[q][2][b]*gradphi[q][0][a] + mu*gradphi[q][0][b]*gradphi[q][2][a] ) * w[q];
                    
                    aloc[a][1][b][0] += ( lambda*gradphi[q][0][b]*gradphi[q][1][a] + mu*gradphi[q][1][b]*gradphi[q][0][a] ) * w[q];
                    
                    aloc[a][1][b][1] += ( gradphi[q][1][a]*(lambda*gradphi[q][1][b] + 2.0*mu*gradphi[q][1][b]) + mu*gradphi[q][0][b]*gradphi[q][0][a] + mu*gradphi[q][2][b]*gradphi[q][2][a] ) * w[q];
                    
                    aloc[a][1][b][2] += ( lambda*gradphi[q][2][b]*gradphi[q][1][a] + mu*gradphi[q][1][b]*gradphi[q][2][a] ) * w[q];
                    
                    aloc[a][2][b][0] += ( lambda*gradphi[q][0][b]*gradphi[q][2][a] + mu*gradphi[q][2][b]*gradphi[q][0][a] ) * w[q];
                    
                    aloc[a][2][b][1] += ( lambda*gradphi[q][1][b]*gradphi[q][2][a] + mu*gradphi[q][2][b]*gradphi[q][1][a] ) * w[q];
                    
                    aloc[a][2][b][2] += ( gradphi[q][2][a]*(lambda*gradphi[q][2][b] + 2.0*mu*gradphi[q][2][b]) + mu*gradphi[q][0][b]*gradphi[q][0][a] + mu*gradphi[q][1][b]*gradphi[q][1][a] ) * w[q];

                }
            }
        }

        double modified_detJ = detjac[ie] * pow( detjac_ref / detjac[ie], Stiffening_power );
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < 3; i_c = i_c + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < 3; j_c = j_c + 1 )
                    {
                        myArows[ie*nln2*9+iii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                        myAcols[ie*nln2*9+iii] = elements[b+ie*numRowsElements] + j_c * NumNodes;
                        myAcoef[ie*nln2*9+iii] = aloc[a][i_c][b][j_c]*modified_detJ;
                        iii = iii + 1;
                    }
                }
            }
        }
    }
    
}
/*************************************************************************/
void SEMMTMaterial_jacobianFast2D(mxArray* plhs[], const mxArray* prhs[])
{
    
    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;
    
    plhs[0] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    
    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    int NumNodes          = (int)(mxGetM(prhs[3]) / dim);
    
    double* U_h   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phi = mxGetPr(prhs[9]);
    double* gradrefphi = mxGetPr(prhs[10]);
    
    double* elements  = mxGetPr(prhs[4]);
    
    int d1,d2;
    
    double* material_param = mxGetPr(prhs[2]);
    double Young = material_param[0];
    double Poisson = material_param[1];
    double Stiffening_power = material_param[2];
    double mu = Young / (2 + 2 * Poisson);
    double lambda =  Young * Poisson /( (1+Poisson) * (1-2*Poisson) );
    
    double detjac_ref = detjac[0];

    /* Assembly: loop over the elements */
    int ie;
    
#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef,U_h) private(ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,mu,lambda,detjac_ref,Stiffening_power)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        double gradphi[NumQuadPoints][dim][nln];
        
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* Compute Gradient of Basis functions*/
            for (k = 0; k < nln; k = k + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphi[q][d1][k] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphi[q][d1][k] = gradphi[q][d1][k] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                    }
                }
            }
        }
        
        int iii = 0;
        int a, b, i_c, j_c;
        
        double aloc[nln][dim][nln][dim];
        /* loop over test functions --> a */
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < 2; i_c = i_c + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < 2; j_c = j_c + 1 )
                    {
                        aloc[a][i_c][b][j_c]  = 0.0;
                    }
                }
            }
        }
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            /* loop over test functions --> a */
            for (a = 0; a < nln; a = a + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    
                    aloc[a][0][b][0] += ( gradphi[q][0][a]*(lambda*gradphi[q][0][b] + 2.0*mu*gradphi[q][0][b]) + mu*gradphi[q][1][b]*gradphi[q][1][a] ) * w[q];
                    
                    aloc[a][0][b][1] += ( lambda*gradphi[q][1][b]*gradphi[q][0][a] + mu*gradphi[q][0][b]*gradphi[q][1][a] ) * w[q];
                    
                    aloc[a][1][b][0] += ( lambda*gradphi[q][0][b]*gradphi[q][1][a] + mu*gradphi[q][1][b]*gradphi[q][0][a] ) * w[q];
                    
                    aloc[a][1][b][1] += ( gradphi[q][1][a]*(lambda*gradphi[q][1][b] + 2.0*mu*gradphi[q][1][b]) + mu*gradphi[q][0][b]*gradphi[q][0][a] ) * w[q];
                }
            }
        }
        
        double modified_detJ = detjac[ie] * pow( detjac_ref / detjac[ie], Stiffening_power );
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < 2; i_c = i_c + 1 )
            {
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < 2; j_c = j_c + 1 )
                    {
                        myArows[ie*nln2*4+iii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                        myAcols[ie*nln2*4+iii] = elements[b+ie*numRowsElements] + j_c * NumNodes;
                        myAcoef[ie*nln2*4+iii] = aloc[a][i_c][b][j_c]*modified_detJ;
                        iii = iii + 1;
                    }
                }
            }
        }
    }
    
}

/*************************************************************************/
void SEMMTMaterial_jacobian(mxArray* plhs[], const mxArray* prhs[])
{
    
    double* dim_ptr = mxGetPr(prhs[0]);
    int dim     = (int)(dim_ptr[0]);
    int noe     = mxGetN(prhs[4]);
    double* nln_ptr = mxGetPr(prhs[5]);
    int nln     = (int)(nln_ptr[0]);
    int numRowsElements  = mxGetM(prhs[4]);
    int nln2    = nln*nln;
    
    plhs[0] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nln2*noe*dim*dim,1, mxREAL);
    
    double* myArows    = mxGetPr(plhs[0]);
    double* myAcols    = mxGetPr(plhs[1]);
    double* myAcoef    = mxGetPr(plhs[2]);
    
    int k,l;
    int q;
    int NumQuadPoints     = mxGetN(prhs[6]);
    int NumNodes          = (int)(mxGetM(prhs[3]) / dim);
    
    double* U_h   = mxGetPr(prhs[3]);
    double* w   = mxGetPr(prhs[6]);
    double* invjac = mxGetPr(prhs[7]);
    double* detjac = mxGetPr(prhs[8]);
    double* phi = mxGetPr(prhs[9]);
    double* gradrefphi = mxGetPr(prhs[10]);
    
    double gradphi[dim][nln][NumQuadPoints];
    double* elements  = mxGetPr(prhs[4]);
    
    double GradV[dim][dim];
    double GradU[dim][dim];
    double GradUh[dim][dim][NumQuadPoints];
    
    double Id[dim][dim];
    int d1,d2;
    for (d1 = 0; d1 < dim; d1 = d1 + 1 )
    {
        for (d2 = 0; d2 < dim; d2 = d2 + 1 )
        {
            Id[d1][d2] = 0;
            if (d1==d2)
            {
                Id[d1][d2] = 1;
            }
        }
    }
    
    double F[dim][dim];
    double EPS[dim][dim];
    double dP[dim][dim];
    double P_Uh[dim][dim];
    
    double* material_param = mxGetPr(prhs[2]);
    double Young = material_param[0];
    double Poisson = material_param[1];
    double Stiffening_power = material_param[2];
    double mu = Young / (2 + 2 * Poisson);
    double lambda =  Young * Poisson /( (1+Poisson) * (1-2*Poisson) );
    
    /* Assembly: loop over the elements */
    int ie;
    
#pragma omp parallel for shared(invjac,detjac,elements,myAcols,myArows,myAcoef,U_h) private(gradphi,F,EPS,dP,P_Uh,GradV,GradU,GradUh,ie,k,l,q,d1,d2) firstprivate(phi,gradrefphi,w,numRowsElements,nln2,nln,NumNodes,Id,mu,lambda)
    
    for (ie = 0; ie < noe; ie = ie + 1 )
    {
        for (k = 0; k < nln; k = k + 1 )
        {
            for (q = 0; q < NumQuadPoints; q = q + 1 )
            {
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    gradphi[d1][k][q] = 0;
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        gradphi[d1][k][q] = gradphi[d1][k][q] + INVJAC(ie,d1,d2)*GRADREFPHI(k,q,d2);
                    }
                }
            }
        }
        
        for (q = 0; q < NumQuadPoints; q = q + 1 )
        {
            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
            {
                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                {
                    GradUh[d1][d2][q] = 0;
                    for (k = 0; k < nln; k = k + 1 )
                    {
                        int e_k;
                        e_k = (int)(elements[ie*numRowsElements + k] + d1*NumNodes - 1);
                        GradUh[d1][d2][q] = GradUh[d1][d2][q] + U_h[e_k] * gradphi[d2][k][q];
                    }
                }
            }
        }
        
        int iii = 0;
        int ii = 0;
        int a, b, i_c, j_c;
        
        /* loop over test functions --> a */
        for (a = 0; a < nln; a = a + 1 )
        {
            /* loop over test components --> i_c */
            for (i_c = 0; i_c < dim; i_c = i_c + 1 )
            {
                /* set gradV to zero*/
                for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                {
                    for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                    {
                        GradV[d1][d2] = 0;
                    }
                }
                
                /* loop over trial functions --> b */
                for (b = 0; b < nln; b = b + 1 )
                {
                    /* loop over trial components --> j_c */
                    for (j_c = 0; j_c < dim; j_c = j_c + 1 )
                    {
                        /* set gradU to zero*/
                        for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                        {
                            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                            {
                                GradU[d1][d2] = 0;
                            }
                        }
                        
                        double aloc = 0;
                        for (q = 0; q < NumQuadPoints; q = q + 1 )
                        {
                            for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                            {
                                GradV[i_c][d2] = gradphi[d2][a][q];
                                GradU[j_c][d2] = gradphi[d2][b][q];
                            }
                            
                            
                            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                            {
                                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                                {
                                    F[d1][d2] = Id[d1][d2] + GradU[d1][d2];
                                }
                            }
                            
                            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                            {
                                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                                {
                                    EPS[d1][d2] = 0.5 * ( F[d1][d2] + F[d2][d1] ) - Id[d1][d2];
                                }
                            }
                            
                            
                            double trace = Trace(dim, EPS);
                            for (d1 = 0; d1 < dim; d1 = d1 + 1 )
                            {
                                for (d2 = 0; d2 < dim; d2 = d2 + 1 )
                                {
                                    dP[d1][d2] = 2 * mu * EPS[d1][d2] + lambda * trace * Id[d1][d2];
                                }
                            }
                            aloc  = aloc + Mdot( dim, GradV, dP) * w[q];
                        }
                        myArows[ie*nln2*dim*dim+iii] = elements[a+ie*numRowsElements] + i_c * NumNodes;
                        myAcols[ie*nln2*dim*dim+iii] = elements[b+ie*numRowsElements] + j_c * NumNodes;
                        myAcoef[ie*nln2*dim*dim+iii] = aloc * detjac[ie] * pow( detjac[0] / detjac[ie], Stiffening_power );
                        
                        iii = iii + 1;
                    }
                }
                
            }
        }
    }
        
}
/*************************************************************************/