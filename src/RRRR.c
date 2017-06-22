/*Kun Chen*/
/*Last modified: 08/30/2012*/
/*Coordinate descent algorithms*/
/*pcd,rssvd_orth,kronecker*/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>


int Sign(double beta);
double Max(double a1, double a2);
double Abs(double a1);
double pow( double base, double exp );
double vSum (double *X, int nX);
void matMply (double *X, int nrX , int ncX, double *Y, int nrY , int ncY, double *ans);
void vecSum (double *X, int nX, double *Y, int nY, double *ans);
void vecMinus (double *X, int nX, double *Y, int nY, double *ans);
void vecDiv (double *X, int nX, double *Y, int nY, double *ans);
void vecMply (double *X, int nX, double *Y, int nY, double *ans);
void vsMply (double *X, int nX, double scalar, double *ans);
void vPow (double *X, int nrX , double k, double *ans);

void vecThresh (double *X, int nX, double *Y, int nY, double *lambda, double *ans);
void vecStz (double *X, int nX, double *ans1 ,double *ans2);
void vecAbs (double *X, int nX, double *ans);



SEXP pcd(SEXP XX, SEXP Y, SEXP Lambda, SEXP Beta0, SEXP Conv, SEXP Maxiter)
{

  /*Read X info*/
  int *dimX;
  double *xxptr;  
  dimX = INTEGER( coerceVector(getAttrib(XX, R_DimSymbol), INTSXP) );
  PROTECT(XX = coerceVector(XX, REALSXP) ); 
  int nrX = dimX[0];
  int ncX = dimX[1];
  xxptr = REAL(XX);


  /*Read Y and other*/
  double *yptr;
  PROTECT(Y = coerceVector(Y, REALSXP) );
  int nY = length(Y);
  yptr = REAL(Y);

  double *lambdaptr, *convptr; 
  int *miterptr;
  lambdaptr = REAL(Lambda);
  convptr = REAL(Conv);
  miterptr = INTEGER(Maxiter);

  double *beta0ptr;  
  PROTECT(Beta0 = coerceVector(Beta0, REALSXP) );  
  beta0ptr = REAL(Beta0);  


  /*Check dimension*/
  if (nrX != nY) error("Dim not match");


  /*Standardize X*/
  SEXP X;
  double *xptr;
  PROTECT(X = allocMatrix(REALSXP, nrX, ncX) );
  xptr = REAL(X);    
  SEXP msqX;
  double *msqptr;
  PROTECT( msqX = allocVector(REALSXP, ncX) );
  msqptr = REAL(msqX);  
  int ii, jj;
  double sum;
  for ( jj =0; jj <ncX; jj ++){
    sum = 0;
    for ( ii = 0; ii <nrX; ii++){
      sum = sum + pow(xxptr[jj*nrX + ii],2);
    } 
    sum = pow(sum,0.5);
    msqptr[jj] = sum;
    for ( ii = 0; ii <nrX; ii++){
      xptr[jj*nrX + ii] = xxptr[jj*nrX + ii]/sum;
    }         
  }


  /* Create beta and its initial value*/
  SEXP beta;
  double *betaptr;
  PROTECT( beta = allocVector(REALSXP, ncX)) ;
  betaptr = REAL(beta); 
  vecMply(beta0ptr, ncX, msqptr, ncX, betaptr);  

  
  /* Create lambda vector*/
  SEXP lambda;
  double *lamptr;
  PROTECT( lambda = allocVector(REALSXP, ncX) ) ;
  lamptr = REAL(lambda);    
  for(ii=0; ii<ncX; ii++){
    lamptr[ii] = *lambdaptr/msqptr[ii];
  }     


  int counter=0;
  double diff=0;
  double *diffptr=&diff;
  *diffptr = *convptr*2;
   
  SEXP res;
  double *resptr;
  PROTECT(res = allocVector(REALSXP, nY));
  resptr = REAL(res);
  SEXP temp1;
  double *temp1ptr;
  PROTECT(temp1 = allocVector(REALSXP, nY));
  temp1ptr = REAL(temp1);
  matMply (xptr, nrX , ncX, betaptr, ncX, 1, temp1ptr);
  vecMinus (yptr, nY, temp1ptr, nrX, resptr);
  
  SEXP betap;
  double *betapptr;
  PROTECT( betap = allocVector(REALSXP, ncX) );
  betapptr = REAL(betap);
  SEXP temp2;
  double *temp2ptr;
  PROTECT( temp2 = allocVector(REALSXP, ncX) );
  temp2ptr = REAL(temp2);  

  double fit, new;
  double *fitptr = &fit;
  double *newptr = &new;  
  
  while ((counter < *miterptr) & (*diffptr > *convptr)){     
    for(ii=0; ii<ncX; ii++){   
        betapptr[ii] = betaptr[ii];    
        matMply(xptr+ii*nrX, 1, nrX, resptr, nrX, 1, fitptr);
        *fitptr += betaptr[ii];
        *newptr = Sign(*fitptr)*Max(0,Abs(*fitptr)-lamptr[ii]);        
        vsMply (xptr+ii*nrX, nrX, *newptr-betaptr[ii], temp1ptr);                
        vecMinus(resptr, nrX, temp1ptr, nrX, resptr);                        
        betaptr[ii] = *newptr;       
    }
    counter++;        
        
    vecMinus(betaptr,ncX,betapptr,ncX,temp2ptr);
    vPow(temp2ptr,ncX,2,temp2ptr);
    vPow(betapptr,ncX,2,betapptr);
    *diffptr = pow(vSum(temp2ptr,ncX)/vSum(betapptr,ncX),0.5);    
    
  }  
  
 
  UNPROTECT(11);
  
  vecDiv(betaptr,ncX,msqptr,ncX,betaptr);    
  return(beta);
}



SEXP rssvd_orth(SEXP YX, SEXP Wu, SEXP Wv, SEXP Wd, SEXP U0, SEXP Lambda, SEXP Conv, SEXP Maxiter)
{

  /*input Y%*%t(X), m by n matrix*/
  int *dimYX;
  double *yxptr;
  dimYX = INTEGER( coerceVector(getAttrib(YX, R_DimSymbol), INTSXP) );
  PROTECT(YX = coerceVector(YX, REALSXP) );
  int nrYX = dimYX[0];
  int ncYX = dimYX[1];
  yxptr = REAL(YX);


  double *wuptr;
  PROTECT(Wu = coerceVector(Wu, REALSXP) );
  wuptr = REAL(Wu);  

  double *wvptr;
  PROTECT(Wv = coerceVector(Wv, REALSXP) );
  wvptr = REAL(Wv);
  
  double *wdptr;
  PROTECT(Wd = coerceVector(Wd, REALSXP) );
  wdptr = REAL(Wd);

  double *u0ptr;
  PROTECT(U0 = coerceVector(U0, REALSXP) );
  u0ptr = REAL(U0);
  
  double *lambdaptr, *convptr;
  int *miterptr;
  lambdaptr = REAL(Lambda);
  convptr = REAL(Conv);
  miterptr = INTEGER(Maxiter);


  /* Create V and its initial value*/
  SEXP V;
  double *vptr;
  PROTECT( V = allocVector(REALSXP, ncYX)) ;
  vptr = REAL(V);


  /* Create U and its initial value*/
  SEXP U;
  double *uptr;
  PROTECT( U = allocVector(REALSXP, nrYX)) ;
  uptr = REAL(U);
  int ii;
  for(ii=0; ii<nrYX; ii++){
    uptr[ii] = u0ptr[ii];
  }



  /*Create D and its initial value*/
  double D=1;
  double *dptr=&D;
  
  
  
 
  int counter=0;
  double diff=0;
  double *diffptr=&diff;
  *diffptr = *convptr*2;

  SEXP Vtemp;
  double *vtptr;
  PROTECT(Vtemp = allocVector(REALSXP, ncYX));
  vtptr = REAL(Vtemp);
  
  SEXP Utemp;
  double *utptr;
  PROTECT(Utemp = allocVector(REALSXP, nrYX));
  utptr = REAL(Utemp); 

  SEXP Vp;
  double *vpptr;
  PROTECT(Vp = allocVector(REALSXP, ncYX));
  vpptr = REAL(Vp);
  
  SEXP Up;
  double *upptr;
  PROTECT(Up = allocVector(REALSXP, nrYX));
  upptr = REAL(Up); 

  double lamtemp=0;
  double *lamtempptr = &lamtemp;

  while ((counter < *miterptr) & (*diffptr > *convptr)){

    for(ii=0; ii<ncYX; ii++)
    {
      vpptr[ii] = vptr[ii];
    }
    for(ii=0; ii<nrYX; ii++)
    {
      upptr[ii] = uptr[ii];
    }  

    /*  V1  */
    matMply(uptr, 1, nrYX, yxptr, nrYX, ncYX, vtptr);
    /*  lambda_u  */
    vecAbs(uptr,nrYX,utptr);     
    matMply(wuptr,1,nrYX,utptr,nrYX,1,lamtempptr);
    lamtempptr[0] = lamtempptr[0]*lambdaptr[0]*wdptr[0];
    /*  D*V  */
    vecThresh(vtptr, ncYX, wvptr, ncYX, lamtempptr,vptr);
    /*  standardization  */
    vecStz (vptr, ncYX, vptr, dptr);  
  
    /*  U1  */
    matMply(yxptr, nrYX, ncYX, vptr, ncYX, 1, utptr);
    /*  lambda_v  */
    vecAbs(vptr,ncYX,vtptr);                
    matMply(wvptr,1,ncYX,vtptr,ncYX,1,lamtempptr);
    lamtempptr[0] = lamtempptr[0]*lambdaptr[0]*wdptr[0];
    /*  D*U  */
    vecThresh(utptr, nrYX, wuptr, nrYX, lamtempptr,uptr);
    /*  standardization  */
    vecStz (uptr, nrYX, uptr, dptr);
    
    counter++;

    vecMinus(vptr,ncYX,vpptr,ncYX,vtptr);
    vPow(vtptr,ncYX,2,vtptr);
    vPow(vpptr,ncYX,2,vpptr);
    vecMinus(uptr,nrYX,upptr,nrYX,utptr);
    vPow(utptr,nrYX,2,utptr);
    vPow(upptr,nrYX,2,upptr);    
    *diffptr = pow(vSum(vtptr,ncYX)/vSum(vpptr,ncYX),0.5)+pow(vSum(utptr,nrYX)/vSum(upptr,nrYX),0.5);

    /*printf("%f\n",diff);*/
  }

  SEXP Result;
  double *resptr;
  PROTECT(Result = allocVector(REALSXP, nrYX+ncYX+1));
  resptr = REAL(Result);
  
  for(ii=0; ii < nrYX; ii++)
  {
    resptr[ii] = uptr[ii];
  }
  for(ii=0; ii < ncYX; ii++)
  {
    resptr[ii+nrYX] = vptr[ii];
  }    
  resptr[nrYX+ncYX] = *dptr;  

  UNPROTECT(12);

  return(Result);
  
}







int Sign(double beta){
  int sign=0;
  if(beta >0){
    sign = 1;
  }
  if(beta<0){
    sign = -1;
  }
  if(beta==0){
    sign = 0;
  }
  return(sign);
}
  
    
    
    
double Max(double a1, double a2){
  double a=0; 
  if(a1 >= a2){
    a=a1;
  }
  if(a1 < a2){
    a=a2;
  }
  return(a);
}




double Abs(double a1)
{
  double a=a1; 
  if(a1 < 0){
    a=-a1;
  }
  return(a);
}




void matMply (double *X, int nrX , int ncX, double *Y, int nrY , int ncY, double *ans)
{
  double sum;
  int ii, jj, kk;
  for ( ii =0; ii <nrX; ii ++){
    for ( jj =0; jj <ncY; jj++){
      sum = 0 ;
      for (kk=0; kk<ncX; kk++){
        sum = sum + X[ii+nrX*kk]*Y[kk+nrY*jj];
      }  
      ans[ii+nrX*jj] = sum;
    }
  }
}




void vPow (double *X, int nrX , double k, double *ans)
{
  int ii;
  for ( ii=0; ii<nrX; ii++){
    ans[ii] = pow(X[ii],k);
  }
}




double vSum (double *X, int nX)
{
  int ii;
  double sum=0;
  for ( ii=0; ii<nX; ii++){
    sum = sum + X[ii];
  }
  return(sum);
}




void vecSum (double *X, int nX, double *Y, int nY, double *ans)
{
  if (nX != nY) error("Dim not match");
  int ii;
  for (ii=0; ii<nX; ii++){
    ans[ii] = X[ii] + Y[ii];
  }
}




void vecMinus (double *X, int nX, double *Y, int nY, double *ans)
{
  if (nX != nY) error("Dim not match");
  int ii;
  for (ii=0; ii<nX; ii++){
    ans[ii] = X[ii] - Y[ii];
  }
}




void vecDiv (double *X, int nX, double *Y, int nY, double *ans)
{
  if (nX != nY) error("Dim not match");
  int ii;
  for (ii=0; ii<nX; ii++){    
    ans[ii] = X[ii] / Y[ii];
  }
}




void vecMply (double *X, int nX, double *Y, int nY, double *ans)
{
  if (nX != nY) error("Dim not match");
  int ii;
  for (ii=0; ii<nX; ii++){    
    ans[ii] = X[ii] * Y[ii];
  }
}




void vsMply (double *X, int nX, double scalar, double *ans)
{
  int ii;
  for (ii=0; ii<nX; ii++){    
    ans[ii] = X[ii] * scalar;
  }
}







void vecAbs (double *X, int nX, double *ans)
{
  int ii; 
  for (ii=0; ii<nX; ii++){
     ans[ii] = fabs(X[ii]);
  }
}



void vecThresh (double *X, int nX, double *Y, int nY, double *lambda, double *ans)
{
  if (nX != nY) error("Dim not match");
  int ii; 
  double diff=0;
  for (ii=0; ii<nX; ii++){
    diff = fabs(X[ii]) - (*lambda)*Y[ii];
    if(diff < 0)
    {
      ans[ii] = 0;
    }
    else
    {
      ans[ii] = copysign(diff,X[ii]);
    }
  }
}



void vecStz (double *X, int nX, double *ans1 ,double *ans2)
{
  int ii; 
  double norm=0;
  for (ii=0; ii<nX; ii++){
     norm = norm + pow(X[ii],2);
  }
  *ans2 = sqrt(norm);
  for (ii=0; ii<nX; ii++){
     ans1[ii] = X[ii]/(*ans2);
  }    
}





SEXP kronecker(SEXP A, SEXP B)
{

  /*Read A info*/
  int *dimA;
  double *aptr;  
  dimA = INTEGER( coerceVector(getAttrib(A, R_DimSymbol), INTSXP) );
  PROTECT(A = coerceVector(A, REALSXP) ); 
  int nrA = dimA[0];
  int ncA = dimA[1];
  aptr = REAL(A);

  /*Read B info*/
  int *dimB;
  double *bptr;  
  dimB = INTEGER( coerceVector(getAttrib(B, R_DimSymbol), INTSXP) );
  PROTECT(B = coerceVector(B, REALSXP) ); 
  int nrB = dimB[0];
  int ncB = dimB[1];
  bptr = REAL(B);

  SEXP C;
  double *cptr;
  PROTECT(C = allocMatrix(REALSXP, nrA*nrB, ncA*ncB) );
  cptr = REAL(C);

  
    int i,j,k,l;

    for(i=0; i<nrA; i++)
    {
        for(k=0; k<nrB; k++)
            {
                for(j=0; j<ncA; j++)
                    {
                        for(l=0; l<ncB; l++)
                            {
                                 
                                cptr[i*nrB+k+nrA*nrB*(j*ncB+l)] = aptr[i+nrA*j]*bptr[k+nrB*l];
                                
                            }
                    }
            } 
    }
     
    UNPROTECT(3); 
    return(C);
} 




/*multiply each column of a matrix*/
SEXP matcolMply(SEXP A, SEXP b)
{

  /*Read A info*/
  int *dimA;
  double *aptr;  
  dimA = INTEGER( coerceVector(getAttrib(A, R_DimSymbol), INTSXP) );
  PROTECT(A = coerceVector(A, REALSXP) ); 
  int nrA = dimA[0];
  int ncA = dimA[1];
  aptr = REAL(A);

  /*Read b */
  double *bptr;
  PROTECT(b = coerceVector(b, REALSXP) );
  int nb = length(b);
  bptr = REAL(b);

  /*Check dimension*/
  if (ncA != nb) error("Dim not match");


  SEXP C;
  double *cptr;
  PROTECT(C = allocMatrix(REALSXP, nrA, ncA) );
  cptr = REAL(C);


  int i,j;

  for(i=0; i<nrA; i++)
  {
    for(j=0; j<ncA; j++)
    {
        cptr[i+nrA*j] = aptr[i+nrA*j]*bptr[j];
    }
  }  
    
  UNPROTECT(3); 
  return(C);
  
}

