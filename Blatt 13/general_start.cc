#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stddef.h>
#include <math.h>

using namespace std;

// Number of data points and parameters
const int ndata=10, ma=3;

/*********************************************************************/
/*                                                                   */
/* Evaluate base functions f_m(x) (associated with each paramter a_m)*/
/* for a single value of x                                           */
/* where the fit function f(x) = sum_m a_m * f_m(x)                  */
/* and m=0,....,ma-1                                                 */
/*                                                                   */
/* afuncs: vector containing values of f_m(x)                        */
/* ma: number of parameters                                          */
/*                                                                   */
/*********************************************************************/

void base_func(double x, double afuncs[], int ma){

  //?????????????????????
}

/*********************************************************************/
/*                                                                   */
/* Performs the general linear fit of function                       */
/* f(x) = sum_m a_m * f_m(x), where m=0,....,ma-1                    */
/* to the data points (x_i,y_i) contained in the vectors x[] and y[] */
/*                                                                   */
/* sig: contains the standard deviations of the y_i                  */
/* a[]: the ma fit parameters                                        */
/* covar[][]: the covariance matrix (diagonal elements give errors   */
/*                                                 on the parameters)*/
/* chisq: value of (minimal) chisq goodness of fit                   */
/*                                                                   */
/*                                                                   */
/*********************************************************************/

void lfit(double x[ndata], double y[ndata], double sig[ndata],  double a[ma], 
	  double covar[ma][ma], double &chisq)
{
  /* Construct the matrix alpha = A^T A 
     from the elements of the design matrix A */

  /* Construct the vector beta = A^T b */

  /* Use Gauss-Jordan elimination to invert alpha (to provide the
     covariance matrix covar) and provide the vector of parameters a=
     alpha^{-1} b */

  /* Compute the chisq */
}

int main(){

  double x[ndata], y[ndata], sig[ndata], a[ma];
  double covar[ma][ma], chi2;
  
  // Read in synthetic data with errors stored in vector sig
  ifstream in("data_gen.txt");
  if(in.is_open()){
    for(int i=0;i<ndata;i++){
      string line;
      getline(in,line);
      sscanf(line.c_str(),"%lf %lf %lf",&x[i],&y[i],&sig[i]);
      printf("%f %f %f\n",x[i],y[i],sig[i]);
    }
  }else{
    cerr<<"File data_gen.txt not opened"<<endl;
    exit(1);
  }

  //Perform the fit
  lfit(x,y,sig,a,covar,chi2);

  //print out the parameters, their errors and the chisq per degree of freedom
  
  return 0;
}
