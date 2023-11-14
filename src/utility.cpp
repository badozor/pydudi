/* ==============================================================================
# Title: utility.cpp
# Description:
# Author:  Pierre Bady <pierre.bady@unil.ch>
# Date : Jun 5, 2015 10:24:19 AM
# Version: 0.1
# Revision:  2023-11-10 modified by pbady
# Comments: RAS
# License: GPL version 2 or newer
# Copyright (C) 2015-2023  Pierre Bady
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# ============================================================================== */
 
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <cmath>
#include <assert.h>
#define TOL 1e-7
 
using namespace Rcpp;
using namespace arma;
 
 // additional test functions
//[[Rcpp::export]]
arma::vec randvec (arma::vec y){
  return shuffle(y);
}
//[[Rcpp::export]]
arma::mat randmat(arma::mat y,int dim){
  return shuffle(y,dim);
}
//[[Rcpp::export]]
arma::colvec Cppfunc(arma::vec x,arma::vec y, Function f){
 List fit = f(x,y);
 arma::colvec  yfit = as<arma::colvec>(fit[1]);
  return yfit;
}
//[[Rcpp::export]]
double AUC(arma::vec x, arma::vec y){
  int len,i;
  double a;
  len= x.size();
  a = 0.0;
  for(i=0;i< (len-1);i++){
    a+=fabs(x(i)-x(i+1))*(fabs(y(i))+fabs(y(i+1)))/2;
  }
  return a;
}
// test coefficient de correlation
// in arma => cor(X,Y, 1); */
//https://arma.sourceforge.net/docs.html#cor
//https://arma.sourceforge.net/docs.html
 
//[[Rcpp::export]]
double CORR(arma::vec x,arma::vec y){
  int i,n;
  n= x.size();
  double sumx,sumy,sumxx,sumyx,sumyy,r;
  sumx = 0.0;
  sumy = 0.0;
  sumxx = 0.0;
  sumyy = 0.0;
  sumyx = 0.0;
  r=0.0;
  for(i=0;i < n;i++){
    sumx += x(i);
    sumy += y(i);
    sumxx += x(i)*x(i);
    sumyx += y(i)*x(i);
    sumyy += y(i)*y(i);
    }
  r = (n*sumyx-sumx*sumy)/(sqrt(n*sumxx-sumx*sumx)*sqrt(n*sumyy - sumy*sumy));
  return r;
}
//[[Rcpp::export]]
List cdudi(arma::mat X,arma::vec cw,arma::vec lw,int nf){
  int i,j,n,p,nf0;
  n = X.n_rows;
  p = X.n_cols;
  nf0=nf-1; 
  //arma::uvec colindex=linspace(0,nf0,1);
  //weighting
  arma::mat D=diagmat(sqrt(lw));
  arma::mat Q=diagmat(sqrt(cw));
  arma::vec eigval;
  arma::mat eigvec;
  //XtDXQ
  arma::mat XD = X.t()*D;
  arma::mat XDQ = XD.t()*Q;
  arma::mat XTX = XDQ.t()*XDQ;
  
  //decompoistion
  arma::eig_sym(eigval,eigvec,XTX,"std");
  // sort and reduce
  arma::uvec index = sort_index(eigval,"descend");
  eigval = eigval.elem(index);
  eigvec = eigvec.cols(index);
  //coordinates
  arma::mat sqrtQ = diagmat(1/sqrt(cw));
  arma::mat c1 = sqrtQ*eigvec.cols(0,nf0);
  arma::mat li = (X*Q)*c1;
  arma::mat l1 = li*diagmat(1/sqrt(eigval.subvec(0,nf0)));
  arma::mat co = c1*diagmat(sqrt(eigval.subvec(0,nf0)));
  return List::create(_["X"]=X,_["cw"]=cw,_["lw"]=lw,_["eig"]=eigval,_["nf"]=nf,_["c1"]=c1,_["co"]=co,_["l1"]=l1,_["li"]=li);
}
//[[Rcpp::export]]
List cpca(arma::mat X,arma::vec cw,arma::vec lw,int nf,int center,int scale){
  int i,j,n,p,nf0;
  n = X.n_rows;
  p = X.n_cols;
  nf0=nf-1; 
  arma::vec meanX(p);
  arma::vec sdX(p);
  if(center==1){
    for(i=0;i<p;i++){
      meanX(i) = arma::sum(lw % X.col(i))/arma::sum(lw);
      X.col(i) = X.col(i)-meanX(i);
     }
  }
  if(scale==1){
    for(i=0;i<p;i++){
      sdX(i) = sqrt(arma::sum(lw % X.col(i) % X.col(i))/arma::sum(lw));
      X.col(i) = X.col(i)/sdX(i);
    }
  }
  return cdudi(X,cw,lw,nf);
}


