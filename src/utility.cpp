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
 List scalewt(){
   
   
   
   
   
 }