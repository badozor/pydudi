#-------------------------------------------
# COA returning object Dudi
# 2023-11-02 modified by pbady
#------------------------------------------
import pandas as pd
import numpy as np

def pyCOA(X,nf=2):
  sumX = X.sum().sum()
  sumCol =  X.sum(axis=0)
  sumRow = X.sum(axis=1)
  pij = X/sumX
  pi = sumRow/sumX
  pj = sumCol/sumX
  Dj = np.diag(1/pj)
  Di = np.diag(1/pi)
  Z = np.dot(Di,pij)
  Z = np.dot(Z,Dj)
  Z = Z - 1
  Z = np.nan_to_num(Z)
  dudi = pyDudi(Z,pj,pi,nf)  
  return dudi;
