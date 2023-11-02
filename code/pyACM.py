#-------------------------------------------
# ACM (MCA) returning object Dudi
# 2023-11-02 modified by pbady
#------------------------------------------
import pandas as pd
import numpy as np

def disjonctif(X):
  X_cat = X.astype("category")
  X_dis = pd.get_dummies(X_cat)
  X_dis = X_dis*1
  return X_dis;

def pyACM(X,lw=None,nf=2):
  Xdis = disjonctif(ours)
  m = Xdis.shape[1]
  n = Xdis.shape[0]
  v = ours.shape[1]
  if lw==None :
    lw = pd.DataFrame(np.repeat(1/n,n))[0]
  D = np.diag(lw)
  cw = np.dot(np.dot(Xdis.T,D), np.ones(n)) 
  Dm = np.diag(cw)
  X = np.dot(Xdis,np.diag(1/cw))-1
  cw = cw/v
  dudi = pyDudi(X,cw,lw,nf)    
  return dudi;
