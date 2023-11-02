#-------------------------------------------
# PCA returning object Dudi
# 2023-11-02 modified by pbady
#------------------------------------------
import pandas as pd
import numpy as np

def pyPCA(X,cw=None,lw=None,nf=2,center=True,scale=True):
  dim = X.shape
  n = dim[0]
  p = dim[1]  
  if center:
    X = X-X.mean()
  if scale:
    X = X/X.std(ddof=0)
  if lw==None :
    lw = pd.DataFrame(np.repeat(1/n,n))[0]
  if cw==None :
    cw = pd.DataFrame(np.repeat(1,p))[0]
  dudi = pyDudi(X,cw,lw,nf)  
  return dudi;

  
