#-------------------------------------------
# Class Dudi and function to compute object Dudi
# 2023-11-02 modified by pbady
#------------------------------------------
import pandas as pd
import numpy as np

class Dudi:
  def __init__(self,tab,eig,rank,nf,c1,co,l1,li):
    self.tab = tab
    self.eig = eig
    self.rank = rank
    self.nf = nf
    self.c1 = c1
    self.co = co
    self.l1 = l1
    self.li = li

def pyDudi(X,cw,lw,nf):
  dim = X.shape
  n = dim[0]
  p = dim[1]
  nf0 = nf-1
  # n=len(X)
  # p=len(X.columns)
  D = np.diag(np.sqrt(lw))
  Q = np.diag(np.sqrt(cw))
  
  # XtDXQ => problem  with Q !!!
  XD = np.dot(X.T,D).T
  XD = np.dot(XD,Q)
  XtX = np.dot(XD.T,XD)
  
  # decomposition
  eigenvalues, eigenvectors = np.linalg.eig(XtX)
  index = np.argsort(eigenvalues)[::-1]
  
  # np.nonzero(eigenvalues)[0]
  eigenvalues = eigenvalues.real
  eigenvectors = eigenvectors.real
  eigenvalues=eigenvalues[index]
  eigenvectors=eigenvectors[:,index]
  # results
  rank = len(np.nonzero(eigenvalues)[0])
  C1 = np.dot(np.diag(1/np.sqrt(cw)),eigenvectors[:,0:nf])
  #C1 = eigenvectors[:,0:nf]
  XQ = np.dot(X,np.diag(cw))
  Li = np.dot(XQ, C1)
  
  # need to adjust the weighting (problem with sqrt)
  L1 = np.dot(Li,np.diag(1/np.sqrt(eigenvalues[0:nf])))
  Co = np.dot(C1,np.diag(np.sqrt(eigenvalues[0:nf])))
  return Dudi(X,eigenvalues,rank,nf,C1,Co,L1,Li);
