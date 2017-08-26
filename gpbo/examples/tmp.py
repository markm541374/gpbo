import gpbo
from gpbo.core import GPdc as GPdc
import numpy as np

X = np.array([[0.],[1.]])
Y = np.array([[1.],[1.]])
S = np.array([[1e-6],[1e-6]])
D = [[np.NaN],[np.NaN]]
g = GPdc.GPcore(X,Y,S,D,GPdc.kernel(GPdc.SQUEXP,1,np.array([2.,0.5])))
print(GPdc.GP_LKonly(X,Y,S,D,GPdc.kernel(GPdc.SQUEXP,1,np.array([2.,0.5]))).llk())
print(g.llk())