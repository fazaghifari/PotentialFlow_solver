import numpy as np

def params (airfoilcoord,aoa):
    m = len(airfoilcoord)-1
    theta = np.zeros(shape=[m]); sine = np.zeros(shape=[m]); cosine = np.zeros(shape=[m]); rhs = np.zeros(shape=[m+1])
    for i in range(0,m):
        inext = i+1 #next panel index
        theta[i] = np.arctan2((airfoilcoord[inext,1]- airfoilcoord[i,1]),(airfoilcoord[inext,0]- airfoilcoord[i,0]))
        sine [i] = np.sin(theta[i])
        cosine [i] = np.cos(theta[i])
        rhs [i] = np.sin(theta[i]-aoa)
    return (theta,sine,cosine,rhs)
