import numpy as np

def params (airfoilcoord,aoa):
    m = len(airfoilcoord)-1
    theta = np.zeros(shape=[m])
    for i in range(0,m):
        inext = i+1 #next panel index
        theta[i] = np.arctan2((airfoilcoord[inext,1]- airfoilcoord[i,1]),(airfoilcoord[inext,0]- airfoilcoord[i,0]))
        
