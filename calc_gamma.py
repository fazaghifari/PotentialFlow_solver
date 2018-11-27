import numpy as np
from numpy.linalg import solve
import matplotlib.pyplot as plt

def calgamma(airfoil,mid_panel,length_panel,theta,sine, cosine, rhs, aoa ):
    m = len(airfoil)-1
    mplus1 = m+1
    midx = mid_panel[:,0]; midy = mid_panel[:,1]
    foilx = airfoil[:,0]; foily = airfoil[:,1]
    cn1 = np.zeros(shape=[m, m])
    cn2 = np.zeros(shape=[m, m])
    ct1 = np.zeros(shape=[m, m])
    ct2 = np.zeros(shape=[m, m])
    v = np.zeros(shape=[m])
    cp = np.zeros(shape=[m])

    for i in range (0,m):
        for j in range (0,m):
            if i == j:
                #kutta condition
                cn1[i, j] = -1
                cn2[i, j] = 1
                ct1[i, j] = np.pi / 2
                ct2[i, j] = np.pi / 2
            else:
                a = -(midx[i]-foilx[j])*cosine[j] - (midy[i]-foily[j])*sine[j]
                b = (midx[i]-foilx[j])**2 + (midy[i]-foily[j])**2
                c = np.sin(theta[i]-theta[j])
                d = np.cos(theta[i]-theta[j])
                e = (midx[i]-foilx[j])*sine[j] - (midy[i]-foily[j])*cosine[j]
                f = np.log(1+length_panel[j]*(length_panel[j]+2*a)/b)
                g = np.arctan2(e*length_panel[j],b+a*length_panel[j])
                p = (midx[i]-foilx[j])*np.sin(theta[i]-2*theta[j]) + (midy[i]-foily[j])*np.cos(theta[i]-2*theta[j])
                q = (midx[i]-foilx[j])*np.cos(theta[i]-2*theta[j]) - (midy[i]-foily[j])*np.sin(theta[i]-2*theta[j])
                cn2[i,j] = d + 0.5*q*f/length_panel[j] - (a*c+d*e)*g/length_panel[j]
                cn1[i,j] = 0.5*d*f + c*g - cn2[i,j]
                ct2[i,j] = c + 0.5*p*f/length_panel[j] + (a*d-c*e)*g/length_panel[j]
                ct1[i,j] = 0.5*c*f - d*g - ct2[i,j]

    # Compute coefficients influence eq(5.47) and (5.49)
    an = np.zeros(shape=[mplus1,mplus1]); at = np.zeros(shape=[m,mplus1])
    for i in range(0,m):
        an[i,0] = cn1[i,0]
        an[i,m] = cn2[i,m-1]
        at[i,0] = ct1[i,0]
        at[i,m] = ct2[i,m-1]
        for j in range(1,m):
            an[i,j] = cn1[i,j] + cn2[i,j-1]
            at[i,j] = ct1[i,j] + ct2[i,j-1]
    an[m,0] = 1
    an[m,m] = 1
    for j in range(1,m):
        an[m,j] = 0
    rhs[m]=0

    #solve eq(5.47)
    #compute dimensionless velo and press coeff at control points
    gamma = solve(an,np.transpose(rhs))
    for i in range(0,m):
        v[i] = np.cos(theta[i]-aoa)
        for j in range(0,mplus1):
            v[i] = v[i] + at[i,j]*gamma[j]
            cp[i] = -(1-v[i]**2)  #minus sign is for flipping the plot

    plt.plot(midx,cp)
    plt.show()
    return (gamma,cp)
