import numpy as np

def calcv (grid,mid_panel,airfoil,theta,thetat,sine,cosine,l_pan,gamma,aoa,vfree):
    npoints =  len(grid)-1
    m = len(airfoil)-1
    midx = mid_panel[:, 0];midy = mid_panel[:, 1]
    foilx = airfoil[:, 0];foily = airfoil[:, 1]
    gridx = grid[:,0]; gridy = grid[:,1]
    Vn = np.zeros(shape=[len(grid)])
    Vn_gamma = np.zeros(shape=[len(grid)])
    Vn_inf = np.zeros(shape=[len(grid)])
    Vt = np.zeros(shape=[len(grid)])
    Vt_gamma = np.zeros(shape=[len(grid)])
    Vt_inf = np.zeros(shape=[len(grid)])

    for i in range(0,npoints):
        # calculate potential
        tempn = 0
        tempt = 0
        for j in range(0,m):
            tempn1 = 0
            tempt1 = 0
            tempn2 = 0
            tempt2 = 0
            a = -(gridx[i] - foilx[j]) * cosine[j] - (gridy[i] - foily[j]) * sine[j]
            b = (gridx[i] - foilx[j]) ** 2 + (gridy[i] - foily[j]) ** 2
            c = np.sin(thetat[i] - theta[j])
            d = np.cos(thetat[i] - theta[j])
            e = (gridx[i] - foilx[j]) * sine[j] - (gridy[i] - foily[j]) * cosine[j]
            f = np.log(1 + l_pan[j] * (l_pan[j] + 2 * a) / b)
            g = np.arctan2(e * l_pan[j], b + a * l_pan[j])
            p = (gridx[i]-foilx[j])*np.sin(thetat[i]-2*theta[j])+(gridy[i]-foily[j])*np.cos(thetat[i]-2*theta[j])
            q = (gridx[i]-foilx[j])*np.cos(thetat[i]-2*theta[j])-(gridy[i]-foily[j])*np.sin(thetat[i]-2*theta[j])
            tempn2 = tempn2 + d + 0.5*q*f/l_pan[j] - (a*c+d*e)*g/l_pan[j]
            tempn1 = tempn1 + 0.5*d*f + c*g - tempn2
            tempt2 = tempt2 + c + 0.5*p*f/l_pan[j] + (a*d-c*e)*g/l_pan[j]
            tempt1 = tempt1 + 0.5*c*f - d*g - tempt2
            tempt = tempt + tempt1*gamma[j] + tempt2*gamma[j+1]
            tempn = tempn + tempn1*gamma[j] + tempn2*gamma[j+1]

        # Calculate Velocity
        Vt_gamma[i] = tempt
        Vn_gamma[i] = tempn
        Vt_inf[i] = np.cos(thetat[i]-aoa)
        Vn_inf[i] = -np.sin(thetat[i]-aoa)
        Vn[i] = vfree*(Vn_inf[i]+Vn_gamma[i])
        Vt[i] = vfree*(Vt_inf[i]+Vt_gamma[i])
        V = np.sqrt(Vn**2 + Vt**2)

    return (Vn,Vt,V)