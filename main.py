import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import dataimport as di
import geometrical_param as geompar
import  calc_gamma as cgm
import gridgen_fcn as gg
import dom_velo as dvel

## Input variables ##
print("Input Variables")
# input grid
domsize = 1
npoints = 20

#Freestream input
vfree = int(input("Freestream velocity (1-75 m/s) : "))
while vfree < 1 or vfree > 75:
    print("Invalid input !")
    vfree = int(input("Freestream velocity (1-75 m/s) : "))

#Angle of Attack Input
aoadeg = int(input("Input angle of attack (-3 to 10) : "))
while aoadeg < -3 or aoadeg > 10:
    print("Invalid input !")
    aoadeg = int(input("Input angle of attack (-3 to 10) : "))
aoa = aoadeg*np.pi/180

#Input airfoil file
filename = str(input("Input airfoil coordinates eg.(2410.txt) : "))

## Main Program ##
mid_panel, normal_panel, length_panel, airfoil = di.importer(filename)
theta, sine, cosine, rhs = geompar.params(airfoil,aoa)
gamma, cp = cgm.calgamma(airfoil,mid_panel,length_panel,theta,sine,cosine,rhs,aoa)
grid,thetat = gg.pointsgen(filename, domsize, npoints)
vy,vx,V = dvel.calcv(grid,mid_panel,airfoil,theta,thetat,sine,cosine,length_panel,gamma,aoa,vfree)

plt.figure(1)
plt.scatter(grid[:,0],grid[:,1],s=1)

plt.figure(2)
plt.plot(mid_panel[:,0],cp,'bo-')

plt.figure(3)
plt.quiver(grid[:,0],grid[:,1],vx,vy,V)
plt.plot(airfoil[:, 0], airfoil[:, 1])
plt.show()

