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

#Calculate lift coefficient
cp_u = 0
cp_l = 0
#upper surface
for i in range(int((len(mid_panel)+1)/2),int(len(mid_panel)-1)):
    cp_u = cp_u + ((cp[i]+cp[i+1])*(mid_panel[i+1,0]-mid_panel[i,0]))/2
#lower surface
for i in range(int((len(mid_panel)+1)/2 - 1), 0, -1):
    cp_l = cp_l + ((cp[i]+cp[i-1])*(mid_panel[i-1,0]-mid_panel[i,0]))/2

Cl = (cp_l-cp_u)*np.cos(aoa)

print("Lift Coefficient = ",Cl)
zz = np.transpose(np.vstack((mid_panel[:,0],cp)))

plt.figure(1)
plt.subplot(211)
plt.plot(mid_panel[:,0],-cp,'bo-',label = "Pressure Coefficient")
plt.xlabel("Chord(x/X)");plt.ylabel("-Cp");plt.legend()
plt.subplot(212)
plt.plot(airfoil[:, 0], airfoil[:, 1],label = "Airfoil Geometry")
plt.xlabel("Chord(x/X)");plt.ylabel("Thickness");plt.legend()

plt.figure(2)
plt.quiver(grid[:,0],grid[:,1],vx,vy,V)
plt.plot(airfoil[:, 0], airfoil[:, 1])
plt.colorbar(label="Velocity (m/s)")
plt.show()

