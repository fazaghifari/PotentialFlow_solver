import numpy as np
from copy import deepcopy
import dataimport as di
import geometrical_param as geompar
import  calc_gamma as cgm
import gridgen_fcn as gg

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
grid = gg.pointsgen(filename, domsize, npoints)
print("cp")
