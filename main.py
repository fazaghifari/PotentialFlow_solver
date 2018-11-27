import numpy as np
from copy import deepcopy
import dataimport as di
import geometrical_param as geompar

## Input variables ##
print("Input Variables")

#Freestream input
vfree = int(input("Freestream velocity (1-75 m/s) : "))
while vfree < 1 or vfree > 75:
    print("Invalid input !")
    vfree = int(input("Freestream velocity (1-75 m/s) : "))

#Angle of Attack Input
aoa = int(input("Input angle of attack (-3 to 10) : "))
while aoa < -3 or aoa > 10:
    print("Invalid input !")
    aoa = int(input("Input angle of attack (-3 to 10) : "))

#Input airfoil file
filename = str(input("Input airfoil coordinates eg.(2410.txt) : "))

## Main Program ##
mid_panel, normal_panel, length_panel, airfoil = di.importer(filename)
print("hi!")
