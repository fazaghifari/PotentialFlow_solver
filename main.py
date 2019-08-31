from analyzer import Airfoil

vfree = 1
aoa = 3
inputfile = "2410.txt"

airfoil = Airfoil(vfree, aoa, inputfile)
cl = airfoil.analyze()
airfoil.visualize()
