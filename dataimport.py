import numpy as np
import matplotlib.pyplot as plt

# function for airfoil plotting
def plotter(file):
    airfoil_coord = np.loadtxt(file, dtype=float)
    plt.plot(airfoil_coord[:, 0], airfoil_coord[:, 1])
    plt.show()
    print(airfoil_coord)


# function for inputting airfoil data to main program
def importer(file):
    airfoil_coord = np.loadtxt(file, dtype=float)
    foilx_coord = airfoil_coord[:, 0]
    foily_coord = airfoil_coord[:, 1]
    airfoil_coord = airfoil_coord[::-1]

    #variable initialization
    mid_panel = np.zeros(shape=[np.size(airfoil_coord, axis=0)-1, 2]) #location of panel's midpoint
    normal_panel = np.zeros(shape=np.shape(mid_panel))  #normal vector
    length_panel = np.zeros(shape=[np.size(mid_panel, 0),]) #length of each panel

    for i in range(0, len(mid_panel)):
        mid_panel[i, :] = (airfoil_coord[i, :] + airfoil_coord[i+1, :]) / 2

        # Calculate normal vector for each panel
        vec_panel = airfoil_coord[i+1, :] - airfoil_coord[i, :]
        res = np.cross(np.append(vec_panel, 0), np.array([0, 0, 1]))
        normal_panel[i, :] = res[0:2]
        length_panel[i] = np.linalg.norm(vec_panel)

    # normalization of normal vector
    norm_panel = np.linalg.norm(normal_panel, axis=1)/(0.01*length_panel)
    temp_norm = np.column_stack([norm_panel, norm_panel])
    normal_panel[:, :] = np.divide(normal_panel, temp_norm)


    del norm_panel, temp_norm
    return (mid_panel, normal_panel, length_panel, airfoil_coord)