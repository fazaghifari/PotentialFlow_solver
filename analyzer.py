import numpy as np
import matplotlib.pyplot as plt
import dataimport as di
import geometrical_param as geompar
import calc_gamma as cgm
import gridgen_fcn as gg
import dom_velo as dvel


class Airfoil:
    """
    Create airfoil model analysis

    Args:
        vfree (float/int): Freestream velocity in m/s
        aoadeg (float/int): Angle of attack in degree.
        filename (str): Path to airfoil coordinates.

    """

    def __init__(self, vfree, aoadeg, filename):
        """
        Initialize model

        Args:
            vfree (float/int): Freestream velocity in m/s
            aoadeg (float/int): Angle of attack in degree.
            filename (str): Path to airfoil coordinates.
        """
        self.vfree = vfree
        self.aoa = aoadeg*np.pi/180
        self.filename = filename
        self.domsize = 1
        self.npoints = 20
        self.mid_panel, self.normal_panel, self.length_panel, self.airfoil = di.importer(self.filename)
        self.cp = None
        self.grid = None
        self.vy = None
        self.vx = None
        self.V = None

    def analyze(self):
        """
        Performs main calculation

        Returns:
             Cl (float): Lift coefficient
        """
        theta, sine, cosine, rhs = geompar.params(self.airfoil, self.aoa)
        gamma, self.cp = cgm.calgamma(self.airfoil, self.mid_panel, self.length_panel,
                                      theta, sine, cosine, rhs, self.aoa)
        self.grid, thetat = gg.pointsgen(self.filename, self.domsize, self.npoints)
        self.vy, self.vx, self.V = dvel.calcv(self.grid, self.mid_panel, self.airfoil, theta, thetat, sine, cosine,
                                              self.length_panel, gamma, self.aoa, self.vfree)

        # Calculate lift coefficient
        cp_u = 0
        cp_l = 0
        # upper surface
        for i in range(int((len(self.mid_panel) + 1) / 2), int(len(self.mid_panel) - 1)):
            cp_u = cp_u + ((self.cp[i] + self.cp[i + 1]) * (self.mid_panel[i + 1, 0] - self.mid_panel[i, 0])) / 2
        # lower surface
        for i in range(int((len(self.mid_panel) + 1) / 2 - 1), 0, -1):
            cp_l = cp_l + ((self.cp[i] + self.cp[i - 1]) * (self.mid_panel[i - 1, 0] - self.mid_panel[i, 0])) / 2

        cl = (cp_l - cp_u) * np.cos(self.aoa)

        return cl

    def visualize(self):
        """
        Visualize the Flow field and the pressure distribution

        Return:
             None
        """

        plt.figure(1)
        plt.subplot(211)
        plt.plot(self.mid_panel[:, 0], -self.cp, 'bo-', label="Pressure Coefficient")
        plt.xlabel("Chord(x/X)")
        plt.ylabel("-Cp")
        plt.legend()
        plt.subplot(212)
        plt.plot(self.airfoil[:, 0], self.airfoil[:, 1], label="Airfoil Geometry")
        plt.xlabel("Chord(x/X)")
        plt.ylabel("Thickness")
        plt.legend()

        plt.figure(2)
        plt.quiver(self.grid[:, 0], self.grid[:, 1], self.vx, self.vy, self.V)
        plt.plot(self.airfoil[:, 0], self.airfoil[:, 1])
        plt.colorbar(label="Velocity (m/s)")
        plt.show()
