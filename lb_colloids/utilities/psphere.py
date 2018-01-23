import numpy as np
import random
import math
import PIL
import matplotlib.pyplot as plt


class PSphere(object):
    """

    """
    def __init__(self, radius=20, porosity=0.5, dimension=256, sensitivity=0.08):
        self.radius = radius
        self.porosity = porosity
        self.sensitivity = sensitivity
        self.dimension = dimension
        self.matrix = np.ones((dimension, dimension), dtype=bool)
        self.matrix_porosity = 0.
        self.particle_space = False
        self.pore_space = True
        self.generate_plane()
        self.check_porosity()
        self.check_percolation()


    def generate_plane(self):
        """

        :param porosity:
        :param radius:
        :param dimension:
        :return:
        """
        porosity = self.porosity
        slice_location = self.dimension / 2
        low_bound = slice_location - int(self.radius)
        up_bound = slice_location + int(self.radius)

        if low_bound <= 0 or up_bound > self.dimension:
            raise AssertionError("Radius to large or slice location incorrect")

        relative_radius = self.radius / float(self.dimension)
        number_of_spheres = self.iround(-3.0 * np.log(self.porosity) /
                                        (4 * np.pi * relative_radius ** 3))

        for i in range(number_of_spheres):
            z = 1 + random.uniform(0, self.dimension)
            if up_bound > z > low_bound:
                x = 1 + int(random.uniform(0, self.dimension))
                y = 1 + int(random.uniform(0, self.dimension))
                slice_distance = abs(z - slice_location)
                slice_radius = np.sqrt(self.radius ** 2 - slice_distance ** 2) - 0.5
                if slice_radius < 0 or np.isnan(slice_radius):
                    slice_radius = 0
                # todo: continue with patch function()
                self.patch(x, y, slice_radius)

    def patch(self, x, y, radius):
        """

        :param x:
        :param y:
        :param radius:
        :return:
        """
        circumference = 2 * np.pi * radius + 1
        nsteps = int(math.ceil(circumference))/4

        if nsteps == 1:
            nsteps = 2

        for i in range(nsteps):
            angle = 0.5 * np.pi * float(i)/float(nsteps)
            ixpos = self.iround(radius * math.cos(angle))
            iypos = self.iround(radius * math.sin(angle))

            for ix in range(-ixpos, ixpos + 1):
                newxpos = x + ix
                newypos = y + iypos

                if newxpos >= self.matrix.shape[1]:
                    newxpos -= self.matrix.shape[1]

                elif newxpos < 0:
                    newxpos += self.matrix.shape[1]

                if newypos >= self.matrix.shape[0]:
                    newypos -= self.matrix.shape[0]

                elif newypos < 0:
                    newypos += self.matrix.shape[0]

                self.matrix[newypos][newxpos] = self.particle_space

                newypos = y - iypos

                if newypos >= self.matrix.shape[0]:
                    newypos -= self.matrix.shape[0]

                elif newypos < 0:
                    newypos += self.matrix.shape[0]

                self.matrix[newypos][newxpos] = self.particle_space



    def burn_2d(self):
        """

        :return:
        """
        return

    def iround(self, val):
        """

        :param val:
        :return:
        """
        if val - math.floor(val) >= 0.5:
            return int(math.ceil(val))
        else:
            return int(math.floor(val))

    def check_percolation(self):
        sweep = np.zeros((self.matrix.shape), dtype=bool)
        sweep[0] = self.matrix[0]

        for i, line in enumerate(self.matrix[1:]):
            temp = []
            for j, val in enumerate(line):
                if not val:
                    pass
                else:
                    if val:
                        if sweep[i, j]:
                            sweep[i + 1, j] = val
                            for n in range(j-1, 0, -1):
                                if n in temp:
                                    temp.pop()
                                    sweep[i + 1, n] = self.pore_space
                                else:
                                    break

                        elif sweep[i + 1, j-1]:
                            sweep[i + 1, j] = self.pore_space

                        else:
                            temp.append((j))

            if np.any(sweep[i]):
                pass

            else:
                self.generate_plane()
                self.check_porosity()


    def check_porosity(self):
        porosity = 1. - (np.count_nonzero(self.matrix)/float(self.dimension * self.dimension))
        if abs(porosity - self.porosity) > self.sensitivity:
            self.generate_plane()
            self.check_porosity()
        else:
            self.check_percolation()
            self.matrix_porosity = porosity

if __name__ == "__main__":
    psphere = PSphere()
    print(psphere.matrix_porosity)
    plt.imshow(psphere.matrix, interpolation="None")
    plt.show()
