"""
The Pshpere module contains a class named PShpere which allows the user to
generate synthetic porous media, and to get information about that porous
media

A user can instantiate and use the PSphere() object as follows:

>>> from lb_colloids import PSphere
>>> img = Psphere(dimension=200, radius=20, porosity=0.375, sensitivity=0.01)
>>> # hydraulic radius can be calculated
>>> rh = img.calculate_hydraulic_radius(resolution=1e-06)
>>> # to get a copy of the porous media use
>>> matrix = img.matrix
>>> # save the image
>>> img.save("test_image.png")

"""

import numpy as np
import random
import math
from PIL import Image
import matplotlib.pyplot as plt


class PSphere(object):
    """
    Pshpere is a class that allows for the automated generation of
    synthetic porous media in two-dimensions. This approach can be
    expanded to three dimensions with some effort.


    Parameters:
    ----------
    :param int radius: grain size radius
    :param float porosity: target porosity for porous media
    :param int dimension: the x and y dimension in pixels for the domain
    :param float sensitivity: a porosity sensitivity target. This is the allowable range of error for PShpere
    """
    def __init__(self, radius=20, porosity=0.5, dimension=256, sensitivity=0.08):
        self.radius = radius
        self.porosity = porosity
        self.sensitivity = sensitivity
        self.dimension = dimension
        self.matrix = np.ones((dimension, dimension), dtype=bool)
        self.matrix_porosity = 0.
        self.matrix_rh = 0.
        self.particle_space = False
        self.pore_space = True
        self.percolates = False
        good = False

        while not good:
            self.generate_plane()
            self.check_percolation()
            self.check_porosity()
            print self.matrix_porosity
            if abs(self.matrix_porosity - self.porosity) <= sensitivity:
                if self.percolates:
                    good = True

            else:
                print("Regenerating porous media")
                self.matrix = np.ones((dimension, dimension), dtype=bool)
                # self.percolates = False

    def generate_plane(self):
        """
        Main method used to generate a porous media plane by PSphere,
        this should not be called by the user
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
                self.patch(x, y, slice_radius)

    def patch(self, x, y, radius):
        """
        The patch method is used to set grains into a porous media.
        Not to be called by the user!

        Parameters:
        ----------
        :param int x: x index location
        :param int y: y index location
        :param int radius: grain radius
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

    def iround(self, val):
        """
        Rounding routine to set index locations.

        :param float val: floating point value
        :return int:
        """
        if val - math.floor(val) >= 0.5:
            return int(math.ceil(val))
        else:
            return int(math.floor(val))

    def check_percolation(self):
        """
        Modified sweep line technique that checks the porous media percolation
        Pshpere automatically calls this during porous media creation.
        """
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
                self.percolates = True

            else:
                self.percolates = False
                return
                # self.generate_plane()
                # self.check_percolation()

    def calculate_hydraulic_radius(self, resolution):
        """
        Calculates the hydraulic radius of the porous medium

        Parameters:
        ----------
        :param float resolution: model resolution applied to image

        :return: hydraulic radius of the image
        """
        surface = 0
        for i, xdim in enumerate(self.matrix):
            for j in range(1, len(xdim)):
                if (xdim[j-1], xdim[j]) == (True, False) or\
                        (xdim[j-1], xdim[j]) == (False, True):
                    surface += 1

        for i, ydim in enumerate(self.matrix.T):
            for j in range(1, len(ydim)):
                if (ydim[j-1], ydim[j]) == (True, False) or\
                        (ydim[j-1], ydim[j]) == (False, True):
                    surface += 1

        pore = np.count_nonzero(self.matrix) * (resolution ** 3)

        sa0 = surface * (resolution ** 2)

        return pore / sa0

    @staticmethod
    def static_hydraulic_radius(matrix, invert=True):
        """
        Static method to calculate the hydraulic radius of a given porous medium

        Parameters:
        ----------
        :param np.ndarray matrix: boolean array corresponding to porous media
        :param bool invert: inverts model, pore space needs to be set to True

        :return: hydraulic radius of the image
        """

        if invert:
            matrix = np.invert(matrix)

        surface = 0
        for i, xdim in enumerate(matrix):
            for j in range(1, len(xdim)):
                if (xdim[j - 1], xdim[j]) == (True, False) or \
                                (xdim[j - 1], xdim[j]) == (False, True):
                    surface += 1

        for i, ydim in enumerate(matrix.T):
            for j in range(1, len(ydim)):
                if (ydim[j - 1], ydim[j]) == (True, False) or \
                                (ydim[j - 1], ydim[j]) == (False, True):
                    surface += 1

        pore = np.count_nonzero(matrix)

        sa0 = surface

        return pore / float(sa0)

    @staticmethod
    def static_surface_area(matrix, invert=True):
        """
        Static method to calculate the non-dimensional surface
        area of a given porous medium

        Parameters:
        ----------
        :param np.ndarray matrix: boolean array corresponding to porous media
        :param bool invert: inverts model, pore space needs to be set to True

        :return: hydraulic radius of the image
        """

        if invert:
            matrix = np.invert(matrix)

        surface = 0
        for i, xdim in enumerate(matrix):
            for j in range(1, len(xdim)):
                if (xdim[j - 1], xdim[j]) == (True, False) or \
                        (xdim[j - 1], xdim[j]) == (False, True):
                    surface += 1

        for i, ydim in enumerate(matrix.T):
            for j in range(1, len(ydim)):
                if (ydim[j - 1], ydim[j]) == (True, False) or \
                        (ydim[j - 1], ydim[j]) == (False, True):
                    surface += 1

        return surface

    @staticmethod
    def static_porosity(matrix, invert=True):
        """
        Static method to calculate the porosity of a given porous media

        Parameters:
        ----------
        :param np.ndarray matrix: boolean array corresponding to porous media
        :param bool invert: inverts model, pore space needs to be set to True

        :return: porosity of the image
        """
        if invert:
            matrix = np.invert(matrix)

        porosity = np.count_nonzero(matrix)/float(matrix.shape[0] * matrix.shape[1])
        return porosity

    def check_porosity(self):
        porosity = (np.count_nonzero(self.matrix)/float(self.dimension * self.dimension))
        if abs(porosity - self.porosity) > self.sensitivity:
            print('Warning: medium porosity outside sensitivity')

        self.matrix_porosity = porosity

    def save_image(self, image_name):
        """
        Save method, to save an image to file!

        :param str image_name: image path and name
        """
        matrix = np.invert(self.matrix)
        plt.imsave(image_name, matrix, cmap="gray")
        plt.close()


if __name__ == "__main__":
    psphere = PSphere(dimension=200, radius=18, porosity=0.375, sensitivity=0.02)
    print(psphere.matrix_porosity)
    print psphere.calculate_hydraulic_radius(1e-06)
    plt.imshow(psphere.matrix, interpolation="None")
    plt.show()
