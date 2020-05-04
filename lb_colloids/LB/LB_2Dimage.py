"""
Lattice boltzmann image preparation and domain setup module

This module contains an image reading utility and a binarization
utility. The LB_2Dimage module is aliased in LB-Colloids as LBImage.

Example usage of how to create a binary image with boundary conditions applied
and save the base model for usage in with LB_permeability

>>> from lb_colloids import LBImage
>>> image = LBImage.Images("my_thin_section.png")
>>> binary = LBImage.BoundaryCondition(image, fluidvx=[0], solidvx=[233, 255], nlayers=3)
>>> binary.binarized
>>> HDF5_write(binary.binarized, binary.porosity, binary.boundary, 'LBModel.hdf5')
"""
import Image
from scipy.ndimage import imread
import numpy as np
import h5py as H
import optparse
import matplotlib.pyplot as plt
import sys
import LBIO


class Images:
    """
    Convience class to open and adjust RGB images to BW if
    necessary.

    Parameters:
    ----------
    :param string infile:
        binary input file name (image file domain)
    """
    
    def __init__(self, infile):
        try:
            self.image = Image.open(infile)
        except IOError:
            raise IOError('Image type not recognized')

        self.__mode = self.image.mode
        
        if self.__mode in ('RGB', 'RGBA'):
            rgb = np.array(self.image)
            r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
            self.image = 0.2990*r + 0.5870*g + 0.1140*b
        else:
            pass
        
        self.arr = np.array(self.image, dtype=int)


class BoundaryCondition:
    """
    Class to instatiate open boundary layers at the top and bottom of the
    lattice boltzmann image, closed boundary layers at either side of the
    image, and binarize the array into a boolen type.

    Parameters:
    ----------
    :param np.ndarray data:
        NxN array of image data relating to the input file
    :param list fluidvx:
        list of fluid voxel grey values
    :param list solidvx:
        list of solid voxel grey values
    :param int nlayers:
        number of open boundary layers to add to the top and
        bottom of the image array.
    """
    def __init__(self, data, fluidvx, solidvx, nlayers, bottom=False):

        if isinstance(fluidvx, int):
            fluidvx = [fluidvx]

        if isinstance(solidvx, int):
            solidvx = [solidvx]
                
        self.__fluidvx = fluidvx
        self.__solidvx = solidvx
        self.__nlayers = nlayers
        self.__dimx = data.shape[1] + 2
        self.__bottom = bottom
        if bottom:
            self.__lower_layer = 0
            self.__dimy = data.shape[0] + nlayers
        else:
            self.__dimy = data.shape[0] + (nlayers * 2)
        self.image = data
        self.__binarized = None
        self.__set_boundary_conditions()
        print("Porosity: ", self.__porosity())
        
    def __set_boundary_conditions(self):
        if self.__bottom:
            bottom = self.__dimy - 1
        else:
            bottom = self.__dimy - self.__nlayers - 1
        setup_bc = np.zeros((self.__dimy, self.__dimx))
        setup_bc.T[0] = 1.
        setup_bc.T[-1] = 1.
        for i in range(self.__dimy):
            for j in range(self.__dimx):
                if self.__nlayers <= i <= bottom:
                    if 0 < j < self.__dimx - 1:
                        value = self.image[i-self.__nlayers, j-1]
                        if value in self.__fluidvx:
                            setup_bc[i, j] = 0.
                        elif value in self.__solidvx:
                            setup_bc[i, j] = 1.
                        else:
                            print('Grey Values: ', self.grey_values)
                            print('Solid Values: ', self.solid_voxels)
                            print('Fluid Values: ', self.fluid_voxels)
                            raise ValueError('Grey Value not in solid or fluid voxel values')

        self.__binarized = setup_bc.astype(bool)

    def __porosity(self):
        img = self.binarized[self.__nlayers:-self.__nlayers, 1:-1]
        nsolid = np.count_nonzero(img)
        return (img.size - nsolid)/float(img.size)

    @property
    def binarized(self):
        """
        :return: Binary image with boundary conditions applied
        """
        return self.__binarized

    @property
    def grey_values(self):
        """
        :return: unique image grayscale values
        """
        return np.unique(self.image)

    @property
    def solid_voxels(self):
        """
        :return: user defined solid voxels
        """
        return self.__solidvx

    @property
    def fluid_voxels(self):
        """
        :return: user defined fluid volxels
        """
        return self.__fluidvx

    @property
    def porosity(self):
        """
        :return: image porosity
        """
        return self.__porosity()

    @property
    def nlayers(self):
        """
        :return: Number of boundary condition layers
        """
        return self.__nlayers
    

class HDF5_write:
    """
    Write class for LB2d_image. Writes a HDF5 file that includes
    the binary image, porosit &, number of boundary layers,

    Parameters:
    ----------
    :param np.ndarray arr:
        binarized image data
    :param float porosity:
        porosity of the image
    :param int boundary:
        number of boundary layers
    :param str output:
        output hdf5 file name
    """
    def __init__(self, arr, porosity, boundary, output):

        self.__x = None
        with H.File(output, "w") as fi:
            print('[Writing to {}]'.format(output))
            fi.create_dataset('Binary_image', data=arr)
            fi.create_dataset('results/porosity', data=porosity)
            fi.create_dataset('results/boundary', data=boundary)


def run(image, solid, fluid, output, boundary=5):
    """
    Funcitonal approach to run the LB-Colloids script

    Parameters:
    ----------
    :param str image:
        image file name
    :param list solid:
        list of interger grey scale values corresponding
        to the solid phase
    :param list fluid:
        list of interger grey scale values corresponding
        to the fluid phase
    :param str output:
        output hdf5 file name
    :param int boundary:
        number of boundary layers for the model

    >>> LBImage.run('my_thin_section.png', fluid=[0], solid=[233, 255],
    >>>             output="LBModel.hdf5", nlayers=3)
    >>>
    """
    img = Images(image)
    bc = BoundaryCondition(img, fluid, solid, boundary)
    HDF5_write(bc.binarized, bc.porosity, boundary, output)
    plt.imshow(bc.binarized, interpolation = 'nearest')
    plt.show()

    
# Test definition for later use
def HDF5_readarray(filename, data):
    f = H.File(filename, "r+")
    dset = f[data][()]
    f.close()
    return dset


# Facilitates defaults through kwargs
def addIO(defaults, config):
    for key in config:
        defaults[key] = config[key]
    return defaults


def prep_conponents(components):
    components = components.split(',')
    solid = [int(i[1:]) for i in components if 's' in i]
    fluid = [int(i[1:]) for i in components if 'f' in i]
    return solid, fluid

if __name__ == '__main__':
    pass

