"""
Colloid_Setup contains background classes and methods to prepare a
model domain for colloid simulation. The user should not need to
import or call methods directly from this module.
"""

import numpy as np
from scipy import interpolate
from copy import copy
import h5py as H
import sys


class Hdf5Reader(object):
    """
    Hdf5 reader class to grab results from lattice Boltzmann model runs
    to parameterize the colloid simulation. Consider moving this class to
    the Colloid_IO module

    Parameters:
    ----------
    :param str HDF_name: lattice boltzmann hdf file name.

    :ivar np.ndarray imarray: binary image array defining model boundaries
    :ivar np.ndarray uarry: velocity array of [y, x]
    :ivar np.ndarray yu: y velocity array from lattice Boltzmann simulation
    :ivar np.ndarray xu: x velocity array from lattice Boltzmann simulation
    :ivar float mean_yu: mean velocity in the y direction
    :ivar float mean_xu: mean velocity in the x direction
    :ivar float velocity_factor: velocity dimensionalization factor
    """
    def __init__(self, hdf_name):
        hdf = H.File(hdf_name, 'r+')
        self.imarray = hdf['Binary_image'][()]
        self.uarray = hdf['results/uarray'][()]
        self.yu = hdf['results/uarray'][()][0]
        self.xu = hdf['results/uarray'][()][1]
        self.mean_yu = hdf['results/mean_ux'][()]
        self.mean_xu = hdf['results/mean_uy'][()]
        self.velocity_factor = hdf['results/velocity_factor'][()]
        hdf.close()


class GridArray(object):
    """
    Gridarray class creates arrays of distances from pore spaces, corrects for
    interpolation effects at pore boundaries, and creates vector arrays that are
    later used to give direction to forces.

    Parameters:
    ----------
    :param (np.array, bool) arr: Array of boolean porous media (segmented image array)
    :param float gridres: model resolution in meters
    :param int gridsplit: interpolation factor for refining grid mesh
    :param bool solid: solid phase boolean identifier, default=True
    """
    def __init__(self, arr, gridres, gridsplit, solid=True):
        self.yarr = np.copy(arr.T)
        self._vimgx, self._vimgy = self._create_vector_array(arr, solid)
        self.__gridx, self.__vector_x = self._arrx(arr, self._vimgx, gridres, gridsplit)
        self.__gridy, self.__vector_y = self._arry(self.yarr, self._vimgy, gridres, gridsplit)

    @property
    def gridx(self):
        """
        :return: (np.array, np.float) Array of distances from nearest solid phase in the x-direction
        """
        return copy(self.__gridx)

    @property
    def gridy(self):
        """
        :return: (np.array, np.float) Array of distances from nearest solid phase in the y-direction
        """
        return copy(self.__gridy)

    @property
    def vector_x(self):
        """
        :return: (np.array, np.float) Array of specific vector directions in the x-direction (-1 == left, 1 == right)
        """
        return copy(self.__vector_x)

    @property
    def vector_y(self):
        """
        :return: (np.array, np.float) Array of specific vector directions in the y-direction (-1 == down, 1 == up)
        """
        return copy(self.__vector_y)

    def _create_vector_array(self, img, solid):
        """
        creates an x-dir copy and y-dir copy of the image domain for vector
        directions to populate.
        """
        vimgx = np.copy(img)
        vimgy = np.copy(img.T)
        vimgx[vimgx == solid] = np.nan
        vimgy[vimgy == solid] = np.nan
        return vimgx, vimgy

    def _arrx(self, arr, vector_x, gridres, gridsplit):
        """
        Method that handles looping through array and sends it to distance_gridx
        """
        for line in range(len(arr)):
            arr[line], vector_x[line] = self._distance_gridx(arr[line], vector_x[line], gridres, gridsplit)
        return arr, vector_x

    def _arry(self, yarr, vector_y, gridres, gridsplit):
        """
        Method that handles looping through array and sends it to distance_gridy
        """
        for line in range(len(yarr)):
            ylen = len(yarr[0])
            yarr[line], vector_y[line] = self._distance_gridy(yarr[line], vector_y[line], gridres, ylen, gridsplit)
        return yarr.T, vector_y.T
    
    def _distance_gridx(self, line, vline, gridres, gridsplit):
        """
        Method uses linear interpolation to correct pore boundary space

        Follows by creating an array of pore boundaries and then counts the distance
        from the nearest solid. It also creates a a vector direction array that corresponds.
        """

        # Saner(?) method which does not raise warnings from Python 2.7.12
        # Forward compatable with Python 3.5.2
        for i in range(len(line)):
            if line[i] > 0.:
                vline[i] = 1
                line[i] = 1
            else:
                vline[i] = 0
                line[i] = 0
        np.array(line)
        
        # create an array of pore boundaries from binary system
        boundary= np.where(np.abs(np.diff(line)) >= 1)[0]
        try:
            for i in range(1, len(boundary), 2):
                rbound = boundary[i] + 1
                lbound = boundary[i-1]
                gap = rbound - lbound

                if gap % 2 == 0:
                    gap = gap//2

                    if gridres > 1e-9:
                        left = (np.arange(1, gap + 1) * gridres) - (gridres - 1e-9)
                    else:
                        left = np.arange(1, gap + 1) * gridres

                    right = left[::-1]
                    line[lbound:rbound] = np.append(left, right)

                    left = np.ones(gap) * -1
                    right = np.ones(gap)
                    vline[lbound:rbound] = np.append(left, right)

                else:
                    gap = gap//2

                    if gridres > 1e-9:
                        left = (np.arange(1, gap + 2) * gridres) - (gridres - 1e-9)
                        right = (np.arange(1, gap + 1) * gridres) - (gridres - 1e-9)
                    else:
                        left = np.arange(1, gap + 2) * gridres
                        right = np.arange(1, gap + 1) * gridres

                    line[lbound:rbound] = np.append(left, right)

                    left = np.ones(gap + 1) * -1
                    right = np.ones(gap)
                    vline[lbound:rbound] = np.append(left, right)
                    
        except IndexError:
            print 'volume does not percolate'
            sys.exit()
        return line, vline

    def _distance_gridy(self, line, vline, gridres, ylen, gridsplit):
        '''
        Method uses linear interpolation to correct pore boundary space

        Follows by creating an array of pore boundaries and then counts the distance
        from the nearest solid. It also creates a a vector direction array that corresponds.
        '''
        
        for i in range(len(line)):
            if line[i] > 0:
                vline[i] = 1
                line[i] = 1
            else:
                vline[i] = 0
                line[i] = 0
        np.array(line)
        
        # create an array of pore boundaries from binary system
        boundary = np.where(np.abs(np.diff(line)) >= 1)[0]

        if len(boundary) > 0:
            rbound = boundary[0] + 1 
            lbound = 0

            if gridres > 1e-9:
                # Use this statement to enforce nm scale DLVO @ boundary
                top = (np.arange(rbound, lbound, -1) * gridres) - (gridres - 1e-9)
            else:
                top = np.arange(rbound, lbound, -1) * gridres

            line[lbound:rbound] = top
            vtop = np.ones(len(top)) * -1
            vline[lbound:rbound] = vtop
            
            for i in range(2, len(boundary), 2):
                rbound = boundary[i] + 1
                lbound = boundary[i - 1] + 1
                gap = rbound - lbound
                if gap % 2 == 0:
                    gap = gap // 2

                    if gridres > 1e-9:
                        # use this statement to enforce nm scale DLVO @ boundaries
                        left = (np.arange(1, gap + 1) * gridres) - (gridres - 1e-9)
                    else:
                        left = np.arange(1, gap + 1) * gridres

                    right = left[::-1]
                    line[lbound:rbound] = np.append(left, right)

                    left = np.ones(gap)
                    right = np.ones(gap) * -1
                    vline[lbound:rbound] = np.append(left, right)
                    
                else:
                    gap = gap // 2

                    if gridres > 1e-9:
                        left = (np.arange(1, gap + 2) * gridres) - (gridres - 1e-9)
                        right = (np.arange(1, gap + 1) * gridres) - (gridres - 1e-9)
                    else:
                        left = np.arange(1, gap + 2) * gridres
                        right = np.arange(1, gap + 1) * gridres

                    line[lbound:rbound] = np.append(left, right)

                    left = np.ones(gap + 1)
                    right = np.ones(gap) * -1
                    vline[lbound:rbound] = np.append(left, right)
                    
            rbound = ylen
            lbound = boundary[-1] + 1
            gap = rbound - lbound

            if gridres > 1e-9:
                bottom = (np.arange(1, gap + 1) * gridres) - (gridres - 1e-9)
            else:
                bottom = np.arange(1, gap + 1) * gridres

            line[lbound:rbound] = bottom

            bottom = np.ones(gap)
            vline[lbound:rbound] = bottom
        else:
            pass
        return line, vline


def LBVArray(LBv, img):
    """
    Method to create a velocity array for use in the colloid simulation model

    Parameters:
    ----------
    :param np.ndarray LBv: lattice boltzmann velocity array
    :param np.ndarray img: boolean image array

    Returns:
    -------
    :return: np.ndarray of velocity
    """
    vel = np.zeros((len(LBv), len(LBv[0])))
    invert = np.invert(img.astype(bool))
    vel = np.array([LBv[i] * invert[i] for i in range(len(LBv))])
    print(vel.shape)
    return vel


def InterpV(LBv, gridsplit, img=False):
    """
    Interpolation method for the lattice Boltzmann velocity array

    Parameters:
    ----------
    :param np.ndarray LBv: lattice boltzmann velocity array with pore boundaries enforced
    :param float gridsplit: interpolation factor
    :param bool img: flag to indicate boolean image interpolation or velocity interpolation

    Returns:
    -------
    :return: (np.ndarray) interpolated velocity array
    """
    ylen = len(LBv)
    xlen = len(LBv[0])
    xindex = np.arange(0, xlen)
    yindex = np.arange(0, ylen)
    ifactor = 1. / gridsplit
    xnew = np.arange(0, xlen - 1 + (ifactor), ifactor)
    ynew = np.arange(0, ylen - 1 + (ifactor), ifactor)
    f = interpolate.interp2d(xindex, yindex, LBv, kind='linear')
    znew = f(xnew, ynew)

    if img:
        # correct pore boundaries from interpolation
        znew[znew >= ifactor * gridsplit / 2.] = 1  # 0] = 1# ifactor*gridsplit/2.] = 1
        znew[znew < ifactor * gridsplit / 2.] = 0  # 0] = 0# ifactor*gridsplit/2.] = 0

    return znew



