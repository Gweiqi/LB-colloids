import numpy as np
import h5py as H
import LB_pretty as pretty
import optparse
import LB2D as LB
import LBIO
import os


def HDF5_readarray(filename, data):
    # import saved array data from LB_2Dimages
    f = H.File(filename, "r+")
    dset = f[data][()]
    f.close()
    return dset


def initial_distribution(q, ny, nx, rho, rhoI, vis, img, wi):
    # initiate distribution and apply Zho he boundary to top
    fd = np.zeros((q, ny, nx))
    for j in range(q):
        for i in range(ny):
            fd[j, i, :] = (rho - (rhoI * (i / (ny - 1.)))) * wi[j]
    return fd


def initiate_model(fd, rho, vis):
    # function solves Zho He boundary condition against a top boundary###
    # gives the model it's initial kick to start running! Could develop a fortran version
    # depreciated, not used with gravity driven flow
    fd[6, 0, :] = fd[2, 0, :] * 2. / 3. * rho * vis
    fd[5, 0, :] = fd[1, 0, :] + 0.5 * rho * vis
    fd[7, 0, :] = fd[3, 0, :] + 0.5 * rho * vis
    return fd


# def set_solids(fs, img):
    # sets/reinforces solid boundaries from image data, not used
#    solids = np.zeros((q, ny, nx))
#    for i in range(q):
#        fs[i, img] = solids[i, img]
#    return fs


def py_rho(f):
    # function calculates density
    rho = np.sum(f, axis=0)
    return rho


def py_u(f, rho):
    # calculates velocity in the y-direction and x-direction
    uy = ((f[1] + f[2] + f[3]) - (f[5] + f[6] + f[7])) / rho
    ux = ((f[3] + f[4] + f[5]) - (f[0] + f[1] + f[7])) / rho
    return uy, ux


def py_usqr(uy, ux):
    # calculates velocity squared
    usqr = uy * uy + ux * ux
    return usqr


def py_eu(uy, ux, tau, g, ny, nx):
    # calculate the einstein velocities of lattice boltzmann domain
    eu = np.zeros((9, ny, nx))
    uy = uy - (tau * g)
    eu[0, :, :] = ux
    eu[1, :, :] = ux + uy
    eu[2, :, :] = uy
    eu[3, :, :] = uy - ux
    eu[4, :, :] = -ux
    eu[5, :, :] = -ux - uy
    eu[6, :, :] = -uy
    eu[7, :, :] = ux - uy
    eu[8, :, :] = ux * 0.

    return eu


def py_feq(eu, rho, usqr, wi, csqr, ny, nx):
    # calculate the lattice boltzmann equalibrium distribution function
    feq = np.zeros((9, ny, nx))
    for i in range(9):
        feq[i, :, :] = wi[i] * rho * (1 + eu[i] / csqr + 0.5 * (eu[i] / csqr) ** 2 - usqr / (2 * csqr))
    return feq


def py_collision(f, feq, tau):
    # calculate the LB collision operator
    fcol = f - ((f - feq) / tau)
    return fcol


def py_zhohe(f, rho, ny, nx):
    # calculate zho he boundary condions, not used at present.
    fzhe = np.zeros((9, ny, nx))

    vy_lb = 1. - ((f[8] + f[0] + f[4] + 2.*(f[5] + f[6] + f[7])) / rho)
    vy_ub = -(1. - ((f[8] + f[0] + f[4] + 2.*(f[1] + f[2] + f[3])) / rho))

    # compute the unkown distributions on the upper domain
    # naming based on python zero based indexing
    f1 = (1. / 6.) * rho * vy_lb + (f[4] - f[0]) / 2. + f[5]
    f2 = (2. / 3.) * rho * vy_lb + f[6]
    f3 = (1. / 6.) * rho*vy_lb + (f[0] - f[4]) / 2. + f[7]
    
    # compute the unkown distribution on the lower domain
    f5 = -(1. / 6.) * rho * vy_ub + (f[0] - f[4]) / 2. + f[1]
    f6 = -(2. / 3.) * rho * vy_ub + f[2]
    f7 = -(1. / 6.) * rho * vy_ub + (f[4] - f[0]) / 2. + f[3]

    for j in range(ny):
        for k in range(nx):
            
            if j == 0:
                
                fzhe[0, j, k] = f[0, j, k]
                fzhe[1, j, k] = f1[j, k]
                fzhe[2, j, k] = f2[j, k]
                fzhe[3, j, k] = f3[j, k]
                fzhe[4, j, k] = f[4, j, k]
                fzhe[5, j, k] = f[5, j, k]
                fzhe[6, j, k] = f[6, j, k]
                fzhe[7, j, k] = f[7, j, k]
                fzhe[8, j, k] = f[8, j, k]

            elif j == ny - 1:
                
                fzhe[0, j, k] = f[0, j, k]
                fzhe[1, j, k] = f[1, j, k]
                fzhe[2, j, k] = f[2, j, k]
                fzhe[3, j, k] = f[3, j, k]
                fzhe[4, j, k] = f[4, j, k]
                fzhe[5, j, k] = f5[j, k]
                fzhe[6, j, k] = f6[j, k]
                fzhe[7, j, k] = f7[j, k]
                fzhe[8, j, k] = f[8, j, k]

            else:
                
                fzhe[0, j, k] = f[0, j, k]
                fzhe[1, j, k] = f[1, j, k]
                fzhe[2, j, k] = f[2, j, k]
                fzhe[3, j, k] = f[3, j, k]
                fzhe[4, j, k] = f[4, j, k]
                fzhe[5, j, k] = f[5, j, k]
                fzhe[6, j, k] = f[6, j, k]
                fzhe[7, j, k] = f[7, j, k]
                fzhe[8, j, k] = f[8, j, k]

    return fzhe


def py_bounceback(f, fcol, image, ny, nx):
    # apply bounceback conditions to LB-model collision function
    fbounce = np.zeros((9, ny, nx))
    # set bounceback indicies as local variable
    # bounce = [4, 5, 6, 7, 0, 1, 2, 3, 8]
    
    for j in range(ny):
        for k in range(nx):
            
            if image[j, k] == True:
                fbounce[0, j, k] = f[4, j, k]
                fbounce[1, j, k] = f[5, j, k]
                fbounce[2, j, k] = f[6, j, k]
                fbounce[3, j, k] = f[7, j, k]
                fbounce[4, j, k] = f[0, j, k]
                fbounce[5, j, k] = f[1, j, k]
                fbounce[6, j, k] = f[2, j, k]
                fbounce[7, j, k] = f[3, j, k]
                fbounce[8, j, k] = f[8, j, k]

            else:
                fbounce[:, j, k] = fcol[:, j, k]
    return fbounce


def py_streaming(fcol, ny, nx):
    # stream the LB-model
    fstream = np.zeros((9, ny, nx))

    for j in range(ny):
        for k in range(nx):

            if k == 0:
                kp = k + 1
                kn = nx - 1
            elif k == nx - 1:
                kp = 0
                kn = k - 1
            else:
                kp = k + 1
                kn = k - 1

            if j == 0:
                jp = j + 1
                jn = ny - 1
            elif j == ny - 1:
                jp = 0
                jn = j - 1
            else:
                jp = j + 1
                jn = j - 1

            fstream[0, j, kp] = fcol[0, j, k]
            fstream[1, jp, kp] = fcol[1, j, k]
            fstream[2, jp, k] = fcol[2, j, k]
            fstream[3, jp, kn] = fcol[3, j, k]
            fstream[4, j, kn] = fcol[4, j, k]
            fstream[5, jn, kn] = fcol[5, j, k]
            fstream[6, jn, k] = fcol[6, j, k]
            fstream[7, jn, kp] = fcol[7, j, k]
            fstream[8, j, k] = fcol[8, j, k]
    return fstream

    
def mean_u(x, img):
    # calculates mean velocity in the y and x directions
    # remeber to pop off ghost layers
    uy = np.ma.masked_where(x[0] == img, x[0])
    ux = np.ma.masked_where(x[1] == img, x[1])
    uy = np.ma.mean(uy)
    ux = np.ma.mean(ux)
    return uy, ux


def mean_rho(rho, img):
    # calculates a mean rho over the entire volume
    rho = np.ma.masked_where(rho == img, rho)
    mrho = np.ma.mean(rho)
    return mrho


def get_mean_pore_size(img, nx):
    """
    Finds the mean pore diameter of the domain

    Parameters:
        img: (np.ndarray) binary image array of domain
        nx: (int) number of pixels in x direction of fluid domain
    """
    pores = []
    for line in img:
        pore = 0
        previous_value = True
        for value in line:
            if value and previous_value:
                pass

            elif value and not previous_value:
                if pore > nx * 0.9:
                    pore = 0
                else:
                    pores.append(pore)

            elif not value and previous_value:
                pore = 1

            elif not value and not previous_value:
                pore += 1

            previous_value = value

    mean_pore_size = sum(pores) / float(len(pores))
    return mean_pore_size


def get_reynolds_number(pore_diameter, uy, porosity, rho, viscosity):
    """
    Calculate the model's reynolds number after simulation
    based off of mean velocity and mean fluid density

    Parameters:
        pore_diameter: (float) calculated lb pore diameter
        uy: (float) lb mean fluid velocity in y direction
        porosity: (float) porosity of the medium
        rho: (float) lb mean fluid density
        viscosity: (float) lb fluid viscosity
    """
    reynolds = (pore_diameter * abs(uy) * porosity * rho) / viscosity
    return reynolds


def get_velocity_conversion(reynolds_number, uy, rho,
                            pore_diameter, viscosity):
    """
    Calculates the physical velocity conversion factor from
    LB reynolds number and user supplied physical parameters.

    Parameters:
        reynolds_number: (float) Model reynolds number
        uy: (float) mean y velocity from LB model
        rho: (float) physical density of the fluid
        pore_diameter: (float) physical pore diameter
        viscosity: (float) physical fluid viscosity
    """
    phys_u = (reynolds_number * viscosity)/(rho * pore_diameter)
    factor = phys_u / abs(uy)
    return factor


def check_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)


def addIO(defaults, config):
    for key in config:
        defaults[key] = config[key]
    return defaults


class HDF5_write:
    def __init__(self, mrho, tau, u, f, rho, output,
                 mean_uy, mean_ux, pore_diameter, reynolds_number,
                 velocity_factor, img=None, porosity=None, boundary=None):
        """
        Hdf5 model write class

        Parameters:
        -----------
            mrho: (float) mean fluid density
            tau: (float) lb relaxation time
            u: (np.ndarray) fluid velocity array
            f: (np.ndarray) distribution function
            rho: (np.ndarray) density array
            output: (str) hdf file name
            mean_uy: (float) mean fluid velocity y direction
            mean_ux: (float) mean fluid velocity x direction
            pore_diameter: (float) mean pore diameter
            reynolds_number: (float) calculated reynolds number
            velocity_factor: (float) calculation velocity
                conversion factor
            img: (np.ndarray) binary image array
            porosity: (float) porosity
            boundary: (int) nlayers boundaty condition
        
        """
        self.__x = None
        print '[Writing to: %s]' % output
        try:
            with H.File(output, "r+") as fi:
                fi.create_dataset('results/mrho', data=mrho)
                fi.create_dataset('results/tau', data=tau)
                fi.create_dataset('results/uarray', data=u)
                fi.create_dataset('results/f', data=f)
                fi.create_dataset('results/rho', data=rho)
                fi.create_dataset('results/mean_uy', data=mean_uy)
                fi.create_dataset('results/mean_ux', data=mean_ux)
                fi.create_dataset('results/pore_diameter', data=pore_diameter)
                fi.create_dataset('results/reynolds_number', data=reynolds_number)
                fi.create_dataset('results/velocity_factor', data=velocity_factor)

        except:
            if os.path.isfile(output):
                os.remove(output)

            with H.File(output, "w") as fi:
                # need to add more to this to include output from LB2D_Image
                fi.create_dataset('results/mrho', data=mrho)
                fi.create_dataset('results/tau', data=tau)
                fi.create_dataset('results/uarray', data=u)
                fi.create_dataset('results/f', data=f)
                fi.create_dataset('results/rho', data=rho)
                fi.create_dataset('Binary_image', data=img)
                fi.create_dataset('results/porosity', data=porosity)
                fi.create_dataset('results/boundary', data=boundary)
                fi.create_dataset('results/mean_uy', data=mean_uy)
                fi.create_dataset('results/mean_ux', data=mean_ux)
                fi.create_dataset('results/pore_diameter', data=pore_diameter)
                fi.create_dataset('results/reynolds_number', data=reynolds_number)
                fi.create_dataset('results/velocity_factor', data=velocity_factor)


class LB2DModel(object):
    """
    object oriented method to instantiate and run a Two-Dimensional
    lattice boltzmann model. Calls upon either fortran or python
    kernals to run a model.

    Use protective programming to ensure data fits within normal model
    parameters. 

    Parameters:
    -----------
    img: (ndarray) binarized image array
    kernal: (str) the simulation kernal. Default is fortran
    
    Attributes:
    -----------
    # todo:

    Methods:
    --------
    
    run:  method to run the lb model and return a distribution function
    
    """
    def __init__(self, img, kernal='fortran'):
        self.__img = img
        self.__nlayers = None
        self.__porosity = None
        self.__kernal = kernal
        self.__gravity = 0.001
        self.__tau = 1.0
        self.__rho = 1.0
        self.__resolution = 1e-6
        self.__physical_rho = 1000.
        self.__physical_viscosity = 8.9e-4
        self.__niters = None
        self.__cs = 0.577350269
        self.__cs2 = self.cs * self.cs
        self.__ny = len(img)
        self.__nx = len(img[0])
        self.__wi = np.array([1./9., 1./36., 1./9., 1./36., 1./9.,
                              1./36., 1./9., 1./36., 4./9.])
        self.mean_rho = False
        self.mean_uy = False
        self.mean_ux = False

    def __setattr__(self, obj, value):

        if obj == 'img':
            raise NotImplementedError('Please re-instantiate LB2D to change images')

        if obj == 'kernal':
            if value.lower() not in ('python', 'fortran'):
                raise AssertionError('kernal type not recognized')
            super(LB2DModel, self).__setattr__('_LB2DModel__kernal', value.lower())

        elif obj == 'gravity':
            super(LB2DModel, self).__setattr__('_LB2DModel__gravity', float(value))

        elif obj == 'tau':
            if 0.5 > obj < 2.0:
                raise AssertionError('Tau is out of stable bounds')
            super(LB2DModel, self).__setattr__('_LB2DModel__tau', float(value))

        elif obj == 'resolution':
            super(LB2DModel, self).__setattr__('_LB2DModel__resolution', float(value))

        elif obj == 'rho':
            super(LB2DModel, self).__setattr__('_LB2DModel__rho', float(value))

        elif obj == 'physical_rho':
            super(LB2DModel, self).__setattr__('_LB2DModel__physical_rho', float(value))

        elif obj == 'physical_viscosity':
            super(LB2DModel, self).__setattr__('_LB2DModel__physical_viscosity', float(value))

        elif obj == 'niters':
            super(LB2DModel, self).__setattr__('_LB2DModel__niters', int(value))

        elif obj == 'cs':
            raise NotImplementedError('cs is implemented as constant')

        elif obj == 'cs2':
            raise NotImplementedError('cs2 is implemented as constant')

        elif obj == 'nx':
            raise NotImplementedError('nx is a constant based on image properties')

        elif obj == 'ny':
            raise NotImplementedError('ny is a constant based on image properies')
        
        else:
            super(LB2DModel, self).__setattr__(obj, value)
            
    @property
    def img(self):
        return self.__img

    @property
    def kernal(self):
        return self.__kernal

    @property
    def gravity(self):
        return self.__gravity

    @property
    def tau(self):
        return self.__tau

    @property
    def rho(self):
        return self.__rho

    @property
    def resolution(self):
        return self.__resolution

    @property
    def physical_rho(self):
        return self.__physical_rho

    @property
    def niters(self):
        return self.__niters

    @property
    def cs(self):
        return self.__cs

    @property
    def cs2(self):
        return self.__cs2

    @property
    def nx(self):
        return self.__nx

    @property
    def ny(self):
        return self.__ny

    @property
    def viscosity(self):
        return (1. / 3.) * (self.tau - 0.5)

    @property
    def physical_viscosity(self):
        return self.__physical_viscosity

    @property
    def q(self):
        return 9

    @property
    def porosity(self):
        if self.__porosity is None:
            self.__get_image_attributes()
        return self.__porosity

    @property
    def nlayers(self):
        if self.__nlayers is None:
            self.__get_image_attributes()
        return self.__nlayers

    def get_reynolds_number(self):
        """
        Returns the model's reynolds number after simulation
        """
        if not self.mean_rho or not self.mean_uy:
            print('Please run model before calculating reynolds number')
            return
        else:
            pore_diameter = self.get_mean_pore_size()
            return get_reynolds_number(pore_diameter,
                                       self.mean_uy,
                                       self.porosity,
                                       self.mean_rho,
                                       self.viscosity)

    def get_mean_pore_size(self):
        """
        Returns the mean pore diameter of the domain
        """
        return get_mean_pore_size(self.__img, self.__nx)

    def get_velocity_conversion(self):
        """
        Returns the conversion factor from LB velocity to Phys.
        """
        pore_diameter = self.get_mean_pore_size() * self.resolution
        reynolds_number = self.get_reynolds_number()
        return get_velocity_conversion(reynolds_number,
                                       self.mean_uy,
                                       self.physical_rho,
                                       pore_diameter,
                                       self.physical_viscosity)

    def __get_image_attributes(self):
        """
        Method to extract the number of boundary layers and porosity from the
        model domain
        """
        nbound = 0
        for line in self.img:
            if True in line[1:-1]:
                pass
            else:
                nbound += 1
        self.__nlayers = nbound//2

        img = self.img[self.__nlayers:-self.__nlayers, 1:-1]
        nsolid = np.count_nonzero(img)
        self.__porosity = (img.size - nsolid) / float(img.size)

    def run(self, output='LBModel.hdf5', image_int=None, image_folder=None,
            image_name="LB_", vmax=0, vmin=-0.010, verbose=None):
        """
        user method to run the lattice Boltzmann model and return the resulting
        distribution function.

       Parameters:
        -----------
            image_int (int) interval to dump velocity images to a file folder
            image_folder (str) path to folder to dump images to
            image_name (int) base name for images
            vmax (float) matplotlib vmax
            vmin (float) matplotlib vmin
            verbose (int) print interval for iterations
        """
        if self.__kernal == 'fortran':
            f = self.__run_fortran(output=output, image_int=image_int, image_folder=image_folder,
                                   image_name=image_name, vmax=vmax, vmin=vmin, verbose=verbose)
        else:
            f = self.__run_python(output=output, image_int=image_int, image_folder=image_folder,
                                  image_name=image_name, vmax=vmax, vmin=vmin, verbose=verbose)

        return f
    
    def __run_fortran(self, output='LBModel.hdf5', image_int=None, image_folder=None,
                      image_name="LB_", vmax=0, vmin=-0.010, verbose=None):
        """
        Object oriented fortran based D2Q9 LB method, uses the fortran kernal

        Parameters:
        -----------
            image_int (int) interval to dump velocity images to a file folder
            image_folder (str) path to folder to dump images to
            image_name (int) base name for images
            vmax (float) matplotlib vmax
            vmin (float) matplotlib vmin
            verbose (int) print interval for iterations
        """
        if image_int is not None:
            if image_folder is not None:
                if not os.path.exists(image_folder):
                    os.makedirs(image_folder)
                # need to add the 'xxxxx' to comply with naming from config file
                image_name = "/".join([image_folder, image_name + 'xxxxx'])
            else:
                raise AssertionError("image_folder must be supplied")

        print self.__niters
        
        f = initial_distribution(9, self.__ny, self.__nx, self.__rho , 0., self.viscosity,
                                 self.__img, self.__wi)
        for i in range(self.__niters + 1):
            rho = LB.f_rho(f, self.__ny, self.__nx)
            uy, ux = LB.f_u(f, rho, self.__ny, self.__nx)
            eu = LB.f_eu(uy, ux, self.__tau, self.__gravity, self.__ny, self.__nx)
            usqr = LB.f_usqr(uy, ux)
            feq = LB.f_feq(eu, rho, usqr, self.__wi, self.__cs2, self.__ny, self.__nx)
            fcol = LB.f_collision(f, feq, self.__tau, self.__ny, self.__nx)
            fcol = LB.f_bounceback(f, fcol, self.__img, self.__ny, self.__nx)
            f = LB.f_streaming(fcol, self.__ny, self.__nx)

            if verbose is not None and verbose:
                if i > 0:
                    if i % verbose == 0:
                        print("Iter: {:05d}".format(i))

            if image_int is not None:
                if i > 0:
                    if i % image_int == 0:
                        u = [uy[:], ux[:] * -1]
                        pretty.velocity_image(u, self.__img, image_name, i, True,
                                              vmin, vmax)

        # macrho = py_rho(rho) / len(rho)
        self.mean_rho = mean_rho(rho, self.__img)

        u = [uy[:], ux[:]* -1]

        self.mean_uy, self.mean_ux = mean_u(u, self.__img)
        pore_diameter = self.get_mean_pore_size()
        reynolds_number = self.get_reynolds_number()
        velocity_factor = self.get_velocity_conversion()

        HDF5_write(self.mean_rho, self.__tau, u, f, rho, output,
                   self.mean_uy, self.mean_ux, pore_diameter, reynolds_number,
                   velocity_factor, self.img, self.porosity, self.nlayers)

        return f

    def __run_python(self, output='LBModel.hdf5', image_int=None, image_folder=None,
                     image_name="LB_", vmax=0, vmin=-0.010, verbose=None):
        """
        Object oriented python based D2Q9 LB method, uses the python kernal
        
        Parameters:
        -----------
            image_int (int) interval to dump velocity images to a file folder
            image_folder (str) path to folder to dump images to
            image_name (int) base name for images
            vmax (float) matplotlib vmax
            vmin (float) matplotlib vmin
            verbose (int) print interval for iterations
        """
        if image_int is not None:
            if image_folder is not None:
                if not os.path.exists(image_folder):
                    os.makedirs(image_folder)
                # need to add the 'xxxxx' to comply with naming from config file
                image_name = "/".join([image_folder, image_name + 'xxxxx'])
            else:
                raise AssertionError("image_folder must be supplied")
        
        f = initial_distribution(9, self.__ny, self.__nx, self.__rho, 0., self.viscosity,
                                 self.__img, self.__wi)
        for i in range(self.__niters + 1):
            rho = py_rho(f)
            uy, ux = py_u(f, rho)
            eu = py_eu(uy, ux, self.__tau, self.__gravity, self.__ny, self.__nx)
            usqr = py_usqr(uy, ux)
            feq = py_feq(eu, rho, usqr, self.__wi, self.__cs2, self.__ny, self.__nx)
            fcol = py_collision(f, feq, self.__tau)
            fcol = py_bounceback(f, fcol, self.__img, self.__ny, self.__nx)
            f = py_streaming(fcol, self.__ny, self.__nx)

            if verbose is not None and verbose:
                if i > 0:
                    if i % verbose == 0:
                        print("Iter: {:05d}".format(i))
            
            if image_int is not None:
                if i > 0:
                    if i % image_int == 0:
                        u = [uy[:], ux[:] * -1]
                        pretty.velocity_image(u, self.__img, image_name, i, True,
                                              vmin, vmax)
        # macrho = py_rho(rho) / len(rho)
        self.mean_rho = mean_rho(rho, self.__img)

        u = [uy[:], ux[:] * -1]

        self.mean_uy, self.mean_ux = mean_u(u, self.__img)
        pore_diameter = self.get_mean_pore_size()
        reynolds_number = self.get_reynolds_number()
        velocity_factor = self.get_velocity_conversion()

        HDF5_write(self.mean_rho, self.tau, u, f, rho, output,
                   self.mean_uy, self.mean_ux, pore_diameter, reynolds_number,
                   velocity_factor, self.img, self.porosity, self.nlayers)
        return f
    
if __name__ == '__main__':
    pass

