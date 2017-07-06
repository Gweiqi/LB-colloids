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
    -----------
    infile: (str) binary input file name (image file domain)
    """
    
    def __init__(self, infile):
        try:
            self.image = Image.open(infile)
        except IOError:
            raise IOError('Image type not recognized')

        self.__mode = self.image.mode
        
        if self.__mode == 'RGB':
            rgb = np.array(self.image)
            r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
            self.image = 0.2990*r + 0.5870*g + 0.1140*b
        else:
            pass
        
        self.arr = np.array(self.image)

class BoundaryCondition:
    """
    Class to instatiate open boundary layers at the top and bottom of the
    lattice boltzmann image, closed boundary layers at either side of the
    image, and binarize the array into a boolen type.
    
    data: (ndarray) NxN array of image data relating to the input file
    fluidvx: (list) list of fluid voxel grey values
    solidvx: (list) list of solid voxel grey values
    nlayers: (int) number of open boundary layers to add to the top and
                   bottom of the image array. 
    """
    def __init__(self, data, fluidvx, solidvx, nlayers):

        if isinstance(fluidvx, int):
            fluidvx = [fluidvx]

        if isinstance(solidvx, int):
            solidvx = [solidvx]
                
        self.__fluidvx = fluidvx
        self.__solidvx = solidvx
        self.__nlayers = nlayers
        self.__dimx = data.shape[1] + 2
        self.__dimy = data.shape[0] + (nlayers * 2)
        self.image = data
        self.binarized = None
        self.__set_boundary_conditions()
        print("Porosity: ", self.__porosity())
        
    def __set_boundary_conditions(self):
        setup_bc = np.zeros((self.__dimy, self.__dimx))
        setup_bc.T[0] = 1.
        setup_bc.T[-1] = 1.
        for i in range(self.__dimy):
            for j in range(self.__dimx):
                if self.__nlayers <= i <= self.__dimy - self.__nlayers - 1:
                    if 0 < j < self.__dimx - 1:
                        value = self.image[i-self.__nlayers, j-1]
                        if value in self.__fluidvx:
                            setup_bc[i, j] = 0.
                        elif value in self.__solidvx:
                            setup_bc[i, j] = 1.
                        else:
                            print 'Grey Values: ', self.grey_values
                            print 'Solid Values: ', self.solid_voxels
                            print 'Fluid Values: ', self.fluid_voxels
                            raise ValueError('Grey Value not in solid or fluid voxel values')

        self.binarized = setup_bc.astype(bool)

    def __porosity(self):
        img = self.binarized[self.__nlayers:-self.__nlayers, 1:-1]
        nsolid = np.count_nonzero(img)
        return (img.size - nsolid)/float(img.size)

    @property
    def grey_values(self):
        return np.unique(self.data)

    @property
    def solid_voxels(self):
        return self.__solidvx

    @property
    def fluid_voxels(self):
        return self.__fluidvx

    @property
    def porosity(self):
        return self.__porosity()

    @property
    def nlayers(self):
        return self.__nlayers
    

class HDF5_write:
    def __init__(self, arr, porosity, boundary, output):
        """
        Write class for LB2d_image.

        Parameters:
        ----------
        arr (np.ndarray) binarized image data
        porosity (float) porosity of the image
        boundary (int) number of boundary layers
        output (str) output hdf5 file name
        """
        self.__x = None
        with H.File(output, "w") as fi:
            print '[Writing to %s]' % output
            fi.create_dataset('Binary_image', data=arr)
            fi.create_dataset('results/porosity', data=porosity)
            fi.create_dataset('results/boundary', data=boundary)

def run(image, solid, fluid, output, boundary=5):
    """
    Object oriented approach to run the LB-Colloids script

    Parameters:
    -----------
        image: (str) image file name
        solid: (list, int) list of interger grey scale values corresponding
            to the solid phase
        fluid: (list, int) list of interger grey scale values corresponding
            to the fluid phase
        output: (str) output hdf5 file name
        boundary: (int) number of boundary layers for the model
    """
    img = Images(image)
    bc = BoundaryCondition(img, fluid, solid, boundary)
    HDF5_write(bc.binarized, bc.porosity, boundary, output)
    plt.imshow(bc.binarized, interpolation = 'nearest')
    plt.show()

    
####Test definition for later use####
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

######Begin program with parse options######
parser = optparse.OptionParser()
parser.add_option('-i','--input', dest='input', help='Input an image file')
parser.add_option('-v', '--voxels', dest='components', help='set voxel component vales ex. f0,s255 (white=0)')
parser.add_option('-o', '--output', dest='output', help='Please provide an output.hdf5 file name')
parser.add_option('-b', '--boundary', dest='boundary', help='Specify the number of top and bottom boundary layers', default='4')
parser.add_option('-c', '--config', dest='config', help='supply a lb config file')
(opts, args) = parser.parse_args()

# create an input determiner class

if __name__ == '__main__':
    if opts.input is not None:
        img = Images('Synth100_1.png')
        
        if opts.components is None:
            raise AssertionError('-v, --voxels must be supplied')
        if opts.boundary is None:
            raise AssertionError('-b, --boundary must be supplied')

        components = opts.components.split(',')
        boundary = int(opts.boundary)
        # todo: fix this issue with soil, fluid, boundary
        bc = BoundaryCondition(img.arr, solid, fluid, boundary)
        plt.imshow(bc.binarized, interpolation = 'nearest')
        plt.show()
    
    else:
        config = LBIO.Config(opts.config)

        ImageDict = {'BOUNDARY': 10}
        ModelDict = {}

        ImageDict = addIO(ImageDict, config.image_parameters())
        ModelDict = addIO(ModelDict, config.model_parameters())

        infile = ImageDict['IMAGE']

        # components = opts.components.split(',')
        fluid = ImageDict['VOID']
        solid = ImageDict['SOLID']

        nlayers = ImageDict['BOUNDARY']
        img = Images(infile)
        print '[Reading image]'
        bc = BoundaryCondition(img.arr, fluid, solid, nlayers)
        print '[Setting boundary condition]'
        print '[Porosity: %.4f]' % bc.porosity
        plt.imshow(bc.binarized, interpolation = 'nearest')
        plt.show()
        out = HDF5_write(bc.binarized, bc.porosity, nlayers, ModelDict['LBMODEL'])
        print '[Done]'
