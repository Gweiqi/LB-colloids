import Image
from scipy.ndimage import imread
import numpy as np
import h5py as H
import optparse
import matplotlib.pyplot as plt
import sys

class images:

    def __init__(self, infile):
        try:
            self.image = Image.open(infile)
        except IOError:
            sys.exit("I do not recognize this image type phil")
        self.format = self.image.format
        self.size = self.image.size
        self.mode = self.image.mode
        
        if self.mode == 'RGB':
            rgb = np.array(self.image)
            r, g, b = rgb[:,:,0], rgb[:,:,1], rgb[:,:,2]
            self.image = 0.2990*r + 0.5870*g + 0.1140*b
        else:
            pass
        
        self.imarray = np.array(self.image)

class boundary_condition:
    def __init__(self, data, fluidvx, solidvx, blayers):
        self.xlen = len(data[0])
        self.ylen = len(data)
        self.setupxbc = np.zeros((blayers, self.xlen), dtype=np.uint8)
        self.setupxbc[self.setupxbc == 0] = fluidvx
        self.tandbbc = np.concatenate((self.setupxbc, data), axis=0)
        self.tandbbc = np.concatenate((self.tandbbc, self.setupxbc), axis=0)
        self.setupybc = np.zeros((1, len(self.tandbbc)), dtype=np.uint8)
        self.setupybc[self.setupybc == 0] = solidvx
        self.ybc = np.column_stack((self.setupybc[0], self.tandbbc))
        self.bcimarray = np.column_stack((self.ybc, self.setupybc[0]))
        self.boarray = np.where(self.bcimarray == solidvx, True, False)
        self.solids =float(np.count_nonzero(self.boarray[blayers:-blayers]))
        self.all = float((self.xlen-2)*self.ylen)
        self.porosity = 1.-(self.solids/self.all)

class HDF5_write:
    def __init__(self, arr, porosity, boundary, output):
        with H.File(output,"w") as self.fi:
            print '[Writing to %s]' % output
            self.imwrite = self.fi.create_dataset('Binary_image', data=arr)
            self.pwrite = self.fi.create_dataset('results/porosity', data=porosity)
            self.bwrite = self.fi.create_dataset('results/boundary', data=boundary)
####Test definition for later use####
def HDF5_readarray(filename, data):
    f = H.File(filename, "r+")
    dset = f[data][()]
    f.close()
    return dset

######Begin program with parse options######
parser = optparse.OptionParser()
parser.add_option('-i','--input', dest='infile', help='Input an image file')
parser.add_option('-c', '--components', dest='components', help='set component vales ex. f0,s255 (white=0)')
parser.add_option('-o', '--output', dest='output', help='Please provide an output.hdf5 file name')
parser.add_option('-b', '--boundary', dest='boundary', help='Specify the number of top and bottom boundary layers', default='4')
(opts, args) = parser.parse_args()

if opts.infile == None:
    sys.exit('Oh no enter an input file!')
if opts.components== None:
    sys.exit('I don\'t know what you want to do with this image')
if opts.output == None:
    sys.exit('Wouldn\'t it be better if you saved me!')

fil = opts.infile

components = opts.components.split(',')
fluid = int(components[0][1:])
solid = int(components[1][1:])

bl = int(opts.boundary)
im = images(fil)
print '[Reading image]'
bcim = boundary_condition(im.imarray,fluid,solid,bl)
print '[Setting boundary condition]'
print '[Porosity: %.4f]' % bcim.porosity
plt.imshow(bcim.boarray, interpolation = 'nearest')
plt.show()
out = HDF5_write(bcim.boarray, bcim.porosity, opts.boundary, opts.output)
print '[Done]'
