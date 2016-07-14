import numpy as np
import h5py as H
import optparse
import sys
#need to write self.delr in LB_dev
class HDF5_reader:
    def __init__(self,HDF_name):
        hdf = H.File(HDF_name, 'r+')
        self.imarray = hdf['Binary_image'][()]
        self.mrho = hdf['results/mrho'][()]
        self.porosity = hdf['results/porosity'][()]
        self.tau = hdf['results/tau'][()]
        self.uarray = hdf['results/uarray'][()]
        self.yu = hdf['results/uarray'][()][0]
        self.bound = int(hdf['results/boundary'][()])
        self.delr = hdf['results/delr'][()]
        self.rho = hdf['results/rho'][()]
        
def mean_yvel(yvel, image, bound):
    image = np.invert(image)
    uy = np.array([yvel[i]*image[i] for i in range(len(yvel))])
    uy = np.ma.masked_where(uy == 0, uy)
    uy = uy[bound:-bound]
    mean =uy.mean()
    return abs(mean)

def rho_grad(urho, delr, bound, image):
    rho_top = (urho-(delr * (bound/(len(image)-1.))))
    rho_bottom = (urho-(delr * ((len(image)-bound)/(len(image)-1.))))
    delta = rho_top - rho_bottom
    return delta

class HDF5_writer:
    def __init__(self, k, output):
	with H.File(output,"r+") as self.fi:
	    print '[Writing to: %s]' % output
	    self.wk = self.fi.create_dataset('results/permeability', data = k)

parser = optparse.OptionParser()
parser.add_option('-i', '--input', dest='input', help= 'input a LB .hdf5 file',
                  default = None)
parser.add_option('-r', '--resolution', dest='res', help= 'input image resolution',
                  default = None)
parser.add_option('-k', '--savek', dest='savek', help='-k y to save k to hdf5 file',
                  default = None)
(opts,args) = parser.parse_args()

if opts.input == None or opts.res == None:
    sys.exit('\n[Use python LB_results -h for help menu]\n')

results = HDF5_reader(opts.input)
yvel = mean_yvel(results.yu, results.imarray, results.bound)
visc = 1./3.*(results.tau - 0.5)
NL = len(results.imarray)-(results.bound*2)
delr = rho_grad(results.rho, results.delr, results.bound, results.imarray)

permeability = ((results.mrho*visc*yvel)/(delr*(1./3.)*NL))* float(opts.res)**2
#look into how to do perm. and k for 2d LB
pum = permeability*float(opts.res) #permeability**2/10**-12 
ksat = ((pum * 10**-12 * 9.81 * 1000.) / (8.94 * 10**-4)) * 100

print '\n[uy        = %.3f]' % yvel
print '[viscosity = %.3f]'   % visc
print '[NL        = %i  ]'   % NL
print '[rho       = %.3f]'   % results.mrho
print '[Delta rho = %.3f]'   % results.delr
print '[cs^2      = 0.333]'
print '[resolution= %.3f]'   % float(opts.res)
print '[k         = %.5f um^2]' % pum
print '[Ksat      = %.9f cm/s]' % ksat

if opts.savek != None:
    HDF5_writer(pum, opts.input)
print '[Done!]\n'
