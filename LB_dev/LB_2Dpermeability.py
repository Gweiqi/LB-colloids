import numpy as np
import h5py as H
import LB_pretty as pretty
import optparse
import LB2D as LB

def HDF5_readarray(filename, data):
    ###import saved array data from LB_2Dimages####
    f = H.File(filename, "r+")
    dset = f[data][()]
    f.close()
    return dset

####initiate distribution and apply Zho he boundary to top####
def initial_distribution(q,ny,nx,rho,rhoI,vis,img,wi):
    fd = np.zeros((q,ny,nx))
    for j in range(q):
        for i in range(ny):
            fd[j,i,:] = (rho-(rhoI * (i/(ny-1.))))* wi[j]
    fd = Zho_he(fd,rho,vis)    
    ###sets solid boundaries from image data####
    #fd = set_solids(fd,img)
    return fd

def Zho_he(fd,rho,vis):
    ###function solves Zho He boundary condition against a top boundary###
    # gives the model it's initial kick to start running! Could develop a fortran version
    fd[6,0,:] = fd[2,0,:]*2./3.*rho*vis
    fd[5,0,:] = fd[1,0,:]+ 0.5*rho*vis
    fd[7,0,:] = fd[3,0,:]+ 0.5*rho*vis
    return fd

def set_solids(fs,img):
    ###sets/reinforces solid boundaries from image data####
    solids = np.zeros((q,ny,nx))
    for i in range(q):
        fs[i,img] = solids[i,img]
    return fs

def f_rho(f):
    ###function calculates density####
    rho = np.sum(f, axis=0)
    return rho

def mean_u(x):
    ###calculates mean velocity in the y and x directions####
    ###remeber to pop off ghost layers####
    ###remeber to take off -1's also####
    uy = np.average(x[0])
    ux = np.average(x[1])
    return uy,ux

def mean_rho(f, delrho):
    ###calculates a mean rho over the entire volume####
    img = np.ma.masked_where(f <= delrho, f)
    mrho = np.ma.mean(img)
    return mrho

class HDF5_write:
    def __init__(self, mrho, tau, u, f, delrho, rho, output):
        with H.File(output,"r+") as self.fi:
            print '[Writing to: %s]' % output
            self.wmrho = self.fi.create_dataset('results/mrho', data=mrho)
            self.wtau = self.fi.create_dataset('results/tau', data=tau)
            self.wuarray = self.fi.create_dataset('results/uarray', data=u)
            self.wf = self.fi.create_dataset('results/f', data=f)
            self.delr = self.fi.create_dataset('results/delr', data=delrho)
            self.rho = self.fi.create_dataset('results/rho', data=rho)
####Begin program with parser options####
parser = optparse.OptionParser()
parser.add_option('-i', '--input' , dest = 'input', help='Input a .hdf5 file')
parser.add_option('-r', '--rho', dest = 'rho', help='density array ex. 1.001,0.999',
                  default = '1.001,0.999')
parser.add_option('-t', '--tau', dest = 'tau', help='relaxation time for LGBK, default = 1',
                  default= '1')
parser.add_option('-p', '--pretty', dest = 'pretty', help='use to print LB images for animation, set interval',
                  default = None)
parser.add_option('-f', '--pfolder', dest = 'pfolder', help = 'specify a directory for images',
                  default = None)
parser.add_option('-n', '--nts', dest='maxTS', help='Enter number of time steps')
parser.add_option('-u', '--pu', dest='pu', help='Enter uy for velocity images with only y vectors',
		  default = None)
parser.add_option('-v', '--vtuple', dest='vtuple', help='enter vmin and vmax default 0.01,0.00001',
		  default = '0.01,0.00001')
parser.add_option('-q', '--notquiet', dest='q', help='enter print interval for verbose',
		  default = None)
(opts,args)=parser.parse_args()


####Open saved image data; set initial variables####
vtuple = opts.vtuple.split(',')
vmax = float(vtuple[0])
vmin = float(vtuple[1])
image = HDF5_readarray(opts.input, 'Binary_image')
rhol = opts.rho.split(',')
rho = float(rhol[0]) #add through opts.input
delr = float(rhol[0])-float(rhol[1])
rho1 = float(rhol[1])
cs = 0.577350269
cs2 = cs*cs
q = 9
ny = len(image)
nx = len(image[0])
tau = float(opts.tau)
vis = 1./3.*(tau-0.5)
maxTS = int(opts.maxTS)

# weights are consistant with up,right == positive
# down, left == negative. position 9 == 0.
wi = np.array([1./9., 1./36., 1./9., 1./36., 1./9.,
               1./36., 1./9., 1./36., 4./9.])

f = initial_distribution(q,ny,nx,rho,delr,vis,image, wi)

for i in range(maxTS):
    # call fortran subroutines to run lattice boltzmann
    rho = LB.f_rho(f, ny, nx)
    uy, ux = LB.f_u(f, rho, ny, nx)
    eu = LB.f_eu(uy, ux, ny, nx)
    usqr = LB.f_usqr(uy, ux)
    feq = LB.f_feq(eu, rho, usqr, wi, cs2, ny, nx)
    fcol = LB.f_collision(f, feq, tau, ny, nx)
    fcol = LB.f_zhohe(fcol, rho, ny, nx)
    fcol = LB.f_bounceback(f, fcol, image, ny, nx)
    f = LB.f_streaming(fcol, ny, nx)

    if opts.pretty != None:   
        if i%int(opts.pretty) == 0:
            print '[Saving image: %i]' %i
            u = [uy, ux]
            pretty.velocity_image(u, image, opts.pfolder + '/' + opts.input,i, opts.pu, vmin, vmax)

    if opts.q != None:
        if i%int(opts.q) == 0:
            print '[Iter: %i]' % i

macrho = f_rho(rho)/len(rho)
mrho = mean_rho(macrho, rho1)
u = [uy, ux]

output = HDF5_write(mrho, tau, u, f, delr, float(rhol[0]), opts.input)

print "[Done]"
