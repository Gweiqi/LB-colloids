import numpy as np
import h5py as H
import LB_pretty as pretty
import optparse
import LB2D as LB
import LBIO
import os

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
    fd = initiate_model(fd,rho,vis)    
    ###sets solid boundaries from image data####
    #fd = set_solids(fd,img)
    return fd

def initiate_model(fd,rho,vis):
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

def py_rho(f):
    # function calculates density
    rho = np.sum(f, axis=0)
    return rho

def py_u(f, rho):
    # calculates velocity in the y-direction and x-direction
    uy = ((f[1] + f[2] + f[3]) - (f[5] + f[6] + f[7]))/rho
    ux = ((f[0] + f[1] + f[7]) - (f[3] + f[4] + f[5]))/rho
    return uy, ux

def py_usqr(uy, ux):
    # calculates velocity squared
    usqr = uy*uy + ux*ux
    return usqr

def py_eu(uy, ux, ny, nx):
    # calculate the einstein velocities of lattice boltzmann domain
    eu = np.zeros((9, ny, nx))
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
        feq[i, :, :] = wi[i]*rho*(1 + eu[i]/csqr + 0.5*(eu[i]/csqr)**2 - usqr/(2*csqr))
    return feq

def py_collision(f, feq, tau):
    # calculate the LB collision operator
    fcol = f - ((f - feq)/tau)
    return fcol

def py_zhohe(f, rho, ny, nx):
    # calculate zho he boundary condions
    fzhe = np.zeros((9, ny, nx))

    vy_lb = 1 - ((f[8] + f[0] + f[4] + 2*(f[5] + f[6] + f[7]))/rho)
    vy_ub = -(1 - ((f[8] + f[0] + f[4] + 2*(f[1] + f[2] + f[3]))/rho))
    
    # compute the unkown distributions on the upper domain
    # naming based on python zero based indexing
    f1 = 1./6.*rho*vy_lb + (f[4] - f[0])/2 + f[5]
    f2 = 2./3.*rho*vy_lb + f[6]
    f3 = 1./6.*rho*vy_lb + (f[0] - f[4])/2 + f[7]

    # compute the unkown distribution on the lower domain
    f5 = -1./6.*rho*vy_ub + (f[0] - f[4])/2 + f[1]
    f6 = -2./3.*rho*vy_ub + f[2]
    f7 = -1./6.*rho*vy_ub + (f[4] - f[0])/2 + f[3]

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
            
            if image[j,k] == True:
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

            fstream[0, j, kp]  = fcol[0, j, k]
            fstream[1, jp, kp] = fcol[1, j, k]
            fstream[2, jp, k]  = fcol[2, j, k]
            fstream[3, jp, kn] = fcol[3, j, k]
            fstream[4, j, kn]  = fcol[4, j, k]
            fstream[5, jn, kn] = fcol[5, j, k]
            fstream[6, jn, k]  = fcol[6, j, k]
            fstream[7, jn, kp] = fcol[7, j, k]
            fstream[8, j, k]   = fcol[8, j, k]
    return fstream
    
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

def check_directory(path):
    if not os.path.exists(path):
        os.makedirs(path)

def addIO(defaults, config):
    for key in config:
        defaults[key] = config[key]
    return defaults

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
parser.add_option('-c', '--config', dest = 'config', help='configuration file name')
(opts,args)=parser.parse_args()

config = LBIO.Config(opts.config)
# define defaults
ModelDict = {'KERNAL': 'fortan'}
PermeabilityDict = {'RHOT': 1.001, 'RHOB': 0.999, 'TAU': 1.0}
OutputDict = {'VMIN': 0.00001, 'VMAX': 0.01, 'VERBOSE': False, 'IMAGE_SAVE_INTERVAL': None,
              'IMAGE_SAVE_FOLDER': os.path.expanduser('~/Desktop/LBimages'), 'PLOT_Y_VELOCITY': False,
              'SAVE_IMAGE': False}

ModelDict = addIO(ModelDict, config.model_parameters())
PermeabilityDict = addIO(PermeabilityDict, config.permeability_parameters())
OutputDict = addIO(OutputDict, config.output_parameters())

####Open saved image data; set initial variables####
lbmodel = ModelDict['LBMODEL']
vmax = OutputDict['VMAX']
vmin = OutputDict['VMIN']
image = HDF5_readarray(lbmodel, 'Binary_image')

rhot = PermeabilityDict['RHOT'] #add through opts.input
rhob = PermeabilityDict['RHOB']
rhoarray = [rhot, rhob]
delr = rhot - rhob
cs = 0.577350269
cs2 = cs*cs
q = 9
ny = len(image)
nx = len(image[0])
tau = PermeabilityDict['TAU']
vis = 1./3.*(tau-0.5)
niters = PermeabilityDict['NITERS']
kernal = ModelDict['KERNAL'].lower()
if OutputDict['SAVE_IMAGE'] is True:
    check_directory(OutputDict['IMAGE_SAVE_FOLDER'])
    image_name = OutputDict['IMAGE_SAVE_FOLDER'] + '/' + lbmodel
                    
    

# weights are consistant with up,right == positive
# down, left == negative. position 9 == 0.
wi = np.array([1./9., 1./36., 1./9., 1./36., 1./9.,
               1./36., 1./9., 1./36., 4./9.])

f = initial_distribution(q, ny, nx, rhot, delr, vis, image, wi)

if kernal == 'fortran':
    for i in range(niters):
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

        if OutputDict['IMAGE_SAVE_INTERVAL'] != None:   
            if i%OutputDict['IMAGE_SAVE_INTERVAL'] == 0:
                print '[Saving image: %i]' %i
                u = [uy, ux]
                pretty.velocity_image(u, image, image_name, i, OutputDict['PLOT_Y_VELOCITY'],
                                      vmin, vmax)

        if OutputDict['VERBOSE'] is not False:
            if i%OutputDict['VERBOSE'] == 0:
                print '[Iter: %i]' % i

elif kernal == 'python':
    for i in range(niters):
        # call python subroutines to run lattice boltzmann
        rho = py_rho(f)
        uy, ux = py_u(f, rho)
        eu = py_eu(uy, ux, ny, nx)
        usqr = py_usqr(uy, ux)
        feq = py_feq(eu, rho, usqr, wi, cs2, ny, nx)
        fcol = py_collision(f, feq, tau)
        fcol = py_zhohe(fcol, rho, ny, nx)
        fcol = py_bounceback(f, fcol, image, ny, nx)
        f = py_streaming(fcol, ny, nx)

        if OutputDict['IMAGE_SAVE_INTERVAL'] != None:   
            if i%OutputDict['IMAGE_SAVE_INTERVAL'] == 0:
                print '[Saving image: %i]' %i
                u = [uy, ux]
                pretty.velocity_image(u, image, image_name, i, OutputDict['PLOT_Y_VELOCITY'],
                                      vmin, vmax)

        if OutputDict['VERBOSE'] is not False:
            if i%OutputDict['VERBOSE'] == 0:
                print '[Iter: %i]' % i

else:
    raise Exception('Kernal type not supported')
    
macrho = py_rho(rho)/len(rho)
mrho = mean_rho(macrho, rhob)
u = [uy, ux]

output = HDF5_write(mrho, tau, u, f, delr, rhot, lbmodel)

print "[Done]"
