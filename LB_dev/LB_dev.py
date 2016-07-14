import numpy as np
import h5py as H
import LB_pretty as pretty
import optparse

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
    fd[4,0,:] = fd[2,0,:]*2./3.*rho*vis
    fd[7,0,:] = fd[5,0,:]+ 0.5*rho*vis
    fd[8,0,:] = fd[6,0,:]+ 0.5*rho*vis
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

def f_u(f,ev,rho):
    ###function calculates macroscopic velocity for each vx####
    ###u[0] is array of vy; u[1] is array of vx####
    u = np.dot(ev,f.transpose((1,0,2)))/rho
    return u

def f_eu(u,ev):
    ###calculates e*u for feq####
    ###Note transpose axes are 1=y, 0=z, 2=x####
    eu = np.dot(ev.transpose(),u.transpose((1,0,2)))
    return eu

def f_usqr(u):
    ###function creates uqsr for feq####
    ###returns nx by ny array####
    usqr = (u[0]*u[0]+u[1]*u[1])
    return usqr

def f_feq(wi,rho,eu,usqr,q,ny,nx):
    ###calculates feq####
    feq = np.zeros((q,ny,nx))
    for i in range(q):
        feq[i,:,:] = wi[i]*rho*(1.+(3.*eu[i])+ (4.5*eu[i]*eu[i]) - (3./2.*usqr))
    return feq

def f_collison(tau,f,feq):
    ###performs the collision step####
    fcol = f-((f-feq)/tau)
    return fcol

def f_bounceback(fcol,f,bounce,img,q):
    ###facilitates the half-bounceback step####
    ###creates no-slip boundaries####
    for i in range(q):
        fcol[i,img]= f[bounce[i],img]
    return fcol

def f_streaming(f,fcol,ev,q):
    ###propogates the particles according to eigenvectors####
    #fs = np.array([])
    for i in range(q):
        f[i,:,:] = np.roll(np.roll(fcol[i,:,:],ev[0,i],axis=0),ev[1,i],axis=1)
    return f

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
####Create eigenvectors and assign weights####


ev = np.array([[0,0,-1,0,1,-1,-1,1,1],
              [0,1,0,-1,0,1,-1,-1,1]])
wi = np.array([4./9.,1./9.,1./9.,1./9.,1./9.,
               1./36.,1./36.,1./36.,1./36.])

bounce = np.array([0,3,4,1,2,7,8,5,6])

i1 = np.array([2,5,6]) #Unkown on bottom 'wall'
i2 = np.array([0,1,3]) #horizontal middle
i3 = np.array([4,7,8]) #Unknown on top 'wall'

f = initial_distribution(q,ny,nx,rho,delr,vis,image, wi)

for i in range(maxTS):
    f[i1,-1,:]=f[i1,-2,:]
    rho = f_rho(f)
    macrho = f_rho(rho)/len(rho)
    mrho = mean_rho(macrho, rho1)
    u = f_u(f,ev,rho)
    eu = f_eu(u,ev)
    usqr = f_usqr(u)
    feq = f_feq(wi,rho,eu,usqr,q, ny,nx)
    rho = f_rho(rho)/len(rho)
    f[i3,0,:] = f[i1,0,:]+feq[i3,0,:]-f[i1,0,:]
    fcol = f_collison(tau, f, feq)
    fcol = f_bounceback(fcol,f, bounce, image, q)
    f = f_streaming(f, fcol, ev, q) 

    if opts.pretty != None:   
        if i%int(opts.pretty) == 0:
            print 'saving image %i' %i
            pretty.velocity_image(u, image, opts.pfolder + '/' + opts.input,i, opts.pu, vmin, vmax)

    if opts.q != None:
	if i%int(opts.q) == 0:
            print '[Iter: %i]' % i

output = HDF5_write(mrho, tau, u, f, delr, float(rhol[0]), opts.input)

print "[Done]"
