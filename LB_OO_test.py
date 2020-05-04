from lb_colloids import LBImage
from lb_colloids import LB2DModel
import matplotlib.pyplot as plt
from lb_colloids import cIO
from lb_colloids import ColloidModel
import os
import time

from lb_colloids import ColloidMath as cm
import numpy as np
import matplotlib

path = os.path.dirname(os.path.realpath(__file__))
lbname = os.path.join(path, "LBModel.hdf5")

imgname = os.path.join(path, "Synth100_2.png")


nlayers = 5
solidvx = 0
fluidvx = 253

img = LBImage.Images(imgname)
binary = LBImage.BoundaryCondition(img.arr, fluidvx, solidvx, nlayers)
plt.imshow(binary.binarized, interpolation='nearest')
plt.show()

lbmodel = LB2DModel(binary.binarized)
lbmodel.niters = 1000
lbmodel.run(output=lbname, verbose=10, image_int=10, image_folder='test')
x = lbmodel.get_velocity_conversion()

# now we work with LB_colloids.py
# need to think about namespace now.
# now we can use this to do other analysis such as sensitivity analysis
# sensitivity analysis, do this soonish!
t0 = time.time()
tc0 = time.clock()
io = cIO.ColloidsConfig()
print(io.valid_model_parameters)
io['lbmodel'] = lbname
io['ncols'] = 500
io['iters'] = 200000
io['lbres'] = 1e-6
io['gridref'] = 10
io['ac'] = 1e-6
io['timestep'] = 5e-7
io['temperature'] = 298.
# io['multiple_config'] = False
io["col_col_update"] = 10
io['continuous'] = 0

io['i'] = 0.001

io['print_time'] = 10000
io['plot'] = False
io['showfig'] = False
io['endpoint'] = os.path.join(path, 'Synth100_3.endpoint')
io['timeseries'] = os.path.join(path, 'Synth100_3.pathline')
io['store_time'] = 100

print(io.model_parameters)
print(io.chemical_parameters)
print(io.physical_parameters)
print(io.output_control_parameters)

# config = os.path.join(path, 'Synth100_3.config')

config = cIO.Config(io.config)
ColloidModel.run(config)
t1 = time.time()
tc1 = time.clock()
print(t1 - t0)
print(tc1 - tc0)

# import matplotlib.pyplot as plt
"""
d = {}

for key, value in io.chemical_parameters.items():
    d[key.lower()] = value
for key, value in io.model_parameters.items():
    d[key.lower()] = value

print(d)

arr = np.linspace(1e-9, 1e-8, 100)
xarr = np.array([arr])
varr = np.array([np.ones(100)])

i = [0.001]


for j in i:
    d['I'] = j
    dlvo = cm.DLVO(xarr, xarr, xvArr=varr, yvArr=varr, **d)
    plt.plot(np.abs(arr), dlvo.EDLx[0], label='{} M'.format(j), lw=1.5)

plt.ylabel('Force (N)', fontsize=14)
plt.xlabel('X (m)', fontsize=14)
plt.legend(loc=0, fontsize=12)
plt.show()

for j in i:
    d['I'] = j
    dlvo = cm.DLVO(xarr, xarr, xvArr=varr, yvArr=varr, **d)
    plt.plot(np.abs(arr), dlvo.LewisABx[0], label='{} M'.format(j), lw=1.5)

plt.ylabel('Force (N)', fontsize=14)
plt.xlabel('X (m)', fontsize=14)
plt.legend(loc=0, fontsize=12)
plt.show()

for j in i:
    d['I'] = j
    dlvo = cm.DLVO(xarr, xarr, xvArr=varr, yvArr=varr, **d)
    plt.plot(np.abs(arr), dlvo.LVDWx[0], label='{} M'.format(j), lw=1.5)

plt.ylabel('Force (N)', fontsize=14)
plt.xlabel('X (m)', fontsize=14)
plt.legend(loc=0, fontsize=12)
plt.show()

arr = np.linspace(1e-9, 1e-7, 1000)
xarr = np.array([arr])
varr = np.array([np.ones(1000)])

for j in i:
    d['I'] = j
    dlvo = cm.DLVO(xarr, xarr, xvArr=varr, yvArr=varr, **d)
    plt.plot(np.abs(arr), (dlvo.LVDWx[0] + dlvo.EDLx[0] + dlvo.LewisABx[0]),
             label='{} M'.format(j), lw=1.5)

plt.plot([1e-9, 1e-7], [0, 0], ls='--', color='k', lw=1.5)
plt.ylabel('Force (N)', fontsize=14)
plt.xlabel('X (m)', fontsize=14)
plt.legend(loc=4, fontsize=12)
plt.show()

d['lbres'] = 1e-8

for j in i:
    d['I'] = j
    dlvo = cm.ColloidColloid(arr, **d)
    t = dlvo.x
    u = dlvo.y
    cdlvo = -1 * (dlvo.x + dlvo.y)
    plt.imshow(cdlvo, interpolation='None', cmap='jet', vmin=1e-6, vmax=1e-10, norm=matplotlib.colors.LogNorm())
    plt.plot([250], [250], 'ko')
    plt.ylim([0, 500])
    plt.xlim([0, 500])
    plt.colorbar().set_label("Force (N)", rotation=270, labelpad=20)
    plt.xlabel('nm')
    plt.ylabel('nm')
    plt.show()
"""
