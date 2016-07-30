import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import Colloid_Setup as cs
import Colloid_Math as cm
import sys
import optparse

LB = cs.HDF5_reader('Synthetic255.hdf5')

#### gridsplit needs to be one of the many config options!
gridsplit = 10
gridres = 1e-6/gridsplit

LBy = cs.LBVArray(LB.yu, LB.imarray)
LBx = cs.LBVArray(LB.xu, LB.imarray)

#### interpolate over grid array and interpolate veloctity profiles ####
#Col_vp = interp_v(LBvp, gridsplit)
LBy = cs.InterpV(LBy, gridsplit)
LBx = cs.InterpV(LBx, gridsplit)
Col_img = cs.InterpV(LB.imarray, gridsplit)

### Use grids function to measure distance from pore space and correct
### boundaries for interpolation effects. 
Grids = cs.Gridarray(Col_img, gridres, gridsplit)
xArr = Grids.gridx
yArr = Grids.gridy

xvArr = Grids.vector_x
yvArr = Grids.vector_y

xArr[LBx == 0.] = np.nan
yArr[LBy == 0.] = np.nan
xvArr[LBx == 0.] = np.nan
yvArr[LBy == 0.] = np.nan

# Begin calling colloid mathematics #

cfactor = cm.Gap(xArr, yArr)

velocity = cm.Velocity(LBx, LBy, 1., gridres)
LBx = velocity.xvelocity
LBy = velocity.yvelocity

drag_forces = cm.Drag(LBx, LBy, LBx, LBy, cfactor.f1, cfactor.f2,
                      cfactor.f3, cfactor.f4, xvArr, yvArr) 

brownian = cm.Brownian(xArr, yArr, cfactor.f1, cfactor.f4)

dlvo = cm.DLVO(xArr, yArr, xvArr=xvArr, yvArr=yvArr)

gravity = cm.Gravity()
bouyancy = cm.Bouyancy()

physicalx = brownian.brownian_x + drag_forces.drag_x
physicaly = brownian.brownian_y + drag_forces.drag_y + gravity.gravity + bouyancy.bouyancy

dlvox = dlvo.EDLx + dlvo.LVDWx + dlvo.LewisABx
dlvoy = dlvo.EDLy + dlvo.LVDWy + dlvo.LewisABy

fx = dlvox + physicalx
fy = dlvoy + physicaly

vx = cm.ForceToVelocity(fx)
vy = cm.ForceToVelocity(fy)

vx = vx.velocity + LBx
vy = vy.velocity + LBy

plt.imshow(vx, interpolation='nearest', vmin=-1e-13, vmax=1e-13)
plt.colorbar()
plt.show()

plt.imshow(vy, interpolation='nearest', vmin=-1e-13, vmax=1e-13)
plt.colorbar()
plt.show()
