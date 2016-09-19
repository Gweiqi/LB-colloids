import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import Colloid_Setup as cs
import Colloid_Math as cm
import Colloid_IO as IO
import sys
import optparse
import random

class Colloid:
    """
    Wrapper class to initiate and track colloid position through the LB Model

    Inputs:
    -------
    xlen: (int) grid length after interpolation in the x-direction
    resolution: (float) grid resolution after interpolation

    Methods:
    --------
    _append_xposition: method object that appends the list tracking a colloids x position
    _append_yposition: method object that appends the list tracking a colloids y position
    update_position: method that updates postion of colloids by calling the append methods
    strip_positions: method that strips all but last item from colloid position lists
    
    Returns:
    --------
    xposition: (list, float) a list of x-position values normalized to the grid resolution
    yposition: (list, float) a list of y-position values normalized to grid resolution (top left is 0,0)
    """
    def __init__(self, xlen, resolution):
        self.xposition = [random.uniform(0.01,0.99)*xlen*resolution]
        self.yposition = [0.]
        self.resolution = resolution
        self.storey = [self.yposition[0]]
        self.storex = [self.xposition[0]]
        self.time = [0]
        self.flag = [1]

    def _append_xposition(self, item):
        self.xposition.append(item)

    def _append_yposition(self, item):
        self.yposition.append(item)

    def _append_storex(self, item):
        self.storex.append(item)

    def _append_storey(self, item):
        self.storey.append(item)
        
    def _append_time(self, item):
        self.time.append(item)

    def _append_flag(self, item):
        self.flag.append(item)

    def update_position(self, xvelocity, yvelocity, ts):
        """
        grid index method to update continuous colloid system using the discete grid forces.
        idxry must be inverted because we assume (0,0) at top left corner and down is negitive.
        """
        irx = self.xposition[-1]
        iry = self.yposition[-1]
        flag = self.flag[-1]
        # find grid indexes and look up velocity (negative y accounts for grid indexes bc of vector direction).

        if flag == 2:
            # this is a NaN condition
            self.update_special(irx, iry, 2)

        elif flag == 3:
            # this is the breakthrough condition
            irx = float('NaN')
            iry = float('NaN')
            self.update_special(irx, iry, 3)

        else:
            # normal streaming conditions
            
            # if NaN flip flag to 2, we can then debug!
            try:
                idxrx = int(irx//self.resolution)
                idxry = int(iry//-self.resolution)
           
            except ValueError:
                self._append_flag(2)
                self._append_xposition(irx)
                self._append_yposition(iry)
                return

        # if colloid breaks through domain change flag to 3  
            try:
                idxry = int(iry//-self.resolution)

            except IndexError:
                self._append_flag(3)
                self._append_xposition(float("NaN"))
                self._append_yposition(float("NaN"))

            xv = xvelocity[idxry][idxrx]
            yv = yvelocity[idxry][idxrx]
        
            # velocity is L/T therefore multiply by T
            deltarx = xv*ts 
            deltary = yv*ts
            rx = irx + deltarx
            ry = iry + deltary
            # Use a check to make sure colloid does not leave domain on first iteration
            # move colloid to a new location to if it leaves domain in the vertical y-direction
            if ry >= 0:
                rx = random.uniform(0.05,0.95)*xlen*self.resolution
                ry = 0.
            else:
                pass

            self._append_xposition(rx)
            self._append_yposition(ry)
            self._append_flag(1)

    def update_special(self, irx, iry, flag):
        self._append_xposition(irx)
        self._append_yposition(iry)
        self._append_flag(flag)

    def strip_positions(self):
        self.xposition = [self.xposition[-1]]
        self.yposition = [self.yposition[-1]]

    def store_position(self, timer):
        self._append_time(timer)
        self._append_storex(self.xposition[-1])
        self._append_storey(self.yposition[-1])
        
    def print_positions(self):
        print(self.xposition[-1], self.yposition[-1])


class TrackTime:
    """
    TrackTime class is the model timer, enables stripping stored time steps which
    is useful to free memoryafter storing the data externally. Is necessary for
    output class functionality!
    """
    def __init__(self, ts):
        self.ts = ts
        self.timer = [0]
        self.time = self.timer[-1]
        self.totim = [self.time*self.ts]

    def update_time(self):
        self.time = self.time + 1
        self.timer.append(self.time)
        self.totim.append(self.time*self.ts)

    def strip_time(self):
        self.timer = [self.timer[-1]]
        self.totim = [self.totim[-1]]

    def print_time(self):
        print(self.timer[-1], "%.3f" % self.totim[-1])


if __name__ == '__main__':
    config = IO.Config('Synthetic.config')
    ModelDict = config.model_parameters()
    PhysicalDict = config.physical_parameters()
    ChemicalDict = config.chemical_parameters()
    OutputDict = config.output_control()

    # set model variable block
    modelname = ModelDict['lbmodel']
    gridsplit = ModelDict['gridref']
    lbres = ModelDict['lbres']
    gridres = lbres/gridsplit
    ts = ModelDict['ts']
    iters = ModelDict['iters']
    ncols = ModelDict['ncols']

    # setup output and boolean flags.    
    if 'print_time' in OutputDict:
        print_time = OutputDict['print_time']
    else:
        print_time = iters

    isitpathline = False
    isittimeseries = False
    isitendpoint = False
    
    if 'pathline'  in OutputDict:
        pathline = IO.Output(OutputDict['pathline'], **OutputDict)
        isitpathline = True
        if 'store_time' not in OutputDict:
            OutputDict['store_time'] = 100
            
    if 'timeseries' in OutputDict:
        timeseries = IO.Output(OutputDict['timeseries'], **OutputDict)
        isittimeseries = True
        assert 'store_time' in OutputDict, 'please provide STORE_TIME as interval time'
        
    if 'endpoint' in OutputDict:
        endpoint = IO.Output(OutputDict['endpoint'], **OutputDict)
        isitendpoint = True
        if 'store_time' not in OutputDict:
            OutputDict['store_time'] = 100

    # implemented for memory management purposes
    if 'store_time' not in OutputDict:
        store_time = 100
    else:
        store_time = OutputDict['store_time']

    print(store_time)
    # get data from LB Model
    LB = cs.HDF5_reader(modelname)

    LBy = cs.LBVArray(LB.yu, LB.imarray)
    LBx = cs.LBVArray(LB.xu, LB.imarray)

    #### interpolate over grid array and interpolate veloctity profiles ####
    #Col_vp = interp_v(LBvp, gridsplit)
    LBy = cs.InterpV(LBy, gridsplit)
    LBx = cs.InterpV(LBx, gridsplit)
    Col_img = cs.InterpV(LB.imarray, gridsplit)

    # Use grids function to measure distance from pore space and correct
    # boundaries for interpolation effects. 
    Grids = cs.Gridarray(Col_img, gridres, gridsplit)
    xArr = Grids.gridx
    yArr = Grids.gridy

    xvArr = Grids.vector_x
    yvArr = Grids.vector_y

    xArr[LBx == 0.] = np.nan
    yArr[LBy == 0.] = np.nan
    xvArr[LBx == 0.] = np.nan
    yvArr[LBy == 0.] = np.nan

    # Begin calling colloid mathematics 

    cfactor = cm.Gap(xArr, yArr)

    # for calculations we need ts = 1., adjust later
    velocity = cm.Velocity(LBx, LBy, gridres)
    LBx = velocity.xvelocity
    LBy = velocity.yvelocity

    # Initial setup block to estimate colloid velocity for drag_force calc. 
    drag_forces = cm.Drag(LBx, LBy, cfactor.f1, cfactor.f2, cfactor.f3,
                          cfactor.f4, xvArr, yvArr, **PhysicalDict)

    brownian = cm.Brownian(xArr, yArr, cfactor.f1, cfactor.f4, **PhysicalDict)

    dlvo = cm.DLVO(xArr, yArr, xvArr=xvArr, yvArr=yvArr, **ChemicalDict)

    gravity = cm.Gravity(**PhysicalDict)
    bouyancy = cm.Bouyancy(**PhysicalDict)

    physicalx = brownian.brownian_x + drag_forces.drag_x
    physicaly = brownian.brownian_y + drag_forces.drag_y + gravity.gravity + bouyancy.bouyancy

    dlvox = dlvo.EDLx + dlvo.LVDWx + dlvo.LewisABx
    dlvoy = dlvo.EDLy + dlvo.LVDWy + dlvo.LewisABy

    fx = dlvox + physicalx
    fy = dlvoy + physicaly
    
    vx = cm.ForceToVelocity(fx, **PhysicalDict)
    vy = cm.ForceToVelocity(fy, **PhysicalDict)

    # adjust LBvelocity by timestep, move this to cm.velocity
    velocity = cm.Velocity(LBx, LBy, gridres, **PhysicalDict)
    LBx = velocity.xvelocity
    LBy = velocity.yvelocity
    
    vx = vx.velocity + LBx
    vy = vy.velocity + LBy
    
    xlen = len(Col_img[0])
    x = [Colloid(xlen, gridres) for i in range(ncols)]
    
    timer = TrackTime(ts)
    while timer.time <= iters:
        # update colloid position and time
        for col in x:
            col.update_position(vx, vy, ts)
        timer.update_time()

        # check for a printing prompt
        if timer.time%print_time == 0.:
            timer.print_time()
            
        # check store_times and strip younger time steps from memory
        if timer.time%store_time == 0.:
            if isitpathline is True:
                pathline.write_output(timer, x)
                timer.strip_time()
                for col in x:
                    col.store_position(timer) # use this for plotting functionality
                    col.strip_positions()
            
            elif isittimeseries is True:
                for col in x:
                    col.store_position(timer) # use this for plotting functionality
                    col.strip_positions()
                timeseries.write_output(timer, x, pathline=False)

            else:
                col.store_positions(timer)
                col.strip_positions()
                
        # check if user wants an endpoint file
        if timer.time == iters:
            if isitendpoint is True:
                for col in x:
                    col.store_position(timer)
                    col.strip_positions()
                endpoint.write_output(timer, x, pathline=False)

    if OutputDict['plot'] is True:  
        plt.imshow(vy, interpolation='nearest', vmin=-1e-13, vmax=1e-13)
        for col in x:
            plt.plot(np.array(col.storex)/gridres, np.array(col.storey)/-gridres, 'o',
                     ms=8)
        plt.colorbar()
        plt.show()
