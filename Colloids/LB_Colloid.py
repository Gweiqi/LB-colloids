import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import Colloid_Setup as cs
import Colloid_Math as cm
import Colloid_IO as IO
import random
from copy import copy


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
    positions = []

    def __init__(self, xlen, ylen, resolution):
        self.xposition = [random.uniform(0.1, 0.9)*xlen*resolution]
        self.yposition = [-resolution]
        self.resolution = resolution
        self.storey = [self.yposition[0]]
        self.storex = [self.xposition[0]]
        self.idx_rx = [int(self.xposition[-1]//resolution)]
        self.idx_ry = [-int(self.yposition[-1]//resolution)]
        self.__cell_time = [0.]
        self.time = [0]
        self.colloid_start_time = copy(TrackTime.model_time)
        self.colloid_end_time = np.nan
        self.flag = [1]
        self.ylen = ylen
        self.xlen = xlen
        Colloid.positions.append(tuple([self.idx_rx[-1], self.idx_ry[-1]]))

    def reset_master_positions(self):
        Colloid.positions = []

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

    def _append_idx_rx(self, item):
        self.idx_rx.append(item)

    def _append_idx_ry(self, item):
        self.idx_ry.append(item)

    def _append_master_positions(self, idx_rx, idx_ry):
        Colloid.positions.append(tuple([idx_rx, idx_ry]))

    def __update_cell_time(self, item, new_cell=False):
        if new_cell:
            self.__cell_time = [item]
        else:
            cell_time = self.__cell_time[-1] + item
            self.__cell_time.append(cell_time)

    def update_position(self, xvelocity, yvelocity, ts):
        """
        grid index method to update continuous colloid system using the discete grid forces.
        idxry must be inverted because we assume (0,0) at top left corner and down is negitive.
        """

        irx = self.xposition[-1]
        iry = self.yposition[-1]
        flag = self.flag[-1]
        # find grid indexes and look up velocity (negative y accounts for grid indexes bc of vector direction).

        #if flag == 2:
            # this is a NaN condition
        #    self.update_special(irx, iry, 2)

        if flag == 3:
            # this is the breakthrough condition
            irx = float('NaN')
            iry = float('NaN')
            self.update_special(irx, iry, 3)

        elif flag == 4:
            irx = float('NaN')
            iry = float('NaN')
            self.update_special(irx, iry, 4)

        else:
            # normal streaming conditions
            
            # if NaN flip flag to 2, we can then debug!
            try:
                idxrx = int(irx//self.resolution)
                idxry = int(iry//-self.resolution)
           
            except ValueError:
                if not np.isnan(self.xposition[-1]) and not np.isnan(self.yposition[-1]):
                    self._append_xposition(self.xposition[-1])
                    self._append_yposition(self.yposition[-1])
                    self.__update_cell_time(ts, new_cell=True)
                else:
                    self._append_xposition(random.uniform(0.1, 0.9) * self.xlen * self.resolution)
                    self._append_yposition(-self.resolution)
                    self._append_flag(2)
                return

        # if colloid breaks through domain change flag to 3  
            try:
                idxry = int(iry//-self.resolution)
                idxrx = int(irx//self.resolution)
                xv = xvelocity[idxry][idxrx]
                yv = yvelocity[idxry][idxrx]
                
            except IndexError:
                if idxry >= self.ylen:
                    self._append_flag(3)
                    self.colloid_end_time = copy(TrackTime.model_time)
                    self._append_xposition(float("NaN"))
                    self._append_yposition(float("NaN"))
                    return
                else:
                    self._append_xposition(random.uniform(0.1, 0.9) * self.xlen * self.resolution)
                    self._append_yposition(-self.resolution)
                    self._append_flag(2)
                    # self._append_flag(2)
                    # self._append_xposition(irx)
                    # self._append_yposition(iry)
                    return

            # track time each colloid spends in a cell to account for acceleration effects
            if idxrx == self.idx_rx[-1] and idxry == self.idx_ry[-1]:
                self.__update_cell_time(ts)
            else:
                self.__update_cell_time(ts, new_cell=True)

            xv = xvelocity[idxry][idxrx]
            yv = yvelocity[idxry][idxrx]
        
            # velocity is L/T, but since we are dealing with acceleration fields 0.5*(a)t^2 is used
            deltarx = xv * ts  # self.__cell_time[-1]
            deltary = yv * ts  # self.__cell_time[-1]

            rx = irx + deltarx
            ry = iry + deltary
            
            # Use a check to make sure colloid does not leave domain on first iteration
            # move colloid to a new location to if it leaves domain in the vertical y-direction
            if ry >= 0.:
                rx = random.uniform(0.1, 0.9)*self.xlen*self.resolution
                ry = -self.resolution
                self.__update_cell_time(ts, new_cell=True)

            self._append_master_positions(idxrx, idxry)
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
    model_time = 1

    def __init__(self, ts):
        self.ts = ts
        self.timer = [1]
        self.time = self.timer[-1]
        self.totim = [self.time*self.ts]

    def update_time(self):
        self.time += 1
        TrackTime.model_time += 1
        self.timer.append(self.time)
        self.totim.append(self.time*self.ts)

    def strip_time(self):
        self.timer = [self.timer[-1]]
        self.totim = [self.totim[-1]]

    def print_time(self):
        print(self.timer[-1], "%.3f" % self.totim[-1])


def fmt(x, pos):
    # function formatter for matplotlib graphs
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def run_save_model(x, iters, vx, vy, ts, timer, print_time, store_time,
                   colloidcolloid, pathline=None, timeseries=None, endpoint=None):
    """
    definition to allow the use of multiple ionic strengths ie. attachment then flush, etc....
    """
    vx0 = copy(vx)
    vy0 = copy(vy)
    # col_col
    while timer.time <= iters:
        # update colloid position and time
        if x:
            x[0].reset_master_positions()
        for col in x:
            col.update_position(vx, vy, ts)

        timer.update_time()

        # check for a printing prompt
        if timer.time % print_time == 0.:
            timer.print_time()
            
        # check store_times and strip younger time steps from memory
        if timer.time % store_time == 0.:
            if pathline is not None:
                pathline.write_output(timer, x)
                timer.strip_time()
                for col in x:
                    col.store_position(timer)  # use this for plotting functionality
                    col.strip_positions()
            
            elif timeseries is not None:
                for col in x:
                    col.store_position(timer)  # use this for plotting functionality
                    col.strip_positions()
                timeseries.write_output(timer, x, pathline=False)

            else:
                col.store_position(timer)
                col.strip_positions()
                
        # check if user wants an endpoint file
        if timer.time == iters:
            if endpoint is not None:
                for col in x:
                    col.store_position(timer)
                    col.strip_positions()
                endpoint.write_output(timer, x, pathline=False)

    del colloidcolloid

def run(config):
    """
    Model definition to setup and run the LB_Colloids from config file
    or from OO.

    config: (class colloid_IO.Config) object or list of colloid_IO.Config objects
    :return:
    """

    if isinstance(config, list):
        multiple_config = config[1:]
        config = config
    else:
        multiple_config = []

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
    preferential_flow = False

    if 'multiple_config' in ModelDict:
        if ModelDict['multiple_config']:
            assert 'nconfig' in ModelDict
        else:
            ModelDict['multiple_config'] = False
    else:
        ModelDict['multiple_config'] = False

    # setup output and boolean flags.    
    if 'print_time' in OutputDict:
        print_time = OutputDict['print_time']
    else:
        print_time = iters

    pathline = None
    timeseries = None
    endpoint = None
    
    if 'pathline' in OutputDict:
        pathline = IO.Output(OutputDict['pathline'], **OutputDict)
        if 'store_time' not in OutputDict:
            OutputDict['store_time'] = 100
            
    if 'timeseries' in OutputDict:
        timeseries = IO.Output(OutputDict['timeseries'], **OutputDict)
        assert 'store_time' in OutputDict, 'please provide STORE_TIME as interval time'
        
    if 'endpoint' in OutputDict:
        endpoint = IO.Output(OutputDict['endpoint'], **OutputDict)
        if 'store_time' not in OutputDict:
            OutputDict['store_time'] = 100

    # implemented for memory management purposes
    if 'store_time' not in OutputDict:
        store_time = 100
    else:
        store_time = OutputDict['store_time']

    print(ncols)
    # get data from LB Model
    LB = cs.HDF5_reader(modelname)

    LBy = cs.LBVArray(LB.yu, LB.imarray)
    LBx = cs.LBVArray(LB.xu, LB.imarray)
    velocity_factor = LB.velocity_factor

    # interpolate over grid array and interpolate veloctity profiles
    # Col_vp = interp_v(LBvp, gridsplit)
    LBy = cs.InterpV(LBy, gridsplit)
    LBx = cs.InterpV(LBx, gridsplit)
    Col_img = cs.InterpV(LB.imarray, gridsplit, img=True)
    
    # Use grids function to measure distance from pore space and correct
    # boundaries for interpolation effects. 
    Grids = cs.Gridarray(Col_img, gridres, gridsplit)
    xArr = Grids.gridx
    yArr = Grids.gridy
    
    xvArr = Grids.vector_x
    yvArr = Grids.vector_y

    xArr[Col_img == 1.] = np.nan
    yArr[Col_img == 1.] = np.nan
    
    xvArr[Col_img == 1.] = np.nan
    yvArr[Col_img == 1.] = np.nan

    # Begin calling colloid mathematics
    cfactor = cm.Gap(xArr, yArr)

    # for calculations we need ts = 1., adjust later
   
    velocity = cm.Velocity(LBx, LBy, velocity_factor, **PhysicalDict)
    LBx = velocity.xvelocity
    LBy = velocity.yvelocity
    
    # Initial setup block to estimate colloid velocity for drag_force calc. 
    drag_forces = cm.Drag(LBx, LBy, cfactor.f1, cfactor.f2, cfactor.f3,
                          cfactor.f4, **PhysicalDict)

    brownian = cm.Brownian(cfactor.f1, cfactor.f4, **PhysicalDict)

    dlvo = cm.DLVO(xArr, yArr, xvArr=xvArr, yvArr=yvArr, **ChemicalDict)

    gravity = cm.Gravity(**PhysicalDict)
    bouyancy = cm.Bouyancy(**PhysicalDict)

    physicalx = brownian.brownian_x + drag_forces.drag_x
    physicaly = brownian.brownian_y + drag_forces.drag_y + gravity.gravity + bouyancy.bouyancy

    dlvox = dlvo.EDLx + dlvo.LVDWx + dlvo.LewisABx
    dlvoy = dlvo.EDLy + dlvo.LVDWy + dlvo.LewisABy
    colloidcolloid = cm.ColloidColloid(Col_img, **ChemicalDict)

    if preferential_flow is True:
        fx = dlvox
        fy = dlvoy
    else:  
        fx = dlvox + physicalx
        fy = dlvoy + physicaly

    vx = cm.ForceToVelocity(fx, **PhysicalDict)
    vy = cm.ForceToVelocity(fy, **PhysicalDict)

    # get LB velocity to add to the physical forces calculated.
    LBx = velocity.xvelocity
    LBy = velocity.yvelocity
  
    vx = vx.velocity + LBx
    vy = vy.velocity + LBy

    Col_img = cs.InterpV(LB.imarray, gridsplit, img=True)
    vx[Col_img == 1.] = np.nan
    vy[Col_img == 1.] = np.nan
    
    ylen = len(Col_img)
    xlen = len(Col_img[0])
    x = [Colloid(xlen, ylen, gridres) for i in range(ncols)]

    # start model timer
    timer = TrackTime(ts)

    run_save_model(x, iters, vx, vy, ts, timer, print_time,
                   store_time, colloidcolloid,
                   pathline, timeseries, endpoint)

    if ModelDict['multiple_config'] is True:
        for confignumber in range(0, ModelDict['nconfig']-1):
            config = IO.Config(multiple_config[confignumber])
            ModelDict = config.model_parameters()
            PhysicalDict = config.physical_parameters()
            ChemicalDict = config.chemical_parameters()

            iters = ModelDict['iters']

            # recalculate all physical chemical forces and continue running model
            drag_forces = cm.Drag(LBx, LBy, cfactor.f1, cfactor.f2, cfactor.f3,
                                  cfactor.f4, **PhysicalDict)

            brownian = cm.Brownian(cfactor.f1, cfactor.f4, **PhysicalDict)

            dlvo = cm.DLVO(xArr, yArr, xvArr=xvArr, yvArr=yvArr, **ChemicalDict)

            gravity = cm.Gravity(**PhysicalDict)
            bouyancy = cm.Bouyancy(**PhysicalDict)

            physicalx = brownian.brownian_x + drag_forces.drag_x
            physicaly = brownian.brownian_y + drag_forces.drag_y + gravity.gravity + bouyancy.bouyancy

            dlvox = dlvo.EDLx + dlvo.LVDWx + dlvo.LewisABx 
            dlvoy = dlvo.EDLy + dlvo.LVDWy + dlvo.LewisABy

            if preferential_flow is True:
                fx = dlvox
                fy = dlvoy
            else:
                fx = dlvox + physicalx
                fy = dlvoy + physicaly
    
            vx = cm.ForceToVelocity(fx, **PhysicalDict)
            vy = cm.ForceToVelocity(fy, **PhysicalDict)
            
            LBx = velocity.xvelocity
            LBy = velocity.yvelocity
 
            vx = vx.velocity + LBx
            vy = vy.velocity + LBy

            run_save_model(x, iters, vx, vy, ts, timer, print_time,
                           store_time, colloidcolloid,
                           pathline, timeseries, endpoint)

    if OutputDict['plot']:
        # set up option for vy vs. LBy plotting
        # re-interpolate Coloid image object to get binary array
        Col_img = cs.InterpV(LB.imarray, gridsplit, img=True)
        LBy[Col_img == 1] = np.nan
        LBy = np.ma.masked_invalid(LBy)

        # setup meshgrid for precise plotting
        xx = np.arange(xlen + 1)
        xx = np.tile(xx, (ylen + 1, 1))

        yy = np.array([np.arange(ylen + 1)])
        yy = np.tile(yy.T, (1, xlen + 1))
        
        plt.pcolormesh(xx, yy, LBy, cmap='viridis_r') 
        
        for col in x:
            plt.plot(np.array(col.storex)/gridres, np.array(col.storey)/-gridres, 'o',
                     ms=8)
        plt.xlim([0, xlen])
        plt.ylim([ylen, 0])
    
        cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
        cbar.set_label('m/s', rotation=270)
    
        plt.show()

    IO.HDF5WriteArray(velocity.xvelocity,
                      velocity.yvelocity,
                      ModelDict,
                      dlvo.all_chemical_params,
                      drag_forces.all_physical_params)

if __name__ == '__main__':
    # todo: Need to fix this issue to check for multiple config
    # todo: and then send it through as a list of dictionaries
    config_file = IO.Config('Synthetic.config')
    run(config_file)

else:
    pass
