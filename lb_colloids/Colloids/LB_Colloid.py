"""
The LB_Colloid module contains base classes to simulate colloid transport for
colloid simulation. This module acts as the control center. The class Colloid
is the base representation of a colloid and contains
streaming and updating rules. The class TrackTime is the simulation timer, which tracks
both number of time steps and the time step length. Also of importance is the run()
method. This method initiates a colloid simulation from a IO.Config object.

A user can initiate a model run with a Colloid_IO.Config() object.
Please see the Input Output section for details on building the Colloid_IO.Config() object

>>> from lb_colloids import ColloidModel
>>>
>>> config = IO.Config()  # We assume that the Colloid_IO.Config() object is already built. See the Colloid_IO section for details
>>> ColloidModel.run(config)
"""
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
    Primary colloid class to instantiate and track colloids through a colloid simulation

    Parameters:
    ----------
    :param int xlen: grid length after interpolation in the x-direction
    :param float resolution: grid resolution after interpolation

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
        """
        Resets the master position storage mechanism for all colloids. Master position
        storage is used to later generate colloid-colloid DLVO fields.
        """
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
        Colloid.positions.append([idx_rx, idx_ry])

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

        Parameters:
        ----------
        :param np.ndarray xvelocity: colloid simulation velocity array in the x domain
        :param np.ndarray yvelocity: colloid simulation velocity array in the y domain
        :param ts float: physical time step in seconds
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
            # irx = float('NaN')
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
                    self._append_xposition(irx)
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

            # more appropriate use 0.5 * 'velocity' * ts, but separate terms?
            deltarx = xv * ts
            deltary = yv * ts

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
        """
        Special updater class for colloids that have exited the model
        domain or experienced an internal error

        Parameters:
        ----------
        :param int irx: grid index in the x domain
        :param int iry: grid index in the y domain
        :param int flag: number indicator of special condition. 3 = normal breakthrough condition
        """
        self._append_xposition(irx)
        self._append_yposition(iry)
        self._append_flag(flag)

    def strip_positions(self):
        """
        Memory saving function to strip unused, save colloid position information
        """
        self.xposition = [self.xposition[-1]]
        self.yposition = [self.yposition[-1]]

    def store_position(self, timer):
        """
        Method to store colloid position and update the the time of storage

        :param float timer: current model time
        """
        self._append_time(timer)
        self._append_storex(self.xposition[-1])
        self._append_storey(self.yposition[-1])
        
    def print_positions(self):
        print(self.xposition[-1], self.yposition[-1])


class TrackTime:
    """
    TrackTime class is the model timer. This class enables stripping stored time steps which
    is useful to free memory after writing the data to an external file. Is necessary for
    output class functionality!

    :param int ts: model time step set by user (physical time)
    """
    model_time = 1

    def __init__(self, ts):
        TrackTime.model_time = 1
        self.ts = ts
        self.timer = [1]
        self.time = self.timer[-1]
        self.totim = [self.time*self.ts]

    def update_time(self):
        """
        Method to update the current model time and store it
        """
        self.time += 1
        TrackTime.model_time += 1
        self.timer.append(self.time)
        self.totim.append(self.time*self.ts)

    def strip_time(self):
        """
        Memory saving method that removes unused information from timer arrays
        """
        self.timer = [self.timer[-1]]
        self.totim = [self.totim[-1]]

    def print_time(self):
        """
        Method prints time in a standard format to the terminal
        """
        print(self.timer[-1], "%.3f" % self.totim[-1])


def fmt(x, pos):
    # function formatter for matplotlib graphs
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)


def _run_save_model(x, iters, vx, vy, ts, xlen, ylen, gridres,
                    ncols, timer, print_time, store_time,
                    colloidcolloid, ModelDict, pathline=None,
                    timeseries=None, endpoint=None):
    """
    definition to allow the use of multiple ionic strengths ie. attachment then flush, etc....
    """
    continuous = 0
    if 'continuous' in ModelDict:
        continuous = ModelDict['continuous']

    colloidcolloid.update(x)
    conversion = cm.ForceToVelocity(1, **ModelDict).velocity

    while timer.time <= iters:
        # update colloid position and time
        if continuous:
            if timer.time % continuous == 0 and timer.time != 0:
                x += [Colloid(xlen, ylen, gridres) for i in range(ncols)]

        colloidcolloid.update(x)
        cc_vx = colloidcolloid.x_array * conversion  # /1e-6
        cc_vy = colloidcolloid.y_array * conversion  # /1e-6
        Colloid.positions = []
        vx0 = vx + cc_vx
        vy0 = vy + cc_vy

        for col in x:
            col.update_position(vx0, vy0, ts)

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
                for col in x:
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
    Model definition to setup and run the LB_Colloids from config file.
    This is the also the preferred user interaction method to run simulations
    when working with python.

    Parameters:
    ----------
    :param colloid_IO.config config: object or list of colloid_IO.Config objects

    """

    if isinstance(config, list):
        if len(config) > 1:
            multiple_config = config[1:]
            config = config[0]
        else:
            multiple_config = []
            config = config[0]
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
    # call the HDF5 array early to get domian size for output
    LB = cs.Hdf5Reader(modelname)

    OutputDict['xlen'] = LB.imarray.shape[1] * gridsplit
    OutputDict['ylen'] = LB.imarray.shape[0] * gridsplit

    OutputDict['mean_ux'] = LB.mean_xu
    OutputDict['mean_uy'] = LB.mean_yu
    OutputDict['velocity_factor'] = LB.velocity_factor

    if 'continuous' in ModelDict:
        OutputDict['continuous'] = ModelDict['continuous']
    else:
        OutputDict['continuous'] = 0

    # setup output and boolean flags.    
    if 'print_time' in OutputDict:
        print_time = OutputDict['print_time']
    else:
        print_time = iters

    if 'showfig' in OutputDict:
        pass
    else:
        OutputDict['showfig'] = False

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
    LBy = cs.LBVArray(LB.yu, LB.imarray)
    LBx = cs.LBVArray(LB.xu, LB.imarray)
    velocity_factor = LB.velocity_factor

    # interpolate over grid array and interpolate veloctity profiles
    LBy = cs.InterpV(LBy, gridsplit)
    LBx = cs.InterpV(LBx, gridsplit)
    Col_img = cs.InterpV(LB.imarray, gridsplit, img=True)
    
    # Use grids function to measure distance from pore space and correct
    # boundaries for interpolation effects. 
    Grids = cs.GridArray(Col_img, gridres, gridsplit)
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

    dlvox = dlvo.EDLx + dlvo.attractive_x  # + dlvo.LVDWx + dlvo.LewisABx
    dlvoy = dlvo.EDLy + dlvo.attractive_y  # dlvo.LVDWy + dlvo.LewisABy
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

    # start model timer
    timer = TrackTime(ts)

    # generate initial pulse of colloids.
    x = [Colloid(xlen, ylen, gridres) for i in range(ncols)]

    _run_save_model(x, iters, vx, vy, ts, xlen, ylen, gridres,
                    ncols, timer, print_time,
                    store_time, colloidcolloid, ModelDict,
                    pathline, timeseries, endpoint)

    if multiple_config:
        for confignumber in range(len(multiple_config)):
            config = multiple_config[confignumber]
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

            dlvox = dlvo.EDLx + dlvo.attractive_x  # dlvo.LVDWx + dlvo.LewisABx
            dlvoy = dlvo.EDLy + dlvo.attractive_y  # dlvo.LVDWy + dlvo.LewisABy

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

            _run_save_model(x, iters, vx, vy, ts, xlen, ylen, gridres,
                            ncols, timer, print_time,
                            store_time, colloidcolloid, ModelDict,
                            pathline, timeseries, endpoint)

    if OutputDict['plot']:
        # set up option for vy vs. LBy plotting
        # re-interpolate Coloid image object to get binary array
        Col_img = cs.InterpV(LB.imarray, gridsplit, img=True)
        LBy[Col_img == 1] = np.nan
        LBy = np.ma.masked_invalid(LBy)
        LBx[Col_img == 1] = np.nan
        LBx = np.ma.masked_invalid(LBx)

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

        if OutputDict['showfig']:
            plt.show()

        else:
            try:
                plt.savefig(OutputDict['endpoint'].strip('endpoint') + "png")
                plt.close()
            except:
                plt.show()

    else:
        # mask the velocity objects for later output plotting
        LBy[Col_img == 1] = np.nan
        LBy = np.ma.masked_invalid(LBy)
        LBx[Col_img == 1] = np.nan
        LBx = np.ma.masked_invalid(LBx)

    IO.HDF5WriteArray(velocity.xvelocity,
                      velocity.yvelocity,
                      colloidcolloid,
                      ModelDict,
                      dlvo.all_chemical_params,
                      drag_forces.all_physical_params)


if __name__ == '__main__':
    pass

