import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py as H


# todo: return an axis object and dump the data into csv for further processing?

class Breakthrough(object):
    """
    Class to prepare and plot breakthrough curve data from endpoint
    files.
    Parameters:
        filename (str) <>.endpoint file

    Attributes:
        df: (pandas DataFrame) dataframe of endpoint data
        resolution: (float) model resolution
        timestep: (float) model timestep
        ncol: (float) number of colloids in simulation

    Properties:
        breakthrough_curve: (pandas DataFrame) data frame of
            raw breakthrough data.
    """
    def __init__(self, filename):
        if not filename.endswith('.endpoint'):
            raise FileTypeError('.endpoint file must be supplied')
        reader = ASCIIReader(filename)

        self.df = reader.df
        self.resolution = reader.resolution
        self.timestep = reader.timestep
        # todo: replace this call with something from the header later!
        self.ncol = float(self.df.shape[0])
        self.__breakthrough_curve = None
        self.__reader = reader

    @property
    def breakthrough_curve(self):
        """
        Dynamic calculation of breakthrough curve data

        Returns
            self.__breakthrough_curve
        """
        max_ts = self.df['nts'].max()
        if self.__breakthrough_curve is None:
            bt_colloids = self.df.loc[self.df['flag'] == 3]
            bt_colloids = bt_colloids.sort_values('end-ts')

            ncols = []
            nts = []
            ncol = 0
            for index, row in bt_colloids.iterrows():
                ncol += 1
                ncols.append(ncol)
                nts.append(row['end-ts'])

            ncols.append(ncol)
            nts.append(max_ts)

            df = pd.DataFrame({'nts': nts, 'ncol': ncols}).set_index('ncol')

            self.__breakthrough_curve = df

        return self.__breakthrough_curve

    def pore_volume_conversion(self):
        """
        Method to retrieve the pore volume calculation
        conversion for plotting colloids.

        Returns:
            pv_factor (float)
        """
        pv_factor = abs(self.__reader.uy)/(self.__reader.ylen * self.resolution)
        return pv_factor

    def plot(self, time=True, *args, **kwargs):
        """
        Convience method to plot data into a matplotlib
        chart.

        Parameters:
            time: (bool) if true x-axis is time, false is nts
            args: matplotlib args for 1d charts
            kwargs: matplotlib keyword arguments for 1d charts
        """
        if time:
            plt.plot(self.breakthrough_curve['nts'] * self.timestep,
                     self.breakthrough_curve.index.values / self.ncol,
                     *args, **kwargs)

        else:
            plt.plot(self.breakthrough_curve['nts'],
                     self.breakthrough_curve.index.values / self.ncol,
                     *args, **kwargs)
        plt.ylim([0, 1])

    def plot_pv(self, *args, **kwargs):
        """
        Method to plot pdf data with pore volumes (non-dimensional time)

        Parameters:
            args: matplotlib args for 1d plotting
            kwargs: matplotlib kwargs for 1d plotting
        """
        pv_factor = self.pore_volume_conversion()
        plt.plot(self.breakthrough_curve['nts'] * pv_factor * self.timestep,
                 self.breakthrough_curve.index.values / self.ncol,
                 *args, **kwargs)

        plt.ylim([0, 1])


class DistributionFunction(object):
    """
    Class to plot a pdf function of colloid breakthrough
    from endpoint files.

    Parameters:
        filename: (str)
        nbin: number of bins for pdf calculation

    Attributes:
        df: (pandas DataFrame) dataframe of endpoint data
        resolution: (float) model resolution
        timestep: (float) model timestep
        ncol: (float) number of colloids in simulation
        pdf: (np.recarray) colloid pdf

    Methods:
        reset_pdf: Allows user to adjust the bin size and
            recalculate the pdf
    """
    def __init__(self, filename, nbin=1000):
        if not filename.endswith('.endpoint'):
            raise FileTypeError('.endpoint file must be supplied')
        reader = ASCIIReader(filename)
        self.df = reader.df
        self.resolution = reader.resolution
        self.timestep = reader.timestep
        # todo: replace this call with something from the header later!
        self.ncol = float(self.df.shape[0])
        self.bin = nbin
        self.pdf = None
        self.reset_pdf(nbin)
        self.__reader = reader

    def reset_pdf(self, nbin):
        """
        Method to generate a probability distribution function
        based upon user supplied bin size.

        Parameters:
            nbin: (int) number of time steps to base bin on

        Returns:
            () probability distribution function of colloid
            breakthrough.
        """
        self.bin = nbin
        ts = []
        ncols = []
        lower_nts = 0
        max_ts = self.df['nts'].max()
        pdf_colloids = self.df.loc[self.df['flag'] == 3]
        pdf_colloids = pdf_colloids.sort_values('delta-ts')

        for upper_nts in range(0, int(max_ts) + 1, nbin):
            ncol = 0
            for index, row in pdf_colloids.iterrows():
                if lower_nts < row['delta-ts'] <= upper_nts:
                    ncol += 1

            ts.append(upper_nts)
            ncols.append(ncol)
            lower_nts = upper_nts

        arr = np.recarray((len(ts),), dtype=[('nts', np.float),
                                             ('ncol', np.float)])
        for idx, value in enumerate(ts):
            arr[idx] = tuple([value, ncols[idx]])

        self.pdf = arr

    def pore_volume_conversion(self):
        """
        Method to retrieve the pore volume calculation
        conversion for plotting colloids.

        Returns:
            pv_factor (float)
        """
        pv_factor = abs(self.__reader.uy)/(self.__reader.ylen * self.resolution)
        return pv_factor

    def plot(self, time=True, *args, **kwargs):
        """
        Convience method to plot data into a matplotlib
        chart.

        Parameters:
            time: (bool) if true x-axis is time, false is nts
            args: matplotlib args for 1d charts
            kwargs: matplotlib keyword arguments for 1d charts
        """

        if time:
            plt.plot(self.pdf['nts'] * self.timestep,
                     self.pdf['ncol'] / self.ncol,
                     *args, **kwargs)

        else:
            plt.plot(self.pdf['nts'],
                     self.pdf['ncol'] / self.ncol,
                     *args, **kwargs)
        plt.ylim([0, 1])

    def plot_pv(self, *args, **kwargs):
        """
        Method to plot pdf data with pore volumes (non-dimensional time)

        Parameters:
            args: matplotlib args for 1d plotting
            kwargs: matplotlib kwargs for 1d plotting
        """
        pv_factor = self.pore_volume_conversion()
        plt.plot(self.pdf['nts'] * pv_factor * self.timestep,
                 self.pdf['ncol'] / self.ncol,
                 *args, **kwargs)


class ADE(object):
    """
    Class to calculate macroscopic advection dispersion
    equation parameters for field scale model parameterization

    Parameters:
        filename: (str) ascii output file name from colloid model
        nbin: (int) number of timesteps to bin a pdf for calculation


    """
    def __init__(self, filename, nbin=1000):

        if not filename.endswith('.endpoint'):
            raise FileTypeError('<>.endpoint file must be supplied')

        reader = ASCIIReader(filename)

        self.timestep = reader.timestep
        self.resolution = reader.resolution
        self.ylen = reader.ylen
        self.ncol = float(reader.df.shape[0])
        self.uy = reader.uy
        self.pdf = None
        self.__dist_func = DistributionFunction(filename, nbin)

        self.reset_pdf(nbin)

    def __reset(self):
        self.pdf = self.__dist_func.pdf

    def reset_pdf(self, nbin):
        """
        User method to reset values based on changing
        the timestep to bin pdf values

        Parameter:
            nbin: (int) umber of timesteps to bin a pdf for calculation
        """
        self.__dist_func.reset_pdf(nbin)
        self.pdf = self.__dist_func.pdf

    def solve_jury_1991(self):
        """
        Scipy optimize method to solve least sqares
        for jury 1991. Pulse flux.
        """
        # todo: test this method! look up references for clearer examples!
        from scipy.optimize import leastsq
        a = self.ncol
        l = self.ylen * self.resolution
        t = self.pdf['nts'].as_matrix
        v = self.uy
        pdf = self.pdf['ncol'].as_matrix / self.ncol
        x0 = np.array([0., 0.])

        return leastsq(self.__jury_residuals, x0, args=(a, l, t, v, pdf))

    def __jury_residuals(self, vars, A, L, t, v, pdf):
        """
        Method to estimate residuals from jury 1991 equation
        using data

        Parameters
            vars: (np.array) [dispersivity, retardation]
            A: ncol
            l: (float) ylen
            v: (float) mean fluid_velocity
            t: (float) time
            pdf: pd.dataframe c/co of colloid pdf
        """
        return pdf - self.__jury_1991(vars, A, L, t, v)

    def __jury_1991(self, vars, A, L, t, v):
        """
        Equation for Jury 1991 calculation of Dispersivity
        and Retardation

        Parameters
            vars: (np.array) [dispersivity, retardation]
            A: ncol
            l: (float) ylen
            v: (float) mean fluid_velocity
            t: (float) time
        """
        D = vars[0]
        R = vars[1]

        eq0 = (A * L * np.sqrt(R))
        eq1 = 2 * np.sqrt(np.pi * D * t ** 3)
        eq2 = -(R * L - v * t) ** 2
        eq3 = 4 * R * D * t

        return (eq0 / eq1) * np.exp(eq2 / eq3)




class ModelPlot(object):
    """
    Class to retrieve Colloid force arrays
    and plot for data analysis.

    Parameters:
        hdf5: (str) hdf5 file name

    Properties:
        keys: method to retrieve valid data keys for hdf5 file

    Methods:
        get_data: key method to retrieve data from hdf5 file
        get_data_by_path: method to retrieve data by hdf5 file path
        plot: method to plot data by key name
    """
    def __init__(self, hdf5):
        if not hdf5.endswith('hdf') and\
                not hdf5.endswith('hdf5'):
            raise FileTypeError('hdf or hdf5 file must be supplied')

        self.__hdf = Hdf5Reader(hdf5)

    @property
    def keys(self):
        return self.__hdf.keys

    def get_data(self, key):
        """
        Get data method to view and analyze colloid
        force arrays

        Parameters:
            key: (str) valid dictionary key from self.keys

        Returns:
            data <varies>
        """
        return self.__hdf.get_data(key)

    def get_data_by_path(self, path):
        """
        Method to retrieve hdf5 data by specific path

        Parameters
            path: (str) hdf5 directory path to data

        Return
            data <varies>
        """
        return self.__hdf.get_data_by_path(path)

    def plot(self, key, *args, **kwargs):
        """
        Hdf array plotting using Hdf5Reader keys

            key: (str) valid dictionary key from self.keys
            *args: matplotlib plotting args
            **kwargs: matplotlib plotting kwargs
        """
        # todo: create a function_fmt for axis options

        if key in ('lvdw_x', 'lvdw_y',
                   'lewis_x', 'lewis_y',
                   'edl_x', 'edl_y',
                   'dlvo_x', 'dlvo_y'):

            x_axis = self.__hdf.get_data('distance_array')
            arr = self.__hdf.get_data(key)
            plt.plot(x_axis, arr, *args, **kwargs)

        elif key in ('conversion_factor',
                     'gravity',
                     'bouyancy'):
            raise KeyError('{}: key not valid for plotting'.format(key))

        else:
            plt.imshow(self.__hdf.get_data(key), *args, **kwargs)

    def plot_velocity_magnitude(self, nbin=10, *args, **kwargs):
        """
        Method to create a quiver plot to display the
        magnitude and direction of velocity vectors within
        the system.

        Parameters:
            nbin: refinement for quiver plotting
            args: matplotlib plotting args
            kwargs: matplotlib plotting kwargs
        """
        x = self.__hdf.get_data('velocity_x')
        y = self.__hdf.get_data('velocity_y')

        xx = np.arange(0, x.shape[1])
        yy = np.arange(0, x.shape[0])

        xx, yy = np.meshgrid(xx, yy)

        Q = plt.quiver(xx[::nbin, ::nbin], yy[::nbin, ::nbin],
                       x[::nbin, ::nbin], y[::nbin, ::nbin],
                       units='width', *args, **kwargs)
        qk = plt.quiverkey(Q, 0.9, 0.9, 0.01, r'$1 \frac{cm}{s}$',
                           coordinates='figure')
        plt.xlim(0, x.shape[1])
        plt.ylim(0, x.shape[0])


class CCModelPlot(object):
    """
    Class to query colloid-colloid interactions
    and plot data as 1d or as a meshgrid object
    More sophisticated than standard ModelPlot

    Parameters:
        hdf5: (str) hdf5 file name

    Properties:
        keys: returns valid dictionary keys for
            colloid plotting

    Methods:
        get_data: retrieves data by keyword
        get_data_by_path: retrieves data by hdf5 data path
        plot: plot 1d profile of dlvo curves
        plot_mesh: plot 2d dlvo profile.

    """
    data_paths = {'col_col_x': 'colloidcolloid/x',
                  'col_col_y': 'colloidcolloid/y',
                  'col_col': None,
                  'distance_x': 'colloid_colloid/distance/x',
                  'distance_y': 'colloid_colloid/distance/y',
                  'distance_fine_x': 'colloid_colloid/fine/distance/x',
                  'distance_fine_y': 'colloid_colloid/fine/distance/y',
                  'col_col_fine_x': 'colloid_colloid/fine/x',
                  'col_col_fine_y': 'colloid_colloid/fine/y',
                  'col_col_fine': None}

    def __init__(self, hdf5):
        if not hdf5.endswith('hdf') and\
                not hdf5.endswith('hdf5'):
            raise FileTypeError('hdf or hdf5 file must be supplied')

        self.__hdf5 = Hdf5Reader(hdf5)

    @property
    def keys(self):
        """
        Method to return valid keys to obtain data
        """
        return CCModelPlot.keys

    def get_data(self, key):
        """
        Method to return data by key
        """
        return self.__hdf5.get_data(key)

    def get_data_by_path(self, path):
        """
        Method to return data by hdf5 path
        """
        return self.__hdf5.get_data_by_path(path)

    def plot(self, key, *args, **kwargs):
        """
        Plotting method for 1d colloid-colloid dlvo profiles

        Parameters:
            key: (str) valid data key
            args: matplotlib plotting args
            kwargs: matplotlib plotting kwargs
        """
        # todo: store at fine discretization also for nicer plotting!!!!

        if key not in ('col_col_x', 'col_col_y',
                       'col_col_fine_x', 'col_col_fine_y'):
            raise KeyError("{} is not a valid key".format(key))

        colcol = self.__hdf5.get_data(key)
        shape = colcol.shape
        center = shape[0] // 2
        if key == "col_col_x":
            x = self.__hdf5.get_data('distance_x')
            x = x[center, center:]
            y = colcol[center, center:]

        elif key == "col_col_y":
            x = self.__hdf5.get_data('distance_y')
            x = x.T[center, center:]
            y = colcol.T[center, center:]

        elif key == "col_col_fine_x":
            x = self.__hdf5.get_data('distance_fine_x')
            x = x[center, center:]  # * 1e-6
            y = colcol[center, center:]

        else:
            x = self.__hdf5.get_data('distance_fine_y')
            x = x[center, center:]  # * 1e-6
            y = colcol[center, center:]

        plt.plot(x, y * -1, *args, **kwargs)

    def plot_mesh(self, key, *args, **kwargs):
        """
        Plotting method for 2d representation of colloid-colloid
        dlvo profiles.

        Parameters:
            key: (str) valid data key
            args: matplotlib plotting args
            kwargs:  matplotlib plotting kwargs
        """
        from matplotlib.colors import LogNorm
        if key not in ('col_col', 'col_col_fine',
                       'col_col_x', 'col_col_y',
                       'col_col_fine_x', 'col_col_fine_y'):
            raise KeyError("{} is not a valid key".format(key))

        if key == 'col_col':
            ccx = np.abs(self.__hdf5.get_data('col_col_x'))
            ccy = np.abs(self.__hdf5.get_data('col_col_y'))
            mesh = ccx + ccy

        elif key == 'col_col_fine':
            ccx = np.abs(self.__hdf5.get_data('col_col_fine_x'))
            ccy = np.abs(self.__hdf5.get_data('col_col_fine_y'))
            mesh = ccx + ccy

        else:
            mesh = np.abs(self.__hdf5.get_data(key))

        xx, yy = np.meshgrid(np.arange(0, mesh.shape[0]+1),
                             np.arange(0, mesh.shape[1] + 1))

        center = mesh.shape[0] / 2.
        plt.pcolormesh(xx, yy, mesh,
                       norm=LogNorm(vmin=mesh.min(),
                                    vmax=mesh.max()),
                       *args, **kwargs)

        plt.ylim([0, mesh.shape[0]])
        plt.xlim([0, mesh.shape[1]])
        plt.plot([center], [center], 'ko')


class ColloidVelocity(object):
    """
    Method to return colloid velocity and stats.
    relating to colloid velocity for a simulation.

    Parameters:
        filename: (str) endpoint file name

    """
    def __init__(self, filename):

        if not filename.endswith(".endpoint"):
            raise FileTypeError('.endpoint file must be supplied')

        reader = ASCIIReader(filename)
        self.timestep = reader.timestep
        self.resolution = reader.resolution
        self.xlen = reader.xlen
        self.ylen = reader.ylen
        self.df = reader.df
        self.ncol = reader.df.shape[0]
        self.max_time = max(reader.df['nts']) * self.timestep
        self.velocity = None
        self.__get_velocity_array()

    def __get_velocity_array(self):
        """
        Built in method to calculate the mean velocity of
        each colloid in the simulation
        """
        colloid = []
        velocity = []

        for index, row in self.df.iterrows():
            if np.isnan(row['y-position']):
                velocity.append((self.ylen * self.resolution) /
                                (row['delta-ts'] * self.timestep))
            else:
                velocity.append((row['y-position'] * self.resolution) /
                                (row['nts'] * self.timestep))

            colloid.append(index)

        arr = np.recarray(len(colloid,), dtype=[('colloid', np.int),
                                                ('velocity', np.float)])

        for idx, value in enumerate(colloid):
            arr[idx] = tuple([value, velocity[idx]])

        self.velocity = arr

    @property
    def max(self):
        return self.velocity['velocity'].max()

    @property
    def min(self):
        return self.velocity['velocity'].min()

    @property
    def mean(self):
        return self.velocity['velocity'].mean()

    @property
    def var(self):
        return np.var(self.velocity['velocity'])

    @property
    def stdev(self):
        return np.std(self.velocity['velocity'])

    @property
    def cv(self):
        return (self.stdev / self.mean) * 100

    def plot(self, *args, **kwargs):
        """
        Method to plot distribution of velocities by
        colloid for array of velocity.

        Parameters
        ----------
            args: matplotlib plotting args
            kwargs: matplotlib plotting kwargs
        """
        plt.plot(self.velocity['colloid'],
                 self.velocity['velocity'],
                 *args, **kwargs)

    def plot_histogram(self, nbin=10, width=0.01,
                       *args, **kwargs):
        """
        User method to plot a histogram of velocities

        Parameters:
            nbin: (int) number of specific bins for plotting
            width: (float) matplotlib bar width.
        """

        adjuster = 0.00001
        bins = np.linspace(self.min - adjuster, self.max, nbin)
        ncols = []
        velocity = []
        lower_v = self.min - adjuster
        upper_v = 0

        for upper_v in bins:
            ncol = 0
            for v in self.velocity['velocity']:
                if lower_v < v <= upper_v:
                    ncol += 1

            velocity.append((lower_v + upper_v)/2.)
            ncols.append(ncol)
            lower_v = upper_v - adjuster

        velocity.append(upper_v + adjuster)
        ncols.append(0)

        plt.bar(velocity, ncols, width, *args, **kwargs)


# todo: think about this one. Does it belong here?
class LBOutput(object):
    """
    Class to anaylze LB fluid/solid properties

    Parameters:
        hdf: (str) hdf5 output filename

    Attributes:

    Parameters:


    """
    data_paths = {'velocity_x': None,
                  'velocity_y': None,
                  'resolution': None,
                  }

    def __init__(self, hdf5):
        if not hdf5.endswith('.hdf') and not\
                hdf5.endswith('.hdf5'):
            raise FileTypeError('hdf or hdf5 file must be supplied')

        self.__hdf5 = Hdf5Reader(hdf5)

    @property
    def keys(self):
        return LBOutput.data_paths.keys()



class ASCIIReader(object):
    """
    Class to read in text based output files to a
    pandas dataframe

    Parameters:
        filename: (str) output filename (ie. endpoint, timestep, or pathline)

    """

    # todo: use a comment block to and add continuos/pulse and possibly ncol
    # todo: to the header

    dtypes = {'colloid': np.int,
              'flag': np.int,
              'nts': np.int,
              'x-position': np.float,
              'y-position': np.float,
              'x-model': np.float,
              'y-model': np.float,
              'start-ts': np.int,
              'end-ts': np.int,
              'delta-ts': np.int}

    def __init__(self, filename):
        self.timestep = 0
        self.resolution = 0
        self.xlen = 0
        self.ylen = 0
        self.ux = 0
        self.uy = 0
        self.__data_startline = 0
        self.__header = []

        if filename.split('.')[-1] not in ('endpoint', 'timeseries', 'pathline'):
            raise FileTypeError("{}: not in supported filetypes".format(filename))

        else:
            self.read_header(filename)
            self.df = self.read_ascii(filename)

    def read_header(self, filename):
        """
        Method to read the header from ascii output files for LB-Colloids

        Parameters:
            filename: (str) output filename (ie. endpoint, timestep, or pathline)
        """
        with open(filename) as f:
            for idx, line in enumerate(f):
                if line.startswith("Timestep"):
                    t = line.split()
                    self.timestep = float(t[-1].rstrip())

                elif line.startswith('Resolution'):
                    t = line.split()
                    self.resolution = float(t[-1].rstrip())

                elif line.startswith('xlen'):
                    t = line.split()
                    self.xlen = float(t[-1].rstrip())

                elif line.startswith('ylen'):
                    t = line.split()
                    self.ylen = float(t[-1].rstrip())

                elif line.startswith('ux'):
                    t = line.split()
                    self.ux = float(t[-1].rstrip())

                elif line.startswith('uy'):
                    t = line.split()
                    self.uy = float(t[-1].rstrip())

                elif line.startswith("#"*10):
                    self.__data_startline = idx + 1
                    break

                else:
                    pass

    def read_ascii(self, filename):
        """
        Method to read endpoint file data from from ascii files for LB-Colloids
        Sets data to pandas dataframe

        Parameters:
            filename: (str) output filename (ie. endpoint, timestep, or pathline)
        """
        with open(filename) as f:
            t = []
            for idx, line in enumerate(f):
                if idx < self.__data_startline:
                    pass
                elif idx == self.__data_startline:
                    self.__header = [i.rstrip() for i in line.split()
                                     if i not in ('\t', '', ' ', '\n')]
                else:
                    t.append([self.__try_float(i.rstrip()) for i
                              in line.split() if i not in ('\t', '', ' ', '\n')])

        temp = np.array(t).T

        temp = {self.__header[idx]: data for idx, data in enumerate(temp)}
        df = pd.DataFrame(temp)
        df = df.reindex_axis(self.__header, axis=1)
        df = df.set_index('colloid')
        return df

    @staticmethod
    def __try_float(val):
        try:
            return float(val)
        except ValueError:
            return float('nan')


class Hdf5Reader(object):
    """
    Reader object to read in HDF5 stored outputs
    from colloid models.

    Parameters:
        hdf5 (str) LB-Colloid hdf5 file name

    Methods:
        keys: method to retrieve valid data keys for hdf5 file
        get_data: key method to retrieve data from hdf5 file
        get_data_by_path: method to retrieve data by hdf5 file path
    """
    data_paths = {'image': 'Binary_image',
                  'lb_velocity_x': 'results/uarray',
                  'lb_velcoity_y': 'results/uarray',
                  'conversion_factor': 'results/velocity_factor',
                  'brownian_x': 'colloids/brownian/x',
                  'brownian_y': 'colloids/brownian/y',
                  'lvdw_x': 'colloids/lvdw/x',
                  'lvdw_y': 'colloids/lvdw/y',
                  'edl_x': 'colloids/edl/x',
                  'edl_y': 'colloids/edl/y',
                  'lewis_x': 'colloids/lewis_acid_base/x',
                  'lewis_y': 'colloids/lewis_acid_base/y',
                  'velocity_x': 'colloids/ux',
                  'velocity_y': 'colloids/uy',
                  'gravity': 'colloids/gravity',
                  'bouyancy': 'colloids/bouyancy',
                  'distance_array': 'colloids/distance_arr',
                  'dlvo_x': None,
                  'dlvo_y': None,
                  'col_col_x': 'colloid_colloid/x',
                  'col_col_y': 'colloid_colloid/y',
                  'col_col': None,
                  'distance_x': 'colloid_colloid/distance/x',
                  'distance_y': 'colloid_colloid/distance/y',
                  'distance_fine_x': 'colloid_colloid/fine/distance/x',
                  'distance_fine_y': 'colloid_colloid/fine/distance/y',
                  'col_col_fine_x': 'colloid_colloid/fine/x',
                  'col_col_fine_y': 'colloid_colloid/fine/y',
                  'col_col_fine': None}

    def __init__(self, hdf5):
        if not hdf5.endswith('hdf') and\
                not hdf5.endswith('hdf5'):
            raise FileTypeError('hdf or hdf5 file must be supplied')

        self.file_name = hdf5

    @property
    def keys(self):
        return [i for i in Hdf5Reader.data_paths]

    def get_data(self, key):
        """
        Method to retrieve hdf5 data by dict. key

        Parameters:
            key: (str) valid dictionary key from self.keys

        Returns:
            data <varies>
        """
        if key not in Hdf5Reader.data_paths:
            raise KeyError('Dictionary key not in valid keys. Use get_data_by_path')

        hdf = H.File(self.file_name, 'r')

        if key == 'lb_velocity_x':
            data = hdf[Hdf5Reader.data_paths[key]][()][1]

        elif key == 'lb_velocity_y':
            data = hdf[Hdf5Reader.data_paths[key]][()][0]

        elif key == 'dlvo_x':
            data = hdf[Hdf5Reader.data_paths['edl_x']][()] +\
                hdf[Hdf5Reader.data_paths['lewis_x']][()] +\
                hdf[Hdf5Reader.data_paths['lvdw_x']][()]
            data = data[0]

        elif key == 'dlvo_y':
            data = hdf[Hdf5Reader.data_paths['edl_y']][()] +\
                hdf[Hdf5Reader.data_paths['lewis_y']][()] +\
                hdf[Hdf5Reader.data_paths['lvdw_y']][()]
            data = data[0]

        elif key in ('lvdw_x', 'lvdw_y',
                     'lewis_x', 'lewis_y',
                     'edl_x', 'edl_y',
                     'dlvo_x', 'dlvo_y',
                     'distance_array'):

            data = hdf[Hdf5Reader.data_paths[key]][()][0]

        else:
            data = hdf[Hdf5Reader.data_paths[key]][()]

        hdf.close()
        return data

    def get_data_by_path(self, path):
        """
        Method to retrieve hdf5 data by specific path

        Parameters
            path: (str) hdf5 directory path to data

        Return
            data <varies>
        """
        hdf = H.File(self.file_name, 'r')
        data = hdf[path][()]
        hdf.close()
        return data


class FileTypeError(Exception):
    pass
