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

            df = pd.DataFrame({'ts': nts, 'ncol': ncols}).set_index('ncol')

            self.__breakthrough_curve = df

        return self.__breakthrough_curve

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
            plt.plot(self.breakthrough_curve['ts'] * self.timestep,
                     self.breakthrough_curve.index.values / self.ncol,
                     *args, **kwargs)

        else:
            plt.plot(self.breakthrough_curve['ts'],
                     self.breakthrough_curve.index.values / self.ncol,
                     *args, **kwargs)
        plt.ylim([0, 1])


class DistributionFunction(object):
    """
    Class to plot a pdf function of colloid breakthrough
    from endpoint files.

    Parameters:
        filename: (str)

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
    def __init__(self, filename, bin=1000):
        if not filename.endswith('.endpoint'):
            raise FileTypeError('.endpoint file must be supplied')
        reader = ASCIIReader(filename)
        self.df = reader.df
        self.resolution = reader.resolution
        self.timestep = reader.timestep
        # todo: replace this call with something from the header later!
        self.ncol = float(self.df.shape[0])
        self.bin = bin
        self.pdf = None
        self.reset_pdf(bin)

    def reset_pdf(self, bin):
        """
        Method to generate a probability distribution function
        based upon user supplied bin size.

        Parameters:
            bin: (int) number of time steps to base bin on

        Returns:
            () probability distribution function of colloid
            breakthrough.
        """
        self.bin = bin
        ts = []
        ncols = []
        lower_nts = 0
        max_ts = self.df['nts'].max()
        pdf_colloids = self.df.loc[self.df['flag'] == 3]
        pdf_colloids = pdf_colloids.sort_values('delta-ts')

        for upper_nts in range(0, int(max_ts) + 1, bin):
            ncol = 0
            for index, row in pdf_colloids.iterrows():
                x = row['delta-ts']
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


class ADE(object):
    pass


class DLVOPlot(object):
    pass


class CCDLVOPlot(object):
    pass


class CCDLVOMesh(object):
    pass


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
                velocity.append((self.ylen * self.resolution)/
                                (row['delta-ts'] * self.timestep))
            else:
                velocity.append((row['y-position'] * self.resolution)/
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
        bar_width = 0.01
        bins = np.linspace(self.min - adjuster, self.max, nbin)
        ncols = []
        velocity = []
        lower_v= self.min - adjuster
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


class ASCIIReader(object):
    """
    Convience method to read in text based output files to a
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

    def __try_float(self, val):
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
                  'bouyancy': 'colloids/bounancy'}

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
