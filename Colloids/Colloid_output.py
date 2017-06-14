import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# todo: return an axis object and dump the data into csv for further processing?


class Breakthrough(object):
    """
    Class to prepare and plot breakthrough curve data from endpoint
    and possibly timeseries files

    """
    def __init__(self, filename):
        pass


class ADE(object):
    pass


class DistributionFunction(object):
    pass


class DLVOPlot(object):
    pass


class CCDLVOPlot(object):
    pass


class CCDLVOMesh(object):
    pass


class ColloidVelocity(object):
    pass


class ASCIIReader(object):
    """
    Convience method to read in text based output files to a
    pandas dataframe

    Parameters:
        filename: (str) output filename (ie. endpoint, timestep, or pathline)

    """

    # todo: use a comment block to deliniate the ts value and add continuos/pulse
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

        # todo: convert this to a structured array or pandas dataframe
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


class FileTypeError(Exception):
    pass
