from lb_colloids import ColloidOutput
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class OutputProcess(object):
    """
    Simple class to process outputs from PhD validation project
    """
    def __init__(self, filename, ncol_per_ml=0, nts_per_ml=0, ncol_per_second=0):
        t = ColloidOutput.ASCIIReader(filename)
        self.name = filename
        self.dataframe = t.df
        self.ux = t.ux
        self.uy = t.uy
        self.continuous = t.continuous
        self.ncol = t.ncol
        self.velocity_factor = t.velocity_factor
        self.resolution = t.resolution
        self.timestep = t.timestep
        self.ylen = t.ylen
        self.xlen = t.xlen
        self.ncol_per_ml = float(ncol_per_ml)
        self.nts_per_ml = float(nts_per_ml)
        self.ncol_per_second = float(ncol_per_second)
        self.nts_per_second = 1. / self.timestep

    @property
    def get_pore_volume_factor(self):
        pv_factor = (abs(self.uy) * self.velocity_factor) / \
                    (self.ylen * self.resolution)
        return pv_factor

    @property
    def get_ml_data(self):
        if not self.ncol_per_ml:
            raise AssertionError("Must generate object with a value for ncol_per_ml")

        if not self.nts_per_ml:
            raise AssertionError("Must generate object with a value for nts_per_ml")

        # todo: Do stuff!

        return


    @property
    def seconds(self):
        """
        Generates a np.ndarray of seconds for the simulation
        :return: np.ndarray
        """
        niters = max(self.dataframe['nts'])
        n_seconds = np.ceil(niters * self.timestep)

        return np.arange(0, n_seconds + 1)

    @property
    def pore_volumes(self):
        """
        Generates a np.ndarray of pore volumes elution per simulation second
        :return:
        """
        return self.seconds * self.get_pore_volume_factor

    @property
    def ncol_release_per_second(self):
        """
        Generates a np.ndarray of number of colloids released from model
        each second
        :return: np.ndarray
        """
        ncol = []
        nsec = np.array(self.dataframe['nts'] * self.timestep)
        flag = np.array(self.dataframe['flag'])

        for sec in self.seconds:
            n = 0
            for ix, val in enumerate(nsec):
                if flag[ix] == 3 and sec - 1 < val <= sec:
                    n += 1
            ncol.append(n)

        return np.array(ncol)

    @property
    def relative_concentration_per_second(self):
        """
        Generates a np.ndarray of relative concentration of colloids
        released from model each second
        :return: np.ndarray
        """
        if not self.ncol_per_second:
            print("Warning: setting ncol_per_second to continuous value")
            self.ncol_per_second = self.continuous

        if not self.ncol_per_second:
            raise AssertionError("A value for ncol_per_second must be supplied")

        return self.ncol_release_per_second / self.ncol_per_second

    @property
    def get_second_data(self):
        """
        Gets a dictionary of data by the second for further processing
        and visualization of model results
        :return: dict (keys are 'seconds', 'ncol', 'c_c0')
        """
        if not self.ncol_per_second:
            print("Warning: setting ncol_per_second to continuous value")
            self.ncol_per_second = self.continuous

        if not self.ncol_per_second:
            raise AssertionError("A value for ncol_per_second must be supplied")

        seconds = self.seconds
        ncol = self.ncol_release_per_second
        c_c0 = self.relative_concentration_per_second

        return {'seconds': seconds, 'ncol': ncol, 'c_c0': c_c0}

    def plot(self, *args, **kwargs):
        return


if __name__ == "__main__":
    oname = "SyVal200p_2_1e-05.endpoint"

    test = OutputProcess(oname, ncol_per_second=71)
    x = test.relative_concentration_per_second
    s = test.seconds
    pv = test.pore_volumes

    plt.plot(pv, x, "bo-")
    plt.ylim([0, 2.0])
    plt.xlim([0, max(pv)])
    plt.show()


