"""
Basic Input and Output control for Colloid Simulation models are hosted within
Colloid_IO.py. Config and ColloidConfig are classes to set up
model dictionaries can be passed along to the main simulation routines.
ColloidConfig is a backend method for the super user. This provides a simple overridden
dictionary class that is able to build configuration files and a list object that
can be passed directly to the Config class.

Importing classes from this module follows the notation:

>>> from lb_colloids import cIO
>>>
>>> config = cIO.Config("Colloids.config")
>>> cc = ColloidsConfig()
>>> cc['I'] = 0.01
>>> # user must supply all required parameters to the ColloidsConfig dictionary, or an AssertionError will be raised
>>> cc_config = cc.config
>>> config = cIO.Config(cc_config)
"""
import sys
import numpy as np
import Colloid_Math as cm
import h5py as H


class Config:
    """
    Class to open and parse configuration files for LB-Colloids. Many data checks
    have been implemented to look for consistancy in data type and
    configuration variable for each input block.

    Parameters:
    ----------
    :param str fname: Configuration file name ex. Model.config. This class can also accepts
    a list of configuration variables that have been set up by cIO.ColloidsConfig()

    Returns:
    -------
    :return: model_parameters (dict) Dictionary of necessary model parameters
    :return: physical_parameters (dict) Dictionary of optional physical parameters
    :return: chemical_parmaeters (dict) Dictionary of chemical parameter options
    :return: output_control (dict) Dictionary of output control options

    """
    def __init__(self, fname):

        if isinstance(fname, str):
            self.config = self._reader(fname)
        elif isinstance(fname, list):
            self.config = [line.strip('\n') for line in fname]
        elif fname is None:
            pass
        else:
            raise TypeError('input data not a list or file name')

        self._strtype = ('LBMODEL', 'ENDPOINT', 'PATHLINE', 'TIMESERIES')
        self._inttype = ('NCOLS', 'ITERS', 'STORE_TIME', 'PRINT_TIME',
                         'CONTINUOUS', 'COL_COL_UPDATE')
        self._floattype = ('LBRES', 'GRIDREF', 'AC', 'TIMESTEP', 'TEMPERATURE',
                           'RHO_WATER', 'RHO_COLLOID', 'VISCOSITY', 'I_INITIAL',
                           'I', 'EPSILON_R', 'SHEER_PLANE', 'LVDWST_WATER',
                           'LVDWST_COLLOID', 'LVDWST_SOLID', 'ZETA_COLLOID',
                           'ZETA_SOLID', 'PSI+_COLLOID', 'PSI+_WATER', 'PSI+_SOLID',
                           'PSI-_COLLOID', 'PSI-_WATER', 'PSI-_SOLID', 'SCALE_LB')
        self._booltype = ('PLOT', 'ADJUST_ZETA', 'OVERWRITE', 'SHOWFIG')
        self._dicttype = ('CONCENTRATION', 'VALENCE')
        self._required = ('LBMODEL', 'NCOLS', 'ITERS', 'LBRES', 'GRIDREF',
                          'TS')
        self.validmodelparams = ('LBMODEL', 'NCOLS', 'ITERS', 'LBRES',
                                 'GRIDREF', 'AC', 'TIMESTEP', 'TEMPERATURE',
                                 'RHO_COLLOID', 'CONTINUOUS', 'COL_COL_UPDATE')
        self.validphysicalparams = ('RHO_WATER', 'RHO_COLLOID', 'VISCOSITY', "SCALE_LB")
        self.validchemicalparams = ('CONCENTRATION', 'ADJUST_ZETA', 'I_INITIAL',
                                    'I', 'EPSILON_R', 'VALENCE', 'SHEER_PLANE',
                                    'LVDWST_WATER', 'LVDWST_COLLOID', 'LVDWST_SOLID',
                                    'ZETA_COLLOID', 'ZETA_SOLID', 'PSI+_COLLOID',
                                    'PSI+_WATER', 'PSI+_SOLID', 'PSI-_COLLOID',
                                    'PSI-_WATER', 'PSI-_SOLID', 'RHO_COLLOID')
        self.validoutputparams = ('STORE_TIME', 'PRINT_TIME', 'PLOT',
                                  'ENDPOINT', 'PATHLINE', 'TIMESERIES',
                                  'OVERWRITE', 'PLOT', 'SHOWFIG')

    def _reader(self, fname):
        """
        reads in the input config file
        """
        with open(fname, 'r') as f:
            config = [line.strip('\n').strip(' ') for line in f]
        return config

    def model_parameters(self):
        """
        Reads the MODEL PARAMETERS block of the configuration file and creates
        the model_dict which contains the required parameters to run LB-Colloid
        """
        model_dict = {}
        blockname = 'MODEL PARAMETERS'

        modelparams = self.get_block(blockname)
        for parameter in modelparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validmodelparams)
            pname = self.adjust_pname(pname)
            model_dict[pname] = param
        self.check_model_parameters(model_dict)
        return model_dict

    def physical_parameters(self):
        """
        Reads the PHYSICAL PARAMETERS block of the configuration file and creates
        the physics_dict that passes optional parameters to the physical force calculations
        in the Colloid_Math.py module
        """
        physics_dict = {}
        blockname = 'PHYSICAL PARAMETERS'

        physicalparams = self.get_block(blockname)
        for parameter in physicalparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validphysicalparams)
            pname = self.adjust_pname(pname)
            physics_dict[pname] = param

        physics_dict = self.add_universal_parameters(physics_dict)
            
        return physics_dict

    def chemical_parameters(self):
        """
        Reads the CHEMICAL PARAMETERS block of the configuration file and creates
        the chemical_dict that passes optional parameters to the chemical force calculations
        in the Colloid_Math.py module
        """
        chemical_dict = {}
        blockname = 'CHEMICAL PARAMETERS'

        chemicalparams = self.get_block(blockname)
        for parameter in chemicalparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validchemicalparams)
            pname = self.adjust_pname(pname)
            chemical_dict[pname] = param

        chemical_dict = self.add_universal_parameters(chemical_dict)
        return chemical_dict

    def output_control(self):
        """
        Reads the OUTPUT CONTROL block of the configuration file and creates
        the output_dict that passes optional parameters to control model output
        """
        output_dict = {}
        blockname = 'OUTPUT CONTROL'

        outputparams = self.get_block(blockname)
        for parameter in outputparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validoutputparams)
            pname = self.adjust_pname(pname)
            output_dict[pname] = param
        
        output_dict = self.add_universal_parameters(output_dict)
        return output_dict

    def get_block(self, blockname):
        """
        Method to isolate an input block from a configuration file for parsing

        Parameters:
        ----------
        :param str blockname: Blockname of an input block in the configuration file.

        Returns:
        -------
        :return: (list) returns a list of parameters contained within the block
        """
        try:
            idx0 = self.config.index('START %s' % blockname)
            idx1 = self.config.index('END %s' % blockname)
        except ValueError:
            print('Please define a %s block' % blockname)
            sys.exit(-1)
        return self.config[idx0+1:idx1]

    def parametertype(self, parameter):
        """
        Method takes a parameter string, splits it, and sets the parameter type

        Parameters:
        ----------
        :param str parameter: String of "<pname>: <param>"

        returns:
        -------
        :return: (pname, param) (tuple) parameter name, value
        """
        # split the paramter block line, and strip whitespace
        pname, param = parameter.split(':')
        pname = pname.strip(' ').upper()
        param = param.strip(' ')

        if pname in self._strtype:
            param = str(param)
            
        elif pname in self._inttype:
            param = int(param)

        elif pname in self._floattype:
            param = float(param)

        elif pname in self._booltype:
            if param.upper() == 'TRUE':
                param = True
            else:
                param = False

        elif pname in self._dicttype:
            keyvalue = [i for i in param.split(' ') if i != ' ']
            param = {}
            for i in range(0, len(keyvalue), 2):
                param[keyvalue[i]] = float(keyvalue[i+1])              
                
        else:
            print 'Parameter name %s is not valid' % pname
            sys.exit(-1)

        return pname, param

    def check_if_valid(self, blockname, pname, validparams):
        """
        Method checks if a specific parameter is valid for the block it was supplied to
        in the configuration file.

        Parameters:
        ----------
        :param str blockname: Name of the input block
        :param str pname:
        :param tuple validparams: tuple of valid parameter names for the specific input block
        """
        if pname in validparams:
            return
        else:
            print('%s is not a valid parameter for %s' % (pname, blockname))
            print(validparams)
            sys.exit(-1)

    def adjust_pname(self, pname):
        """
        Adjusts parameter name from configuration file name to LB-Colloids name for a limited number of
        parameters. Sets all other parameters as lowercase to follow PEP-8

        :param str pname: parameter name
        """
        if pname in ('I_INITIAL', 'I', 'TEMPERATURE', 'TIMESTEP'):
            if pname == 'I_INITIAL':
                return 'I_initial'
            elif pname == 'I':
                return 'I'
            elif pname == 'TEMPERATURE':
                return 'T'
            elif pname == 'TIMESTEP':
                return 'ts'
            else:
                pass
        else:
            return pname.lower()

    def add_universal_parameters(self, Dict):
        """
        Add common model parameters to other dictionaries if present. Necessary for parameterization by
        kwargs of physics and chemistry.

        :param dict Dict: dictionary to add set of universal paramameters from the required parameters

        :return: Dict (dict)
        """
        modelparams = self.model_parameters()
        for key in modelparams:
            if key in ('ac', 'T', 'visocity', 'ts', 'lbres', 'gridref', 'ncols'):
                Dict[key] = modelparams[key]
        return Dict

    def check_model_parameters(self, ModelDict):
        """
        Check for required parameters. If not present inform the user which parameter(s) are needed

        :param dict ModelDict: Checks the model_dict for all required parameters.
        :raise AssertionError: If all required parameters are not supplied

        """
        for key in self._required:
            if key.lower() not in ModelDict:
                if key.upper() not in ModelDict:
                    if key.lower() == 'ts':
                        key = 'TIMESTEP'
                    else:
                        raise AssertionError('%s is a required MODEL PARAMETER' % key)
                else:
                    pass
            else:
                pass


class Output:
    """
    Output class for writing formatted ASCII files. This class generates
    <.endpoint>, <.timeseries> and <.pathline> files. All keywords are required.

    Parameters:
    ----------
    :param str fi: filename of the output.
    :keyword bool overwrite: Determines if file is overwritten, or appended to.
        useful for generating new pathline and timeseries files
    :keyword float ts: Physical time step
    :keyword float lbres: Lattice Boltzmann resolution
    :keyword float gridref: Grid refinement factor
    :keyword int ncols: number of colloids simulated
    :keyword int xlen: length of the xdomain in pixels
    :keyword int ylen: length of the ydomain in pixels
    :keyword float mean_ux: mean fluid velocity in x-direction
    :keyword float mean_uy: mean fluid velocity in y-direction
    :keyword int continuous: flag for continuous release of colloids
    """
    def __init__(self, fi, **kwargs):

        defaults = {'overwrite': True,
                    'ts': 1e-6}

        for key in kwargs:
            defaults[key] = kwargs[key]

        self.filename = fi
        self.resolution = defaults['lbres'] / defaults['gridref']
        self.header = "LB Colloids output file\nTimestep: {}\n".format(defaults['ts'])
        self.header += "Ncols: {}\n".format(defaults['ncols'])
        self.header += "Resolution: {}\nxlen: {}\n".format(self.resolution, defaults['xlen'])
        self.header += "ylen: {}\nux: {}\nuy: {}\n".format(defaults['ylen'],
                                                             defaults['mean_ux'],
                                                             defaults['mean_uy'])
        self.header += "velocity_factor: {}\n".format(defaults['velocity_factor'])
        self.header += "Continuous: {}\n\n".format(defaults['continuous'])
        self.header += "#" * 136 + "\n"
        self.header += '{:>8}\t{:>5}\t{:>5}\t{:>11}\t{:>12}\t' \
                       '{:>11}\t{:>10}\t{:>10}\t{:>10}\t{:>10}\n'.format(
                       'colloid', 'flag', 'nts', 'x-position', 'y-position',
                       'x-model', 'y-model', 'start-ts', 'end-ts',
                       'delta-ts')

        # checks if file exists, and overwrites if overwrite is True!
        if defaults['overwrite'] is True:
            self._writer(fi, self.header, wtype='w')
        
    def _writer(self, fi, data, wtype='a'):
        with open(fi, wtype) as f:
            for line in data:  # format as 1d list that is str fomated with \n
                f.write(line)

    def write_output(self, timer, colloids, pathline=True):
        """
        Set up and write colloid streaming output to an ASCII file

        Parameters:
        ----------
        :param TrackTime timer: Model TrackTime instance
        :param list colloids: Colloid simulation list containing LB_Colloid.Colloid objects
        :param bool pathline: Flag to indicate if an endpoint or pathline/timeseries is being written.
        """

        time = timer.timer
        output = []
        
        if pathline is not True:
            for number, colloid in enumerate(colloids):
                output_string = '{:>8d}\t{:>5d}\t{:5d}\t{:09.8f}\t' \
                                '{:10.9f}\t ' \
                                '{:09.5f}\t{:09.5f}\t{:6f}\t{:6f}\t{:6f}\n'.format(
                                    colloid.tag, colloid.flag[-1], time[-1],
                                    colloid.xposition[-1],
                                    colloid.yposition[-1],
                                    colloid.xposition[-1]/self.resolution,
                                    colloid.yposition[-1]/self.resolution,
                                    colloid.colloid_start_time,
                                    colloid.colloid_end_time,
                                    colloid.colloid_end_time - colloid.colloid_start_time)
                output.append(output_string)
        
        else:    
            for idx in range(len(time)-1):  # so we don't duplicate
                for number, colloid in enumerate(colloids):
                    output_string = '{:>8d}\t{:>5d}\t{:5d}\t{:09.8f}\t{:09.8f}\t' \
                                    '{:09.8f}\t{:09.8f}\t{:7f}\t{:7f}\t{:7f}\n'.format(
                                        colloid.tag, colloid.flag[-1], time[-1],
                                        colloid.xposition[-1],
                                        colloid.yposition[-1],
                                        colloid.xposition[-1]/self.resolution,
                                        colloid.yposition[-1]/self.resolution,
                                        colloid.colloid_start_time,
                                        colloid.colloid_end_time,
                                        colloid.colloid_end_time - colloid.colloid_start_time)
                    output.append(output_string)

        self._writer(self.filename, output)

    def write_single_colloid(self, timer, colloid):
        """
        Method to write a single colloid to an endpoint file upon breakthrough
        :param TrackTime timer: Model TrackTime instance
        :param LB_Colloid.Colloid colloid:
        """
        time = timer.timer
        output_string = '{:>8d}\t{:>5d}\t{:5d}\t{:09.8f}\t' \
                        '{:10.9f}\t ' \
                        '{:09.5f}\t{:09.5f}\t{:6f}\t{:6f}\t{:6f}\n'.format(
            colloid.tag, colloid.flag[-1], time[-1],
            colloid.xposition[-1],
            colloid.yposition[-1],
            colloid.xposition[-1] / self.resolution,
            colloid.yposition[-1] / self.resolution,
            colloid.colloid_start_time,
            colloid.colloid_end_time,
            colloid.colloid_end_time - colloid.colloid_start_time)

        self._writer(self.filename, [output_string])


class ColloidsConfig(dict):
    """
    OO class to build config files, or build a list that the Config class will recognize and parse
    Recomended setup method for the super user who is looping many models.
    Facilitates easy sensitivity analysis, etc....

    Class uses a dictionary override to set parameters to the class, and writes them out as a list
    or as a configuration file

    example class usage:

    >>> from lb_colloids import cIO
    >>> cconfig = cIO.ColloidsConfig()
    >>> cconfig['I'] = 0.1
    >>> cconfig['ncols'] = 500
    >>> x = cconfig.config  # returns a formatted list that imitates a colloids configuration file
    >>> config = cIO.Config(x)
    """
    def __init__(self):
        self.__formats = Config(None)
        self.__config = []
        super(ColloidsConfig, self).__init__()

    def __setitem__(self, key, value):
        """
        Override method that does an initial parameter check before setting to the dictionary
        """
        parameter = '{}: {}'.format(key, value)
        key, value = self.__formats.parametertype(parameter)
        super(ColloidsConfig, self).__setitem__(key, value)

    def __getitem__(self, key):
        """
        Override method that assures the key is in uppercase notation:
        """
        key = key.upper()
        return super(ColloidsConfig, self).__getitem__(key)

    @property
    def config(self):
        """
        Property method that creates the config list on the fly for the user from the overridden dictionary

        Returns:
        -------
        :return: self.__config (list) formatted configuration file list
        """

        # check for all required parameters before creating configuration list
        self.__formats.check_model_parameters(self)
        self.__format_config(self.valid_model_parameters, "model parameters")
        self.__format_config(self.valid_physical_parameters, 'physical parameters')
        self.__format_config(self.valid_chemical_parameters, 'chemical parameters')
        self.__format_config(self.valid_output_control_parameters, 'output control')
        return self.__config

    def __format_config(self, valid_parameters, block_name):
        """
        Hidden method to easily loop through the Colloids config dict and create a config
        list with proper formatting

        Parameters
        ----------
            valid_parameters: (tuple) tuple of valid parameters for the model block
            block_name: (str) configuration file block name
        """
        self.__config.append('START {}\n'.format(block_name.upper()))
        for key in self:
            if key in valid_parameters:
                self.__config.append('{}: {}\n'.format(key, self[key]))
        self.__config.append('END {}\n\n'.format(block_name.upper()))

    def __get_parameters(self, valid_parameters):
        parameter_dict = {}
        for key in self:
            if key in valid_parameters:
                parameter_dict[key] = self[key]
        return parameter_dict

    @property
    def valid_model_parameters(self):
        """
        List of valid model parameters
        """
        return self.__formats.validmodelparams

    @property
    def valid_chemical_parameters(self):
        """
        List of valid chemical parameters
        """
        return self.__formats.validchemicalparams

    @property
    def valid_physical_parameters(self):
        """
        List of valid physical parameters
        """
        return self.__formats.validphysicalparams

    @property
    def valid_output_control_parameters(self):
        """
        List of valid output control parameters
        """
        return self.__formats.validoutputparams

    @property
    def model_parameters(self):
        """
        Current user supplied model parameters
        """
        return self.__get_parameters(self.__formats.validmodelparams)

    @property
    def chemical_parameters(self):
        """
        Current user supplied chemical parameters
        """
        return self.__get_parameters(self.__formats.validchemicalparams)

    @property
    def physical_parameters(self):
        """
        Current user supplied physical parameters
        """
        return self.__get_parameters(self.__formats.validphysicalparams)

    @property
    def output_control_parameters(self):
        """
        Current user supplied output control parameters
        """
        return self.__get_parameters(self.__formats.validoutputparams)

    def write(self, fname):
        """
        Writes a configuration file with user supplied parameters to file.

        :param str fname: Configuration file name to write
        """
        with open(fname, "w") as f:
            f.writelines(self.config)


class HDF5WriteArray(object):
    """
    Class to write chemical and physical force arrays to the Model HDF5 object
    for later use in data processing and analysis.

    Parameters:
    ----------
    :param np.ndarray ux: Dimensionalized fluid velocity in the x-direction
    :param np.ndarray uy: Dimensionalized fluid velocity in the y-direction
    :param Colloid_Math.ColloidColloid colloidcolloid: ColloidColloid object
    :param dict model_dict: The supplied model dict from parameterization
    :param dict chemical_dict: The supplied chemical dict used for parameterization
    :param dict physical_dict: The supplied physical dict used for parameterization

    """
    def __init__(self, ux, uy, colloidcolloid,
                 model_dict, chemical_dict, physical_dict):

        self.__model = model_dict['lbmodel']
        arr = np.array([np.arange(1, 101) * model_dict['lbres'] / model_dict['gridref']])
        varr = np.array([np.ones(100)])
        gravity = cm.Gravity(**model_dict)
        bounancy = cm.Bouyancy(**model_dict)
        cf = cm.Gap(arr, arr, **model_dict)
        brownian = cm.Brownian(cf.f1, cf.f4, **physical_dict)

        # pop and store the existing vector profiles until dummy data has been
        # used, then add the model profile back to the chemical dict for writing
        xvArr = chemical_dict.pop('xvArr')
        yvArr = chemical_dict.pop('yvArr')
        dlvo = cm.DLVO(arr, arr, xvArr=varr, yvArr=varr, **chemical_dict)
        # set up a fine colloid surface array for plotting
        fine_arr = np.array([np.linspace(1e-9, 1e-7, 1000)])
        varr = np.array([np.ones(1000)])
        fine_dlvo = cm.DLVO(fine_arr, fine_arr, xvArr=varr, yvArr=varr, **chemical_dict)
        chemical_dict['xvArr'] = xvArr
        chemical_dict['yvArr'] = yvArr
        res = chemical_dict['lbres']

        # set up fine grained output for col-col
        fine_res = 1e-9 * chemical_dict['gridref']
        chemical_dict['lbres'] = fine_res
        cc_fine = cm.ColloidColloid(np.ones((1,1)), **chemical_dict)

        chemical_dict['lbres'] = res

        with H.File(self.__model, 'r+') as h:
            for key, value in model_dict.items():
                h.create_dataset('colloids/model_dict/{}'.format(key), data=value)
            for key, value in chemical_dict.items():
                if key in ('concentration', 'valence'):
                    for species, species_value in value.items():
                        h.create_dataset('colloids/chemical_dict/{}/{}'
                                         .format(key, species), data=species_value)
                else:
                    h.create_dataset('colloids/chemical_dict/{}'.format(key), data=value)
            for key, value in physical_dict.items():
                h.create_dataset('colloids/physical_dict/{}'.format(key), data=value)
            h.create_dataset('colloids/mock_arr', data=arr)
            h.create_dataset('colloids/gravity', data=gravity.gravity)
            h.create_dataset('colloids/bouyancy', data=bounancy.bouyancy)

            cf_dict = {'f1': cf.f1, 'f2': cf.f2, 'f3': cf.f3, 'f4': cf.f4}
            for key, value in cf_dict.items():
                h.create_dataset('colloids/cf/{}'.format(key), data=value)

            dlvo_profiles = {'brownian/x': brownian.brownian_x,
                             'brownian/y': brownian.brownian_y,
                             'edl/x': dlvo.EDLx,
                             'edl/y': dlvo.EDLy,
                             'lewis_acid_base/x': dlvo.LewisABx,
                             'lewis_acid_base/y': dlvo.LewisABy,
                             'lvdw/x': dlvo.LVDWx,
                             'lvdw/y': dlvo.LVDWy,
                             'attractive/x': dlvo.attractive_x,
                             'attractive/y': dlvo.attractive_y,
                             'edl_fine': fine_dlvo.EDLx,
                             'attractive_fine': fine_dlvo.attractive_x,
                             'distance_fine': fine_arr,
                             'distance_arr': arr}
            for key, value in dlvo_profiles.items():
                h.create_dataset('colloids/{}'.format(key), data=value)

            h.create_dataset('colloids/ux', data=ux)
            h.create_dataset('colloids/uy', data=uy)
            h.create_dataset('colloid_colloid/x', data=colloidcolloid.x)
            h.create_dataset('colloid_colloid/y', data=colloidcolloid.y)
            h.create_dataset('colloid_colloid/fine/x', data=cc_fine.x)
            h.create_dataset('colloid_colloid/fine/y', data=cc_fine.y)
            h.create_dataset('colloid_colloid/distance/x',
                             data=colloidcolloid.x_distance_array)
            h.create_dataset('colloid_colloid/distance/y',
                             data=colloidcolloid.y_distance_array)
            h.create_dataset('colloid_colloid/fine/distance/x',
                             data=cc_fine.x_distance_array)
            h.create_dataset('colloid_colloid/fine/distance/y',
                             data=cc_fine.y_distance_array)

            # todo: add the drag force arrays to the dataset.
            # todo: add method to convert dlvo forces to kT (first to J then 4.11*10^-21)
            # todo: check the DLVO adjustments that I included on Monday
            h.close()
