"""
lbIO is the lattice Boltzmann input control module. This module contains a Config class that
reads in the lattice Boltzmann configuration file and parses the parameters
contained within that file. HDF5 write options will eventually be moved from other
portions of the lattice Boltzmann modeling system to LBIO.py

Example import of a lattice Boltzmann configuration file is shown below

>>> from lb_colloids import lbIO
>>> config_file = "LB_model.config"
>>> config = lbIO.Config(config_file)
>>> model_dict = config.model_parameters()
>>> image_dict = config.image_parameters()
>>> permeability_dict = config.permeability_parameters()
>>> output_dict = config.output_parameters()
"""

import sys


class Config:
    """
    Class to handle the input output control of the lattice botzmann modeling
    system. Can be added to, and uses consistancy checks to ensure correct
    types and required values are present for LBModel to run.

    Parameters:
    ----------
    :param str fname: configuration file name
    """
    def __init__(self, fname):

        self.config = self._reader(fname)
        self.__strtype = ('LBMODEL', 'KERNEL', 'IMAGE_SAVE_FOLDER',
                         'IMAGE', 'IMAGE_SAVE_NAME')
        self.__inttype = ('SOLID', 'VOID', 'NITERS', 'VERBOSE',
                         'IMAGE_SAVE_INTERVAL', 'BOUNDARY')
        self.__listtype = ('SOLID', 'VOID')
        self.__floattype = ('LBRES', 'RHO', 'TAU', 'VMIN',
                           'VMAX', 'GRAVITY', 'PHYSICAL_VISCOSITY',
                            'PHYSICAL_RHO')
        self.__booltype = ('PLOT',)
        self.__dicttype = ()
        self.__required = ('LBMODEL', 'LBRES')
        self.validmodelparams = ('LBMODEL', 'LBRES', 'KERNEL',
                                 'PHYSICAL_VISCOSITY', 'PHYSICAL_RHO')
        self.validimageparams = ('SOLID', 'VOID', 'IMAGE', 'BOUNDARY', 'PLOT')
        self.validpermeabilityparmas = ('RHO', 'TAU', 'NITERS', 'GRAVITY')
        self.validoutputparams = ('IMAGE_SAVE_NAME', 'IMAGE_SAVE_FOLDER', 'IMAGE_SAVE_INTERVAL',
                                  'VMIN', 'VMAX', 'VERBOSE')

    def _reader(self, fname):
        with open(fname, 'r') as f:
            config = [line.strip('\n').strip(' ') for line in f]
        return config

    def model_parameters(self):
        """
        reads the MODEL PARAMETERS block of the configuration file and creates
        the ModelDict which contains the required parameters to run LB-Colloid

        Returns:
        -------
        :return: (dict) Dictionary containing model parameter block information
        """
        ModelDict = {}
        blockname = 'MODEL PARAMETERS'

        modelparams = self.get_block(blockname)
        for parameter in modelparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validmodelparams)
            # pname = self.adjust_pname(pname)
            ModelDict[pname] = param
        
        self.check_model_parameters(ModelDict)
        return ModelDict

    def image_parameters(self):
        """
        reads the IMAGE PARAMETERS block of the configuration file and creates
        the ModelDict which contains the required parameters to run LB-Colloid

        Returns:
        -------
        :return: (dict) Dictionary containing image paramaters for simulation
        """
        ImageDict = {}
        blockname = 'IMAGE PARAMETERS'

        imageparams = self.get_block(blockname)
        for parameter in imageparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validimageparams)
            ImageDict[pname] = param

        return ImageDict        
    
    def permeability_parameters(self):
        """
        reads the PERMEABILITY PARAMETERS block of the configuration file and creates
        the ModelDict which contains the required parameters to run LB-Colloid

        Returns:
        -------
        :return: (dict) Dictionary containing parameters from the Permeability input block
        """
        PermeabilityDict = {}
        blockname = 'PERMEABILITY PARAMETERS'

        permeabilityparams = self.get_block(blockname)
        for parameter in permeabilityparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validpermeabilityparmas)
            # pname = self.adjust_pname(pname)
            PermeabilityDict[pname] = param
        
        return PermeabilityDict

    def output_parameters(self):
        """
        reads the OUTPUT CONTROL block of the configuration file and creates
        the ModelDict which contains the required parameters to run LB-Colloid

        Returns:
        -------
        :return: (dict) Dictionary containing output control parameters
        """
        OutputDict = {}
        blockname = 'OUTPUT CONTROL'

        outputparams = self.get_block(blockname)
        for parameter in outputparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validoutputparams)
            # pname = self.adjust_pname(pname)
            OutputDict[pname] = param
        
        return OutputDict

    def get_block(self, blockname):
        """
        Method to locate the beginning and end of each input block

        Parameters:
        ----------
        :param str blockname: blockname of a parameters block in the config file.

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
        :param str parameter: string of "<pname>: <param>"

        returns:
        -------
        :return: (tuple) parameter name, parameter associated with pname
        """
        # split the paramter block line, and strip whitespace
        pname, param = parameter.split(':')
        pname = pname.strip(' ')
        param = param.strip(' ')

        if pname in self.__listtype:
            param = self.list_param(pname, param)
        
        elif pname in self.__strtype:
            param = str(param)
            
        elif pname in self.__inttype:
            param = int(param)

        elif pname in self.__floattype:
            param = float(param)

        elif pname in self.__booltype:
            if param.upper() == 'TRUE':
                param = True
            else:
                param = False

        elif pname in self.__dicttype:
            keyvalue = [i for i in param.split(' ') if i != ' ']
            param = {}
            for i in range(0, len(keyvalue), 2):
                param[keyvalue[i]] = float(keyvalue[i+1])              
                
        else:
            print 'Parameter name %s is not valid' % pname
            sys.exit(-1)

        return pname, param

    def list_param(self, pname, param):
        params = param.split(' ')

        if pname in self.__strtype:
            param = [str(i) for i in params if i != ' ']
            
        elif pname in self.__inttype:
            param = [int(i) for i in params if i != ' ']

        elif pname in self.__floattype:
            param = [float(i) for i in params if i != ' ']

        return param

    def check_if_valid(self, blockname, pname, validparams):
        """
        Method checks if a specific parameter is valid for the block it was supplied
        in the configuration file.

        Parameters:
        ----------
        :param str blockname: the name of the configuration block
        :param str pname: parameter name
        :param tupele validparams: tuple of valid parameters for the specific configuration block
        """
        if pname in validparams:
            return
        else:
            print('%s is not a valid parameter for %s' % (pname, blockname))
            print(validparams)
            sys.exit(-1)

    def check_model_parameters(self, ModelDict):
        """
        Check for required parameters. If not present inform the user which parameter(s) are needed

        Parameters:
        ----------
        :param dict ModelDict: Model Input block dictionary
        """
        for key in self.__required:
            if key not in ModelDict:
                print('%s is a required MODEL PARAMETER' % key)
                sys.exit(-1)
            else:
                pass
