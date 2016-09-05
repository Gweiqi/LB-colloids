import sys

class Config:
    def __init__(self, fname):
        '''
        Class to handle configuration files for LB-Colloids. Checks in place
        to look for consistancy in data type and config variable per data block.

        Input:
        ------
        fname: (string) configuration file name ie. Model.config

        Returns:
        --------
        model_parameters: (dict) Dictionary of necessary model parameters
        physical_parameters: (dict) Dictionary of optional physical parameters
        chemical_parmaeters: (dict) Dictionary of chemical parameter options
        output_control: (dict) Dictionary of output control options
        
        '''
        self.config = self._reader(fname)
        self._strtype = ('LBMODEL', 'ENDPOINT', 'PATHLINE', 'TIMESERIES')
        self._inttype = ('NCOLS', 'ITERS', 'STORE_TIME', 'PRINT_TIME')
        self._floattype = ('LBRES', 'GRIDREF', 'AC', 'TIMESTEP', 'TEMPERATURE',
                           'RHO_WATER', 'RHO_COLLOID', 'VISCOSITY', 'I_INITIAL',
                           'I', 'EPSILON_R', 'SHEER_PLANE', 'LVDWST_WATER',
                           'LVDWST_COLLOID', 'LVDWST_SOLID', 'ZETA_COLLOID',
                           'ZETA_SOLID', 'PSI+_COLLOID', 'PSI+_WATER', 'PSI+_SOLID',
                           'PSI-_COLLOID', 'PSI-_WATER', 'PSI-_SOLID')
        self._booltype = ('PLOT', 'ADJUST_ZETA', 'OVERWRITE')
        self._dicttype = ('CONCENTRATION', 'VALENCE')
        self._required = ('LBMODEL', 'NCOLS', 'ITERS', 'LBRES', 'GRIDREF',
                          'TS')
        self.validmodelparams = ('LBMODEL', 'NCOLS', 'ITERS', 'LBRES',
                                 'GRIDREF', 'AC', 'TIMESTEP', 'TEMPERATURE')
        self.validphysicalparams = ('RHO_WATER', 'RHO_COLLOID', 'VISCOSITY')
        self.validchemicalparams = ('CONCENTRATION', 'ADJUST_ZETA', 'I_INITIAL',
                                    'I', 'EPSILON_R', 'VALENCE', 'SHEER_PLANE',
                                    'LVDWST_WATER', 'LVDWST_COLLOID', 'LVDWST_SOLID',
                                    'ZETA_COLLOID', 'ZETA_SOLID', 'PSI+_COLLOID',
                                    'PSI+_WATER', 'PSI+_SOLID', 'PSI-_COLLOID',
                                    'PSI-_WATER', 'PSI-_SOLID')
        self.validoutputparams = ('STORE_TIME', 'PRINT_TIME', 'PLOT',
                                   'ENDPOINT', 'PATHLINE', 'TIMESERIES',
                                  'OVERWRITE', 'PLOT')

    def _reader(self, fname):
        '''
        reads in the input config file
        '''
        with open(fname, 'r') as f:
            config = [line.strip('\n').strip(' ') for line in f]
        return config

    def model_parameters(self):
        '''
        reads the MODEL PARAMETERS block of the configuration file and creates
        the ModelDict which contains the required parameters to run LB-Colloid
        '''
        ModelDict = {}
        blockname = 'MODEL PARAMETERS'

        modelparams = self.get_block(blockname)
        for parameter in modelparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validmodelparams)
            pname = self.adjust_pname(pname)
            ModelDict[pname] = param
        self.check_model_parameters(ModelDict)
        return ModelDict

    def physical_parameters(self):
        '''
        reads the PHYSICAL PARAMETERS block of the configuration file and creates
        the PhysicsDict that passes optional parameters to the physical force calculations
        in the Colloid_Math.py module
        '''
        PhysicsDict = {}
        blockname = 'PHYSICAL PARAMETERS'

        physicalparams = self.get_block(blockname)
        for parameter in physicalparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validphysicalparams)
            pname = self.adjust_pname(pname)
            PhysicsDict[pname] = param

        PhysicsDict = self.add_universal_parameters(PhysicsDict)
            
        return PhysicsDict

    def chemical_parameters(self):
        '''
        reads the CHEMICAL PARAMETERS block of the configuration file and creates
        the ChemicalDict that passes optional parameters to the chemical force calculations
        in the Colloid_Math.py module
        '''
        ChemicalDict = {}
        blockname = 'CHEMICAL PARAMETERS'

        chemicalparams = self.get_block(blockname)
        for parameter in chemicalparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validchemicalparams)
            pname = self.adjust_pname(pname)
            ChemicalDict[pname] = param

        ChemicalDict = self.add_universal_parameters(ChemicalDict)
        return ChemicalDict

    def output_control(self):
        '''
        reads the OUTPUT CONTROL block of the configuration file and creates
        the OutputDict that passes optional parameters to control model output
        '''
        OutputDict = {}
        blockname = 'OUTPUT CONTROL'

        outputparams = self.get_block(blockname)
        for parameter in outputparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validoutputparams)
            pname = self.adjust_pname(pname)
            OutputDict[pname] = param
        
        OutputDict = self.add_universal_parameters(OutputDict)
        return OutputDict

    def get_block(self, blockname):
        '''
        Input:
        ------
        blockname: (str) blockname of a parameters block in the config file.

        Returns:
        --------
        list: (list) returns a list of parameters contained within the block
        '''
        try:
            idx0 = self.config.index('START %s' % blockname)
            idx1 = self.config.index('END %s' % blockname)
        except ValueError:
            print('Please define a %s block' % blockname)
            sys.exit(-1)
        return self.config[idx0+1:idx1]

    def parametertype(self, parameter):
        '''
        Method takes a parameter string, splits it, and sets the parameter type

        Input:
        ------
        parameter: (string) string of "<pname>: <param>"

        returns:
        --------
        pname: (string) parameter name
        param: (type dependent) parameter associated with pname
        '''
        # split the paramter block line, and strip whitespace
        pname, param = parameter.split(':')
        pname = pname.strip(' ')
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
        '''
        Method checks if a specific parameter is valid for the block it was supplied
        in the configuration file.

        Input:
        ------
        blockname: (string) the name of the configuration block
        pname: (string) parameter name
        validparams: (tuple, string) tuple of valid parameters for the specific configuration block
        '''
        if pname in validparams:
            return
        else:
            print('%s is not a valid parameter for %s' % (pname, blockname))
            print(validparams)
            sys.exit(-1)

    def adjust_pname(self, pname):
        '''
        Adjusts parameter name from configuration file name to LB-Colloids name for a limited number of
        parameters. Sets all other parameters as lowercase to follow PEP-8
        '''
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
        '''
        Add common model parameters to other dictionaries if present. Necessary for parameterization by
        kwargs of physics and chemistry.
        '''
        modelparams = self.model_parameters()
        for key in modelparams:
            if key in ('ac', 'T', 'visocity', 'ts', 'lbres', 'gridref'):
                Dict[key] = modelparams[key]
        return Dict

    def check_model_parameters(self, ModelDict):
        '''
        Check for required parameters. If not present inform the user which parameter(s) are needed
        '''
        for key in self._required:
            if key.lower() not in ModelDict:
                if key.lower() == 'ts':
                    key = 'TIMESTEP'
                else:
                    pass
                print('%s is a required MODEL PARAMETER' % key)
                sys.exit(-1)
            else:
                pass


class Output:
    def __init__(self, fi, **kwargs):

        defaults = {'overwrite': True}

        for key in kwargs:
            defaults[key] = kwargs[key]

        self.filename = fi
        self.header = header = '{:>8}{:>9}{:>11}{:>14}{:>14}{:>14}{:>14}{:>14}\n'.format(
                               'colloid', 'ts', 'totim', 'x-position', 'y-position',
                               'resolution', 'x-model', 'y-model')
        self.resolution = defaults['lbres']/defaults['gridref']
        
        #checks if file exists, and overwrites if overwrite is True!
        if defaults['overwrite'] is True:
            self._writer(fi, header, wtype='w')
        
    def _writer(self, fi, data, wtype='a'):
        with open(fi, wtype) as f:
            for line in data: # format as 1d list that is str fomated with \n
                f.write(line)


    def write_output(self, timer, colloids, pathline=True):
        """
        set up and write full output to an ascii file

        [format is [time, totim, col#, xpos, ypos, resolution, modelx, modely] Add flag to this?
        """
        # todo: add support for grid location, model flag.
        time = timer.timer
        totim = timer.totim
        output = []
        
        if pathline is not True:
            for colloid_number in range(len(colloids)):
                output_string = '{:>8d}    {:5d}    {:07.5f}    {:09.8f}    {:09.8f}    {:09.8f}    {:09.8f}    {:09.8f}\n'.format(
                                colloid_number, time[-1] , totim[-1],
                                colloids[colloid_number].xposition[-1],
                                colloids[colloid_number].yposition[-1],
                                self.resolution,
                                colloids[colloid_number].xposition[-1]/self.resolution,
                                colloids[colloid_number].yposition[-1]/self.resolution)
                output.append(output_string)
        
        else:    
            for idx in range(len(time)-1): # so we don't duplicate
                for colloid_number in range(len(colloids)):
                    output_string = '{:>8d}    {:5d}    {:07.5f}    {:09.8f}    {:09.8f}    {:09.8f}    {:09.8f}    {:09.8f}\n'.format(
                                    colloid_number, time[idx] , totim[idx],
                                    colloids[colloid_number].xposition[idx],
                                    colloids[colloid_number].yposition[idx],
                                    self.resolution,
                                    colloids[colloid_number].xposition[idx]/self.resolution,
                                    colloids[colloid_number].yposition[idx]/self.resolution)
                    output.append(output_string)

        self._writer(self.filename, output)
        

            
        
