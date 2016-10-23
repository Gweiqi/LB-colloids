import sys

class Config:
    def __init__(self, fname):
        """
        Class to handle the input output control of the lattice botzmann modeling
        system. Can be added to, and uses consistancy checks to ensure correct
        types and required values are present for LBModel to run.
        """
        self.config = self._reader(fname)
        self._strtype = ('LBMODEL', 'KERNAL', 'IMAGE_SAVE_FOLDER',
                         'IMAGE')
        self._inttype = ('SOLID', 'VOID', 'NITERS', 'VERBOSE',
                         'IMAGE_SAVE_INTERVAL', 'BOUNDARY')
        self._floattype = ('LBRES', 'RHOT', 'RHOB', 'TAU', 'VMIN',
                           'VMAX', 'GRAVITY')
        self._booltype = ('PLOT_Y_VELOCITY', 'SAVE_IMAGE')
        self._dicttype = ()
        self._required = ('LBMODEL', 'LBRES')
        self.validmodelparams = ('LBMODEL', 'LBRES', 'KERNAL')
        self.validimageparams = ('SOLID', 'VOID', 'IMAGE', 'BOUNDARY')
        self.validpermeabilityparmas = ('RHOT', 'RHOB', 'TAU', 'NITERS', 'GRAVITY')
        self.validoutputparams = ('SAVE_IMAGE', 'IMAGE_SAVE_FOLDER', 'IMAGE_SAVE_INTERVAL',
                                  'VMIN', 'VMAX', 'VERBOSE', 'PLOT_Y_VELOCITY')

    def _reader(self, fname):
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
            # pname = self.adjust_pname(pname)
            ModelDict[pname] = param
        
        self.check_model_parameters(ModelDict)
        return ModelDict

    def image_parameters(self):
        '''
        reads the IMAGE PARAMETERS block of the configuration file and creates
        the ModelDict which contains the required parameters to run LB-Colloid
        '''
        ImageDict = {}
        blockname = 'IMAGE PARAMETERS'

        imageparams = self.get_block(blockname)
        for parameter in imageparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validimageparams)
            ImageDict[pname] = param

        return ImageDict        
    
    def permeability_parameters(self):
        '''
        reads the PERMEABILITY PARAMETERS block of the configuration file and creates
        the ModelDict which contains the required parameters to run LB-Colloid
        '''
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
        '''
        reads the OUTPUT CONTROL block of the configuration file and creates
        the ModelDict which contains the required parameters to run LB-Colloid
        '''
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

    def check_model_parameters(self, ModelDict):
        '''
        Check for required parameters. If not present inform the user which parameter(s) are needed
        '''
        for key in self._required:
            if key not in ModelDict:
                print('%s is a required MODEL PARAMETER' % key)
                sys.exit(-1)
            else:
                pass
