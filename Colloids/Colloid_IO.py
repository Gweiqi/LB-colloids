import sys

class Config:
    def __init__(self, fname):
        self.config = self._reader(fname)
        self._strtype = ('LBMODEL', 'ENDPOINT', 'PATHLINE', 'INTERVAL')
        self._inttype = ('NCOLS', 'ITERS', 'STORE_TIME', 'PRINT_TIME')
        self._floattype = ('LBRES', 'GRIDREF', 'AC', 'TIMESTEP', 'TEMPERATURE',
                           'RHO_WATER', 'RHO_COLLOID', 'VISCOSITY', 'I_INITIAL',
                           'I', 'EPSILON_R', 'SHEER_PLANE', 'LVDWST_WATER',
                           'LVDWST_COLLOID', 'LVDWST_SOLID', 'ZETA_COLLOID',
                           'ZETA_SOLID', 'PSI+_COLLOID', 'PSI+_WATER', 'PSI+_SOLID',
                           'PSI-_COLLOID', 'PSI-_WATER', 'PSI-_SOLID')
        self._booltype = ('PLOT', 'ADJUST_ZETA')
        self._dicttype = ('CONCENTRATION', 'VALENCE')
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
                                   'ENDPOINT', 'PATHLINE', 'INTERVAL')

    def _reader(self, fname):
        with open(fname, 'r') as f:
            config = [line.strip('\n').strip(' ') for line in f]
        return config

    def model_parameters(self):
        ModelDict = {}
        blockname = 'MODEL PARAMETERS'

        modelparams = self.get_block(blockname)
        for parameter in modelparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validmodelparams)
            pname = self.adjust_pname(pname)
            ModelDict[pname] = param
            
        return ModelDict

    def physical_parameters(self):
        PhysicsDict = {}
        blockname = 'PHYSICAL PARAMETERS'

        physicalparams = self.get_block(blockname)
        for parameter in physicalparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validphysicalparams)
            pname = self.adjust_pname(pname)
            PhysicsDict[pname] = param    
            
        return PhysicsDict

    def chemical_parameters(self):
        ChemicalDict = {}
        blockname = 'CHEMICAL PARAMETERS'

        chemicalparams = self.get_block(blockname)
        for parameter in chemicalparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validchemicalparams)
            pname = self.adjust_pname(pname)
            ChemicalDict[pname] = param
        return ChemicalDict

    def output_control(self):
        OutputDict = {}
        blockname = 'OUTPUT CONTROL'

        outputparams = self.get_block(blockname)
        for parameter in outputparams:
            pname, param = self.parametertype(parameter)
            self.check_if_valid(blockname, pname, self.validoutputparams)
            pname = self.adjust_pname(pname)
            OutputDict[pname] = param
            
        return OutputDict

    def get_block(self, blockname):
        try:
            idx0 = self.config.index('START %s' % blockname)
            idx1 = self.config.index('END %s' % blockname)
        except ValueError:
            print('Please define a %s block' % blockname)
            sys.exit(-1)
        return self.config[idx0+1:idx1]

    def parametertype(self, parameter):
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
                # maybe use a pass here and define defaults
                param = False

        elif pname in self._dicttype:
            keyvalue = [i for i in param.split(' ') if i != ' ']
            param = {}
            for i in range(0, len(keyvalue), 2):
                param[keyvalue[i]] = float(keyvalue[i+1])              
                
        else:
            print 'Parameter name %s is not valid' % pname
            sys.exit(-1)

        # pname = pname.lower()
        return pname, param

    def check_if_valid(self, blockname, pname, validparams):
        if pname in validparams:
            return
        else:
            print('%s is not a valid parameter for %s' % (pname, blockname))
            print(validparams)
            sys.exit(-1)

    def adjust_pname(self, pname):
        if pname in ('I_INITIAL', 'I', 'TEMPERATURE'):
            if pname == 'I_INITIAL':
                return 'I_initial'
            elif pname == 'I':
                return 'I'
            elif pname == 'TEMPERATURE':
                return 'T'
            else:
                pass
        else:
            return pname.lower()
        
test = Config('Synthetic.config')
print test.model_parameters()
print test.physical_parameters()
print test.output_control()
print test.chemical_parameters()

