from lb_colloids import NamFile
from lb_colloids import lbIO
from lb_colloids import LBImage
from lb_colloids import LB2DModel
from lb_colloids import cIO
from lb_colloids import ColloidModel

import matplotlib.pyplot as plt
import os


def addIO(defaults, config):
    """
    Simple method to update defaults from IO
    """
    for key in config:
        defaults[key] = config[key]
    return defaults


nam_files = [fi for fi in os.listdir(os.getcwd()) if fi.endswith(".nam")]

# trap for common nam file errors
if len(nam_files) > 1:
    raise AssertionError("Only one nam file can be present in a directory")

elif len(nam_files) == 0:
    raise AssertionError("A <>.nam file must be present to run model")

else:
    pass

nam_file = nam_files[0]

configs = NamFile(nam_file)

# set up lattice boltzman model parameters
if configs.lb_config is not None:
    lb_params = lbIO.Config(configs.lb_config)

    ImageDict = {'BOUNDARY': 10}
    ModelDict = {'KERNAL': 'fortan'}
    PermeabilityDict = {'RHO': 1.0, 'TAU': 1.0, 'GRAVITY': 0.001, 'NITERS': 1,
                        'PHYSICAL_VISCOSITY': 8.9e-4, 'PHYSICAL_RHO': 997.}
    OutputDict = {'VERBOSE': 100,
                  'IMAGE_SAVE_INTERVAL': None,
                  'IMAGE_SAVE_FOLDER': os.path.expanduser('~/Desktop/LBimages'),
                  'IMAGE_SAVE_NAME': 'LB',
                  'VMAX': 0.,
                  'VMIN': -0.010}

    ImageDict = addIO(ImageDict, lb_params.image_parameters())
    ModelDict = addIO(ModelDict, lb_params.model_parameters())
    PermeabilityDict = addIO(PermeabilityDict, lb_params.permeability_parameters())
    OutputDict = addIO(OutputDict, lb_params.output_parameters())

    plot_binary = ImageDict.pop('PLOT', False)
    infile = ImageDict['IMAGE']
    fluid = ImageDict['VOID']
    solid = ImageDict['SOLID']
    nlayers = ImageDict['BOUNDARY']


    # run LB_2Dimage functions to binarize the domain
    print('[Reading image]')
    img = LBImage.Images(infile)

    print('[Setting boundary condition]')
    binary = LBImage.BoundaryCondition(img.arr, fluid, solid, nlayers)
    print('[Porosity: %.4f]' % binary.porosity)

    if plot_binary:
        plt.imshow(binary.binarized, interpolation = 'nearest')
        plt.show()

    out = LBImage.HDF5_write(binary.binarized, binary.porosity, nlayers, ModelDict['LBMODEL'])
    print('[Binarization Finished]')

    # run LB_2Dpermeability for computational fluid dynamics
    print('[Building LB Model]')
    lbmodel = LB2DModel(binary.binarized, kernel=ModelDict['KERNEL'])

    lbmodel.niters = PermeabilityDict['NITERS']
    lbmodel.tau = PermeabilityDict['TAU']
    lbmodel.gravity = PermeabilityDict['GRAVITY']
    lbmodel.rho = PermeabilityDict['RHO']
    lbmodel.physical_viscosity = PermeabilityDict['PHYSICAL_VISCOSITY']
    lbmodel.physical_rho = PermeabilityDict['PHYSICAL_RHO']

    print('[Running LB permeability model]')
    lbmodel.run(output=ModelDict["LBMODEL"], verbose=OutputDict['VERBOSE'],
                image_int=OutputDict['IMAGE_SAVE_INTERVAL'], image_folder=OutputDict["IMAGE_SAVE_FOLDER"],
                image_name=OutputDict["IMAGE_SAVE_NAME"], vmax=OutputDict['VMAX'], vmin=OutputDict['VMIN'])
    print('[Normal Termination of Simulation: LB Model Complete]')

# set up and run colloids model!
if configs.colloid_config:
    print('[Loading LB Colloid Config file(s)]')
    configs = [cIO.Config(config) for config in configs.colloid_config]
    print('[Config load Finished]')

    print('[Running LB Colloid Model]')
    ColloidModel.run(configs)
    print('[Normal Termination of Simulation: Colloid Model Complete]')