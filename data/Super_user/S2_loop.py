from lb_colloids import LBImage
from lb_colloids import LB2DModel
import matplotlib.pyplot as plt
from lb_colloids import cIO
from lb_colloids import ColloidModel
import os

path =  os.path.dirname(os.path.realpath(__file__))
imgname = os.path.join(path, "Synth100_2.png")


nlayers = 5
solidvx = 0
fluidvx = 253

# NaCl at a range of 1-10 ds/m
ionics = [i * 680 * 0.001 /(35.5 + 23.) for i in range(1, 11)]

for ions in ionics:

    lbname = os.path.join(path, "s2_i{:.3}.hdf5".format(ions))

    img = LBImage.Images(imgname)
    binary = LBImage.BoundaryCondition(img.arr, fluidvx, solidvx, nlayers)
    lbmodel = LB2DModel(binary.binarized)
    lbmodel.niters = 1000
    lbmodel.run(output=lbname, verbose=10, image_int=10, image_folder='test')

    io = cIO.ColloidsConfig()
    io['lbmodel'] = lbname
    io['ncols'] = 50
    io['iters'] = 100000
    io['lbres'] = 1e-6
    io['gridref'] = 10
    io['ac'] = 1e-6
    io['timestep'] = 5e-7
    io['temperature'] = 298.
    io['continuous'] = 10000

    io['i'] = ions

    io['print_time'] = 10000
    io['plot'] = True
    io['endpoint'] = os.path.join(path, "S2_i{:.3}.endpoint".format(ions))
    io['store_time'] = 100

    print io.model_parameters
    print io.chemical_parameters
    print io.physical_parameters
    print io.output_control_parameters

    config = cIO.Config(io.config)
    ColloidModel.run(config)

