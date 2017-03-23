import LB_dev.LB_2Dimage as LBimage
import LB_dev.LB_2Dpermeability as LB2D
import matplotlib.pyplot as plt
import Colloids.Colloid_IO as cIO
import Colloids.LB_Colloid as Colloid
import os

path =  os.path.dirname(os.path.realpath(__file__))
lbname = os.path.join(path, "LBModel.hdf5")

imgname = os.path.join(path, "Synth100_2.png")
"""

nlayers = 5
solidvx = 0
fluidvx = 253


img = LBimage.Images(imgname)
binary = LBimage.BoundaryCondition(img.arr, fluidvx, solidvx, nlayers)
plt.imshow(binary.binarized, interpolation='nearest')
plt.show()

lbmodel = LB2D.LB2DModel(binary.binarized)
lbmodel.niters = 1000
lbmodel.run(output=lbname, verbose=10, image_int=10, image_folder='test')
"""
# now we work with LB_colloids.py
# need to think about namespace now.
# now we can use this to do other analysis such as sensitivity analysis

io = cIO.ColloidsConfig()
print io.valid_model_parameters
io['lbmodel'] = lbname
io['ncols'] = 500
io['iters'] = 3000
io['lbres'] = 1e-6
io['gridref'] = 10
io['ac'] = 1e-6
io['timestep'] = 1e-5
io['temperature'] = 298.
io['multiple_config'] = False

io['i'] = 1.0

io['print_time'] = 100
io['plot'] = True
io['endpoint'] = os.path.join(path, 'Synth100_3.endpoint')
io['store_time'] = 100

print io.model_parameters
print io.chemical_parameters
print io.physical_parameters
print io.output_control_parameters

# config = os.path.join(path, 'Synth100_3.config')

config = cIO.Config(io.config)
Colloid.run(config)
