#!/bin/bash

#python LB_2Dimage.py -i 2D_256_255.png -c f0,s255 -o 255f-tstr.hdf5 -b 10
#python LB_2Dimage.py -i 2D_256_255.png -c f0,s255 -o 255py-tstr.hdf5 -b 10

#python LB_2Dpermeability.py -i 255f-tstr.hdf5 -n 3001  -q 100 -p 50 -f ftestr -u pu
python LB_2Dpermeability.py -i 255py-tstr.hdf5 -n 3001  -q 100 -p 50 -f pytestr -u pu -k p

