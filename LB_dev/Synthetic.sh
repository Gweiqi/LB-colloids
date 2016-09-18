#!/bin/bash

python LB_2Dimage.py -i synthetic.png -c f255,s0 -o tstr.hdf5 -b 10

python LB_2Dpermeability.py -i tstr.hdf5 -n 3001  -q 100 -p 50 -f testr -u pu

