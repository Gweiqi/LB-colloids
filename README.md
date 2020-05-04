[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3245555.svg)](https://doi.org/10.5281/zenodo.3245555)

# LB-colloids
Lattice Boltzmann Colloids version 0.2 candidate

## LB Colloids is a package to simulate fluid flow in porous media using D2Q9 lattice Boltzmann and perform colloid particle tracking simulations

Installation of LB-Colloids is relatively easy and can be accomplished using the 
python installer pip. Download a zip file of the repository, unzip, and open a 
terminal window in that directory.

Installation can be accomplished using the command

```
python -m pip install .
```

or

```
python -m pip install . -e
```

for an editable version of the installation.

If there are problems with installation or importing LB-Colloids please raise 
an issue or send me an email. It may be a FORTRAN compiler issue that we can 
work out together. Windows support is in progress and will be included in
the 0.2 release.

A user document/package API for the LB-Colloids package is located in the 
Docs/build/latex directory as a PDF file.
[User Guide](https://github.com/jdlarsen-UA/LB-colloids/blob/develop/Docs/build/latex/LB-Colloids.pdf).  

[Examples](https://github.com/jdlarsen-UA/LB-colloids/tree/develop/examples/Notebooks) 
of how to use LB-Colloids using python scripting are located in the 
examples/Notebooks directory. These are jupyter notebooks that can be run 
interactively on the user's machine if Jupyter is installed.

Additional scripts and code samples can be found in the data directory.

## Required software
To install LB-Colloids the following software is required:  
A fortran compiler (ex. mingw, GFORTRAN)

On Windows:  
Microsoft visual C++ for python

Python > 2.7
   * numpy
   * matplotlib
   * h5py
   * pandas
   * scipy
   * jupyter
   
## How to reference this work

Larsen, Joshua. Ph.D., Pore Scale Computational Fluid Dynamic Modeling: 
Approaches for Permeability Modeling and Particle Tracking Using Lattice 
Boltzmann Methods, The University of Arizona 2018, 230 pages; 10978423

Open Access: https://pqdtopen.proquest.com/doc/2139258833.html?FMT=ABS

Software DOI: 10.5281/zenodo.3245555

## Authors
Joshua D. Larsen 