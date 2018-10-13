.. LB-Colloids documentation master file, created by
   sphinx-quickstart on Wed Oct  4 10:53:29 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the LB-Colloids documentation!
=========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:


Introduction to LB-Colloids
===========================

Lattice Boltzmann colloids is an object oriented python package that builds and simulates two dimensional, nine fluid node (D2Q9) computational fluid dynamic models. Lattice Boltzmann is able to represent fluid domains defined by complex geometries ehich are observed in natural porous media. Simulated fluid flow is slightly compressible and has been shown to return an approximation of the Navier-Stokes equation (*Benzi et al 1992*). The lattice Boltzmann code used for this software is based on the single relaxation time Bhantager Gross Krook formulation of lattice Boltzmann. Simple bounceback boundary conditions are implemented at fluid boundaries to preserve mass balance. As a result of the bounceback boundary conditions no-slip conditions develop at pore-solid interfaces. 

Binarized imagery is used to represent the fluid flow domain in this implementation of lattice Boltzmann. Distinct advantages of using binarized imagery are: fluid flow in natural porous media can be simulated through the use of segmented thin sections, and synthetic porous media can easily be generated as an image and used for fluid dynamic simulations. 

Colloid particle tracking simulations can be parameterized and simulated for steady state lattice Boltzmann computational fluid dynamic models using this software. Colloid transport through the soil environment is of great interest and importance to soil development processes through the translocation of clays, contaminant transport (*Saiers 1996, Jaisi et al 2008*), filtration and transport of bio-colloids (*Harter 2000, Redman 2004, Foppen et al 2005*), and soil nutrient dynamics. Microscale colloid particle tracking simulations allow the user to gain additional insight into physio-chemical processes influencing colloid transport. Coupling macroscopic breakthrough curve information with the microscopic distribution of colloids during a simulation provides insight into attachment and colloid retention mechanisms for simulations. These insights have the potential to illuminate retention and transport mechanisms for column and field scale studies of colloid and contaminant transport. 

Parameterization of an LB-Colloids model is possible in two ways. Formatted text documents can be supplied to the LB-Colloids simulation software and results will be stored in an hdf5 file and a human readable text document. Advanced users may prefer to parameterize simulations through a python interface. This gives the advanced user the ability to easily calibrate models by adjusting parameters within python code blocks, run multiple models, and perform sensitivity analysis. Currently output readers and built in analysis methods are only available to the advanced user. API documentation serves as a guide to the advanced user on importing and using these python objects. 

This document serves as the user guide to the LB-Colloids python package. For more details pertaining to the development of this package see *Larsen and Schaap TBD*. The user guide follows the structure:

Chapter 2: Table of mathematical symbols

Chapter 3: Installation of LB-Colloids

Chapter 4: Formatted text input

Chapter 5: Parameterization of LB-Colloids

Chapter 6: Lattice Boltzmann API

Chapter 7: Colloid Simulation API

Chapter 8: Penetrable-Sphere (PSHPERE) API

Table of Mathematical Symbols
=============================
Lattice Boltzmann
-----------------

:math:`\pmb{e_{i}}`: lattice Boltzmann eigenvector array

:math:`f_{i}` :  distribution function

:math:`f_{eq}` : equilibrium distribution function

:math:`\rho` :  macroscopic fluid density

:math:`\tau` :  relaxation time 

:math:`\pmb{u}` :  macroscopic fluid velocity

:math:`v` : fluid viscosity

:math:`w_{i}` : fluid link weights


Colloid Simulation
------------------

:math:`A_{H}` : Hamaker constant

:math:`a_{c}` : colloid radius

:math:`D_{0}` : diffusivity

:math:`\epsilon_{0}` : dielectric permativity of a vacuum

:math:`\epsilon_{r}` : dielectric permativity of water

:math:`e` : electron charge

:math:`F^{b}` : bouyancy force

:math:`F^{B}` : brownian force

:math:`F^{D}` : drag force

:math:`F^{G}` : gravity force

:math:`f_{1}` : hydrodynamic correction factor

:math:`f_{2}` : hydrodynamic correction factor

:math:`f_{3}` : hydrodynamic correction factor

:math:`f_{4}` : hydrodynamic correction factor

:math:`G(0, 1)` : random gaussian distribution 

:math:`g` : acceleration due to gravity

:math:`h` : gap distance

:math:`\bar{h}` : non-dimensional gap distance

:math:`I^{*}` : Two times fluid ionic strength.

:math:`k` : Boltzmann constant

:math:`\kappa` : Debye length

:math:`M` : molarity

:math:`N_{A}` : Avagadro's number

:math:`\Phi^{A}` : Attractive interaction energy

:math:`\Phi^{EDL}` : Electric double layer interaction energy

:math:`\psi_{c}` : colloid potential

:math:`\psi_{s}` : surface potential

:math:`\rho_{c}` : colloid density

:math:`\rho_{w}` : fluid density

:math:`u` : fluid velocity

:math:`\mu` : fluid viscosity

:math:`V` : colloid velocity

:math:`T` : Fluid temperature

:math:`Z` : cation charge

:math:`\xi` : :math:`6\pi \mu a_{c}`

Installation of LB-Colloids
===========================
LB-Colloids is currently avaialable for python 2.7.6 - 2.7.12 installations. Support is currently not available for python 3, however most source code can be converted to python 3 with little trouble. Future support is planned for python 3.5

Recommended installation of LB-Colloids is performed using the python tool pip. If pip is not included with your python distribution an installation script can be found at https://bootstrap.pypa.io/get-pip.py. It is recommended that the user adds pip as a path variable for easy installation of python packages. 

An installation of gfortran must be present on the user's computer to properly install LB-Colloid's. Mathematic modules are compiled locally in FORTRAN upon install and called by python for computational efficiency. Gfortran should be added as a path variable to your computer.

LB-Colloids can be downloade from *pypi* and/ort *tbd*. Move the LB-Colloids package to your prefered source code location. 

Open a terminal and navigate to the base directory of the LB-Colloids source code. You should see setup.py in this directory. Install LB-Colloids using pip.

>>> cd Desktop/LB-Colloids
>>> pip install .
>>> # or for the developer use
>>> pip install -e .

Congratulations! LB-Colloids is now installed on your machine.

.. include:: ./ruid.rst

Lattice Boltzmann
=================
Lattice Boltzmann boundary conditions
-------------------------------------

.. include:: ./lbbc.rst

API documentation
*****************
.. py:module:: lb_colloids
.. autoclass:: LBImage
.. automodule:: lb_colloids.LB.LB_2Dimage
	:members:

Lattice Boltzmann Permeability
------------------------------

.. include:: ./lbperm.rst

API documentation
*****************
.. automodule:: lb_colloids.LB.LB_2Dpermeability
	:members:
.. automodule:: lb_colloids
.. autoclass:: LB2DModel

Lattice Boltzmann Input Output
------------------------------

.. include:: ./lbio.rst

API documentation
*****************
.. py:module:: lb_colloids
.. autoclass:: lbIO
.. automodule:: lb_colloids.LB.LBIO
	:members:

.. autoclass:: lb_colloids.LB.LB_pretty
.. automodule:: lb_colloids.LB.LB_pretty
	:members:

Colloid Simulation
==================
LB Colloid simulations
----------------------

.. include:: ./csim.rst

API documentation
*****************
.. py:module:: lb_colloids
.. autoclass:: ColloidModel
.. automodule:: lb_colloids.Colloids.LB_Colloid
	:members:

LB Colloid mathematics
----------------------

.. include:: ./cmath.rst

API documentation
*****************
.. py:module:: lb_colloids
.. autoclass:: ColloidMath
.. automodule:: lb_colloids.Colloids.Colloid_Math
	:members:

LB Colloid Input Output
-----------------------

.. include:: ./cio.rst

|todo: find a place to include run_model.py

API documentation
*****************
.. py:module:: lb_colloids
.. autoclass:: cIO
.. automodule:: lb_colloids.Colloids.Colloid_IO
	:members:
.. py:module:: lb_colloids
.. autoclass:: ColloidOutput
.. automodule:: lb_colloids.Colloids.Colloid_output
	:members:
.. automodule:: lb_colloids.nam_file
	:members:
	
LB Colloid Setup (background classes)
-------------------------------------

.. include:: cset.rst

API documentation
*****************
.. automodule:: lb_colloids.Colloids.Colloid_Setup
	:members:

Penetrable Sphere (PSHPERE)
===========================

.. include:: ./pshpere.rst

API documentation
*****************
.. py:module:: lb_colloids
.. autoclass:: PSphere
.. automodule:: lb_colloids.utilities.psphere
	:members:
