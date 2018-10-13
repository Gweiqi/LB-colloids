FORMATED TEXT INPUT
===================

*LB Colloids* computational fluid dynamic models can be parameterized
through the use of formatted text files that are parsed by internal
input control modules. Five separate data types are supported in
the input reader modules. Data types used during the
*LB-Colloids* input parameterization are strictly enforced. It is highly
recommended that the user familiarize themselves with data types before
attempting to parameterize and run *LB-Colloids.*

Description of Data Types
-------------------------

**STRING:** String type data within *LB-Colloids* is limited to file
names.

*Example usage:*

LBMODEL: mymodel.png

**INTEGER:** Integer type data.

*Example usage:*

NCOLS: 200

**FLOAT:** A floating point value or real number must be supplied.

*Example usage:*

AC: 1e-6

**BOOLEAN:** A value of True or False must be supplied.

*Example usage:*

PLOT: True

**LIST:** List type data combines data types to allow the user to
specify multiple values in the parameterization process. This data type
is used to set image boundary conditions using the solid and void
keywords.

*Example usage:*

SOLID: 255 223 200

**DICTIONARY:** Dictionary type data combines data types to use related
pairs of data in the parameterization process. This data type is limited
in usage to the optional concentration and valence keywords.

*Example usage:*

VALENCE: Na 1 Ca 2 Mg 2 Al 3

NAM File Inputs
---------------

The *LB-Colloid NAM file* is a simple formatted text file that is read
by run\_model.py. The NAM file contains relative or absolute file paths
to a Lattice Boltzmann configuration file and colloid transport
configuration files. The program run\_model.py looks within its current
directory for files ending with the suffix .nam, reads them, and
delegates the input and output control to the appropriate *LB-Colloids*
modules at runtime.

Input Structure:
****************

Input uses a block type structure. Block keywords begin a
parameterization section. When the parameterization section is complete,
the keyword end must be supplied for *LB-*\ Colloids to close the block
input reader. Within the input structure description, brackets
surrounding a parameter indicates that it is an optional parameter.

**LBMODEL (STRING):** Block keyword to designate the beginning of a
lattice Boltzmann parameterization block. This keyword is used to inform
run\_model.py that a lattice Boltzmann simulation will be performed.

    **LBCONFIG (STRING):** File name of the lattice Boltzmann
    configuration file that is used for simulation.

**END:** Block keyword to inform run\_model.py that it has finished
gathering parameters for the LBMODEL data block.

**[COLLOIDMODEL] (STRING):** Optional block keyword to designate the
beginning of a colloid transport parameterization block. The keyword is
used to inform run\_model.py that a colloid transport simulation will be
performed.

    **[COLLOIDCONFIG]\*n (STRING):** Optional file name of the colloid
    transport configuration file that is used for the simulation.
    Multiple colloid transport configuration files can be supplied for a
    multiple ionic strength model (MISM) and each filename must occupy a
    new line within the NAM file. This parameter is required if the
    COLLOIDMODEL keyword is used

**[END] (STRING):** Optional block keyword to inform run\_model.py that
it has finished gathering parameters for the COLLOIDMODEL data block.
This parameter is required if the COLLOIDMODEL keyword is used

Example NAM file:
*****************

|./images/NAM.png|

Lattice Boltzmann Configuration File
------------------------------------

The lattice Boltzmann configuration file uses a block input format to
parameterize D2Q9 lattice Boltzmann simulations. Four separate
blocks are available to the user. The MODEL PARAMETERS and IMAGE
PARAMETERS block must be supplied for basic model configuration. The
PERMEABILITY PARAMETERS block is optional, however the NITERS value
must be updated for a permeability model to be able to run to steady
state conditions. Configuration blocks can be added to the lattice
Boltzmann configuration file in any order as long as the proper
formatting requirements are fulfilled. Parameters within a configuration
block may be listed in any order as long as parameter requirements are
fulfilled. Only one parameter may be listed per line within the
configuration file.

Input Structure:
****************

**START MODEL PARAMETERS (STRING):** Required block keyword to begin
parameterization of the lattice Boltzmann model parameters data block.

    **LBMODEL (STRING):** Hdf5 file path and name corresponding to the
    lattice Boltzmann model. Lattice Boltzmann simulation results will
    be written to this Hdf5 file.

    **LBRES (FLOAT):** The lattice Boltzmann simulation resolution 
    parameter corresponds to the image resolution, in meters, of the 
    synthetic porous media or soil thin section for D2Q9 lattice 
    Boltzmann simulations.

    **[KERNEL] (STRING):** Lattice Boltzmann debug option that allows
    users to switch between simulating models with Python or FORTRAN as
    the base computational code. This option defaults to FORTRAN
    (recommended) if it is not supplied by the user. The FORTRAN kernel
    is approximately 100x faster than the Python kernel. Valid options
    are ‘python’ and ‘fortran’.

    **[PHYSICAL\_VISCOSITY] (FLOAT):** Optional parameter relating to
    the physical viscosity of the simulated fluid within a model run.
    This parameter is used in the dimensionalization process from
    non-dimensional lattice units to SI units. Default value is 8.9e-4
    Pa S which corresponds to the viscosity of water at 25\ :sup:`o` C.

    **[PHYSICAL\_RHO] (FLOAT):** Optional parameter relating to the
    physical density of the simulated fluid within a model run. This
    parameter is used in the dimensionalization process from
    non-dimensional lattice units to SI units. Default value is 997
    Kg/m\ :sup:`3` which corresponds to the viscosity of water at
    25\ :sup:`o` C.

**END MODEL PARAMETERS (STRING):** Required block keyword to end
parameterization of the lattice Boltzmann model parameters block

**START IMAGE PARAMETERS (STRING):** Required block keyword to begin the
parameterization of model domain boundary conditions using image
properties of the supplied thin section.

    **IMAGE (STRING):** Filename of the image representing the lattice Boltzmann simulation domain

    **SOLID (LIST [INTEGER]):** List of solid phase greyscale values as
    integers corresponding to the supplied model thin section

    **VOID (LIST [INTEGER]):** List of pore/void phase greyscale values
    as integers corresponding to the supplied model thin section

    **[BOUNDARY] (INTEGER):** Number of boundary layers to add to the
    top and bottom of the simulation domain. This parameter helps avoid
    local compression effects due to dead end pores connected to the top
    and bottom of the simulation domain. Default is 10 boundary layers.

    **[PLOT] (BOOLEAN):** Keyword parameter to display a plot of the
    binarized model domain before simulating fluid flow. This option is
    recommended during initial parameterization runs to ensure that
    fluid and solid phases are being properly represented within the
    domain. Default is False.

**END IMAGE PARAMETERS (STRING):** Required block keyword to end the
parameterization of the image parameters data block.

**START PERMEABILITY PARAMETERS (STRING):** Required keyword to begin
the lattice Boltzmann permeability model parameterization block.

    **[NITERS] (INTEGER):** Number of simulation time steps (number of
    iterations) applied to the lattice Boltzmann permeability
    simulation. Defaults to 1. It is highly recommended to change this
    value as one time step will not reach equilibrium conditions.

    **[RHO] (FLOAT):** Initial fluid density for lattice Boltzmann
    simulation. Defaults to 1.

    **[TAU] (FLOAT):** Lattice Boltzmann BGK relaxation time parameter.
    Defaults to 1. The acceptable range of TAU is from 0.5 to 1.5.

    **[GRAVITY] (FLOAT):** Gravity force applied to fluid as a body
    force. Gravity drives fluid flow in the simulation. Default is 1e-3.

**END PERMEABILITY PARAMETERS (STRING):** Required keyword to end the
lattice Boltzmann permeability model parameterization block

**[START OUTPUT CONTROL] (STRING):** Keyword to begin the lattice
Boltzmann output control parameterization block. This block contains
optional parameters that write updates to the terminal and save imagery
of lattice Boltzmann simulation progression.

    **[VERBOSE] (INTEGER):** Integer flag that indicates the number of
    time steps between lattice Boltzmann print statements to the
    terminal. This option updates user of the model’s progression. 
    Default is 100.

    **[IMAGE\_SAVE\_INTERVAL] (INTEGER):** Integer flag that indicates
    the time step interval to save model fluid velocity to a matplotlib 
    image. Default is None and no images are saved to disk.

    **[IMAGE\_SAVE\_NAME] (STRING):** Base name of images to save as
    output. This parameter is ignored unless the IMAGE\_SAVE\_INTERVAL 
    parameter is used. Default IMAGE\_SAVE\_NAME is “LB”.

    **[IMAGE\_SAVE\_FOLDER] (STRING):** Directory location to save
    simulation velocity images. This parameter is ignored unless the
    IMAGE\_SAVE\_INTERVAL parameter is used. Default IMAGE\_SAVE\_FOLDER
    is “~/Desktop/LBimages”.

    **[VMIN] (FLOAT):** VMIN controls the minimum boundary of velocity
    to plot on a lattice Boltzmann output image. This parameter is used
    to adjust the color scale for plotting only. VMIN can only be used
    if IMAGE\_SAVE\_INTERVAL is used. Default value is -0.010.

    **[VMAX] (FLOAT):** VMAX controls the maximum boundary of velocity
    to plot on a lattice Boltzmann output image. This parameter is used
    adjust the color scale for plotting only. VMAX can only be used if
    IMAGE\_SAVE\_INTERVAL is used. Default value is 0.0.

**[END OUTPUT CONTROL] (STRING):** Keyword to end the lattice Boltzmann
output control parameterization block.

Example Lattice Boltzmann configuration file:
*********************************************

|./images/LBCF.png|

Colloids Simulation Control File
--------------------------------

The colloids simulation control file uses a block input structure to
parameterize colloid simulations within the *LB-Colloids* system. Four
model parameter blocks are available to the user to parameterize colloid
models. The MODEL PARAMETERS block must be supplied in the colloid
simulation control file to run a model. The PHYSICAL PARAMETERS and
CHEMICAL PARAMETERS configuration blocks are optional, however it will
be necessary to specify some parameters within each of these blocks to
properly simulate experimental conditions for individual models.
Defaults within these blocks pertain to glass bead media and kaolinite
colloids. The OUTPUT PARAMETERS configuration block is optional and
contains useful parameters for later data analysis. Configuration blocks
can be added to the colloids simulation control file in any order as
long as the proper formatting requirements are fulfilled. Parameters
within a configuration block may be listed in any order as long as
parameter requirements are fulfilled. Only one parameter may be listed
per line within the configuration file.

Input Structure:
****************

**START MODEL PARAMETERS (STRING):** Keyword to begin the model
parameters input block. This input block contains necessary
parameterization information to run a basic colloid simulation.

    **LBMODEL (STRING):** HDF5 file name containing results from a
    steady state lattice Boltzmann simulation

    **LBRES (FLOAT):** Lattice Boltzmann grid resolution.

    **GRIDREF (FLOAT):** Grid refinement option which creates a bilinear
    interpolation of the lattice Boltzmann model domain. Colloid grid
    resolution is defined by

	.. math:: COLRES = \frac{\text{LBRES}}{\text{GRIDREF}}


    **ITERS (INTEGER):** Number of time steps the colloid simulation
    will run for

    **TIMESTEP (FLOAT):** Time step length. Reduction of time step
    length creates increased stability and greater accuracy, but longer 
    model run times. A recommended starting point is 1e-06 s.

    **NCOLS (INTEGER):** Number of colloids to introduce into the system
    in the initial time step.

    **[CONTINUOUS] (INTEGER):** Continuous is a flag that indicates
    multiple releases of colloids throughout the simulation. The value
    of continuous is the interval when additional colloids are released
    into the colloid simulation. Default is 0 (a single pulse of
    colloids released at the beginning of time step 1).

    **[AC] (FLOAT):** The colloid radius parameter defaults to 1e-6 m.

    **[RHO\_COLLOID] (FLOAT):** Colloid density parameter is optional
    and adjustable based on type of colloidal particle being simulated.
    Default value is 2650. Kg/m\ :sup:`3`

    **[TEMPERATURE] (FLOAT):** Fluid temperature parameter defaults to
    298.15 K.

**END MODEL PARAMETERS (STRING):** Keyword to end the model parameters
input block.

**[START PHYSICAL PARAMETERS] (STRING):** Keyword to that indicates the
beginning of the physical parameters input block.

    **[RHO\_WATER] (FLOAT):** Physical density of water. Default is 997.
    Kg/m\ :sup:`3`

    **[RHO\_COLLOID] (FLOAT):** Colloid density parameter. It is
    unnecessary to parameterize if RHO\_COLLOID has been added in the 
    model parameters input block. Defaults to 2650 Kg/m\ :sup:`3`

    **[VISCOSITY] (FLOAT):** Dynamic viscosity of the fluid phase within
    a colloid simulation. Default is 8.9e-4 Pa S.

    **[SCALE\_LB] (FLOAT):** Lattice Boltzmann velocity scaling factor.
    Exercise extreme caution in using this option.

**[END PHYSICAL PARAMETERS] (STRING):** Keyword to end the physical
parameters input block. This is a required parameter if START PHYSICAL
PARAMETERS keyword is used.

**[START CHEMICAL PARAMETERS] (STRING):** Keyword to that indicates the
beginning of the chemical parameters input block.

    **[I] (FLOAT):** The ionic strength of the fluid phase is calculated
    by the equation

	.. math:: I =\sum_{i = 1}^{n}{Z_{i}^{2}M_{i}}


    where Z\ :sub:`i` is the cation charge and M\ :sub:`i` is the
    molarity of each solution component i. The default ionic strength is
    set to 1e-3 M.

    **[ZETA\_SOLID] (FLOAT):** Bulk zeta potential of the solid phase
    within a colloid simulation. Default value is for glass bead media
    -60.9e-3 mV.

    **[ZETA\_COLLOID] (FLOAT):** Bulk zeta potential of colloids
    introduced into a colloid simulation. Default value is for kaolinite
    colloids -40.5e-3 mV.

    **[CONCENTRATION] (DICTIONARY [STRING, FLOAT]):** Optional parameter
    that allows the user to specify cation molarity pairs to
    parameterize the model instead of using ionic strength. This
    parameter must be supplied if VALENCE is used.

    **[VALENCE] (DICTIONARY [STRING, INTEGER]):** Optional parameter
    that allows the user to specify cation valence pairs to parameterize
    the model instead of using ionic strength. This parameter must be
    supplied if CONCENTRATION is used.

    **[LVDWST\_WATER] (FLOAT):** The van der Waals surface tension of
    the simulation fluid which is used to parameterize van der Waals
    interactions within the colloids simulation. Default value is
    21.8e-3 J/m\ :sup:`2` which corresponds to water.

    **[LVDWST\_COLLOID] (FLOAT):** The van der Waals surface tension of
    the simulated colloidal material which is used to parameterize van
    der Waals interactions within the colloids simulation. Default value
    is 39.9e-3 J/m\ :sup:`2` which corresponds to a kaolinite colloid.

    **[LVDWST\_SOLID] (FLOAT):** The van der Waals surface tension of
    the simulated solid phase which is used to parameterize van der
    Waals interactions within the colloids simulation. Default value is
    33.7e-3 J/m\ :sup:`2` for glass bead porous media.

    **[PSI+\_WATER] (FLOAT):** The electron-acceptor parameter of Lewis
    Acid Base surface tension for the simulation fluid. Default is
    25.5e-3 J/m\ :sup:`2` for water

    **[PSI+\_COLLOID] (FLOAT):** The electron-acceptor parameter of
    Lewis Acid Base surface tension for the colloid material. Default is
    0.4e-3 J/m\ :sup:`2` for a kaolinite colloid.

    **[PSI+\_SOLID] (FLOAT):** The electron-acceptor parameter of Lewis
    Acid Base surface tension for the solid phase. Default is 1.3e-3
    J/m\ :sup:`2` for a glass bead porous media.

    **[PSI-\_WATER] (FLOAT):** The electron-donor parameter of Lewis
    Acid Base surface tension for the simulation fluid. Default is
    25.5e-3 J/m\ :sup:`2` for water.

    **[PSI-\_COLLOID] (FLOAT):** The electron-donor parameter of Lewis
    Acid Base surface tension for the colloid material. Default is
    34.3e-3 J/m\ :sup:`2` for a kaolinite colloid.

    **[PSI-\_SOLID] (FLOAT):** The electron-donor parameter of Lewis
    Acid Base surface tension for the solid phase. Default is 62.2e-3
    J/m\ :sup:`2` for glass bead porous media.

    **[SHEER\_PLANE] (FLOAT):** Equivalent to the thickness of one layer
    of water molecules. Also referred to as thickness of the stern
    layer. Default is 3e-10 m.

    **[EPSILON\_R] (FLOAT):** The dielectric constant of water at the
    simulation temperature. The default value is 78.3 which corresponds
    to 25\ :sup:`o` C.

**[END CHEMICAL PARAMETERS]:** Keyword to end the chemical parameters
input block. This is a required parameter if START CHEMICAL PARAMETERS
keyword is used.

**[START OUTPUT CONTROL] (STRING):** Keyword to that indicates the
beginning of the output control input block.

    **[PRINT\_TIME] (INTEGER):** Integer flag that indicates the number
    of time steps between colloid simulations print statements to the
    terminal. Updates user of the model’s progression. Default is equal
    to the parameter ITERS.

    **[STORE\_TIME] (INTEGER):** Integer flag that indicates the number
    of time steps between internal storage functions within the colloid
    model. Increasing the value of this flag reduces memory consumption.
    This flag is also used to specify the interval for saving to a
    TIMESERIES file and output plotting.

    **[ENDPOINT] (STRING):** String flag that indicates an endpoint file
    should be saved. This option is highly recommended for use. The 
    endpoint variable should correspond to the name of the endpoint 
    file the user wishes to save.

    **[TIMESERIES] (STRING):** String flag that indicates a timeseries
    file should be saved. The save interval is indicated by the
    STORE\_TIME variable. The timeseries variable corresponds to the
    name of the timeseries file the user wishes to save.

    **[PATHLINE] (STRING):** String flag that indicates a pathline file
    should be saved. Colloid position is written to output at every time
    step using this option. The pathline variable corresponds to the
    name of the pathline file the user wishes to save.

    **[PLOT] (BOOLEAN):** Boolean flag that indicates a plot of colloid
    positions within the model domain be produced upon successful model
    completion. The plotting interval of colloid positions is set with 
    the STORE\_TIME flag. The plotted image will be saved as a <.png> 
    file with the same base name as was provided for ENDPOINT. If no 
    endpoint file was provided the figure will be displayed for the 
    user to save manually. Default is False.

    **[SHOWFIG] (BOOLEAN):** Boolean flag to indicate that the user
    wants to examine the figure before manually saving to file. SHOWFIG
    will only work when PLOT is True. Default is False.

    **[OVERWRITE] (BOOLEAN):** Boolean flag to overwrite data on
    existing HDF5 file. This flag is useful if the user does not wish to
    re-run lattice Boltzmann simulations while optimizing a colloid
    simulation.

**[END OUTPUT CONTROL] (STRING):** Keyword to end the output control
input block. This is a required parameter if START OUTPUT CONTROL
keyword is used.

Example colloid model configuration file:
*****************************************

|./images/CMCF.png|

.. |./images/NAM.png| image:: images/NAM.png
   :width: 6.50000in
   :height: 5.31250in
.. |./images/LBCF.png| image:: images/LBCF.png
   :width: 6.50000in
   :height: 5.20000in
.. |./images/CMCF.png| image:: images/CMCF.png
   :width: 6.50000in
   :height: 5.17708in
