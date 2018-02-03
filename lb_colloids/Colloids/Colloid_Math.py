"""
ColloidMath is the primary mathematics module for Colloid Simulations.
This module contains both Physical and Chemical formulations of colloid forces within a
porous media. The DLVO and ColloidColloid classes contain complex formulations
of chemical interaction forces. Other classes contain physical force calculations or
provide mathematical conversion from Force to a Velocity like unit that can be used
to recover the change in colloid position. Users should not have to call these classes
directly when running a model.

Basic examples of how these modules are called assume that a user has already provided input to
the lbIO.Config() module and the appropriate dictionaries have been built. For more information
on required parameters and keywords please inspect API documentation for each respective class.
Docstrings also provide basic mathematical relationships for each class.

>>> from lb_colloids import ColloidMath as cm
>>>
>>> grav = cm.Gravity(**PhysicalDict)
>>> grav.gravity  # returns the gravity force on a colloid
>>> bouy = cm.Bouyancy(**PhysicalDict)
>>> bouy.bouyancy  # returns the bouyancy force on a colloid
>>> gap = cm.Gap(xarr, yarr, **PhysicalDict)
>>> brownian = cm.Brownian(gap.f1, gap.f2, **PhysicalDict)
>>> brownian.brownian_x  # returns brownian force in the x-direction on a colloid
>>> brownian.brownian_y
>>> drag = cm.Drag(ux, uy, gap.f1, gap.f2, gap.f3, gap.f4, **PhysicalDict)
>>> drag.drag_x  # returns an array of drag forces in the x-direction
>>> dlvo = cm.DLVO(xarr, yarr, **ChemicalDict)
>>> dlvo.EDLx  # returns an array of electric double layer forces in the x-direction
>>> dlvo.LewisABy  # returns an array of lewis acid base forces in the y-direction
>>> dlvo.LVDWx  # returns an array of lifshitz-van der waals forces in the x-direction
>>> colcol = cm.ColloidColloid(xarr, **ChemicalDict)
>>> colcol.x_array  # returns an array of dlvo forces for colloid-colloid interactions
>>> colcol.update(Colloid.positions)  # updates the class to generate new colloid-colloid interaction arrays
"""

from .LB_Colloid import Colloid
import numpy as np
import sys
import copy


class ForceToVelocity:
    """
    Class that calculates a "velocity-like" value from force arrays

    Parameters:
    ----------
    :param np.ndarray forces: Array of forces felt by a colloid
    :keyword float ts: Physical time step value
    :keyword float rho_colloid: Colloid particle density, default :math:`2650 kg/m^3`
    :keyword float ac: colloid radius, default 1e-6 m

    Returns:
    -------
    :return: velocity (np.array, np.float) Array of "velocities" calculated from forces
    """
    def __init__(self, forces, **kwargs):

        params = {'rho_colloid': 2650., 'ac': 1e-6, 'ts': 1.}
        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]

        rho_colloid = params['rho_colloid']
        ac = params['ac']
        ts = params['ts']
        self.mass_colloid = (4. / 3.) * np.pi * (ac * ac * ac) * rho_colloid
        self.velocity = 0.5 * (forces * ts) / self.mass_colloid
        
        
class Velocity:
    """
    Class that dimensionalizes LB velocity from non-dimensional lattice Boltzmann units

    Parameters:
    ----------
    :param np.ndarray LBx: Array of Lattice Boltzmann velocities in the x-direction
    :param np.ndarray LBy: Array of Lattice Boltzmann velocities in the y-direction
    :keyword float ts: Time step value, default is 1.
    :keyword float scale_lb: Scale the dimensionalized velocity from lattice Boltzmann. Use with caution. Default is 1
    :param float velocity_factor: LB to physical velocity conversion factor. Default is 1

    Returns:
    -------
    :return: xvelocity (np.array, np.float) array of dimensionalized velocities in the x-direction
    :return: yvelocity (np.array, np.float) array of dimensionalized velocities in the y-direction
    """
    def __init__(self, LBx, LBy, velocity_factor, **kwargs):

        params = {'lb_timestep': 1e-5, 'ts': 1, 'scale_lb': 1.}

        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]

        ts = params['ts']

        self.xvelocity = LBx * velocity_factor * params['scale_lb']
        self.yvelocity = LBy * velocity_factor * params['scale_lb']
        

class Gravity:
    """
    Class to generate the estimated gravitational force experienced by a colloid

    .. math::

        F^{G} = \\frac{-4 \pi a_{c}^{3} \\rho_{c} g}{3}

    Parameters:
    ----------
    :keyword float rho_colloid: Particle density of a colloid in :math:`kg/m^3`. Default is 2650.
    :keyword float ac: colloid radius in m. Default is 1e-6

    Returns:
    -------
    :return: gravity (float) Gravitational force that a colloid experiences
    """
    def __init__(self, **kwargs):

        params = {'rho_colloid': 2650., 'ac': 1e-6}
        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]

        ac = params['ac']
        rho_colloid = params['rho_colloid']
        
        self.colloid_mass = (4./3.)*np.pi*(ac*ac*ac)*rho_colloid
        self.gravity = (self.colloid_mass*-9.81)


class Bouyancy:
    """
    Class to estimate the gravitational force experienced by a colloid. Gravity
    is applied as a positive value to maintain vector direction.

    .. math::

        F^{b} = \\frac{4 \pi a_{c}^{3} \\rho_{w} g}{3}

    Parameters:
    ----------
    :keyword flaot rho_water: density of water :math:`kg/m^3`. Default is 997.
    :keyword float rho_colloid: particle density of a colloid in :math:`kg/m^3`. Default is 2650.
    :keyword float ac: colloid radius in m. Default is 1e-6.

    Returns:
    -------
    :return: bouyancy (float) Bouyancy force that a colloid experiences
    """
    def __init__(self, **kwargs):
        params = {'rho_water': 997., 'rho_colloid': 2650., 'ac': 1e-6}
        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]
        
        rho_water = params['rho_water']
        rho_colloid = params['rho_colloid']
        ac = params['ac']
        
        self.water_mass = (4./3.)*np.pi*(ac*ac*ac)*rho_water
        self.bouyancy = self.water_mass * 9.81


class Brownian:
    """
    Class to estimate brownian forces on colloids. Uses the relationships outlined in Qui et. al. 2010
    where

    .. math::
        F_{x}^{B} = \\xi \sqrt{\\frac{2D_{0}}{f_{1}dt}}G(0,1)

        F_{y}^{B} = \\xi \sqrt{\\frac{2D_{0}}{f_{4}dt}}G(0,1)

    Parameters:
    ----------
    :param np.ndarray f1: Drag force correction term [Gao et. al. 2010. Computers and Math with App]
    :param np.ndarray f4: Drag force correction term [Gao et. al. 2010]
    :keyword float ac: Colloid radius. Default 1e-6
    :keyword float viscosity: Dynamic viscosity of water. Default 8.9e-4 Pa S.
    :keyword float T: Absolute temperature in K. Default is 298.15

    Returns:
    -------
    :return: brownian_x: (np.ndarray) array of browian (random)
        forces in the x direction [Qiu et. al 2011.]
    :return: brownian_y: (np.ndarray) array of browian (random)
        forces in the y direction [Qiu et. al 2011.]
    """
    def __init__(self, f1, f4, **kwargs):

        # todo: update brownian motion to include the timestep!!!!
        params = {'viscosity': 8.9e-4, 'ac': 1e-6, 'T': 298.15}
        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]
                
        self.ac = params['ac']
        self.ts = params['ts']
        self.viscosity = params['viscosity']
        self.boltzmann = 1.38e-23
        self.epsilon = 6. * np.pi * self.viscosity * self.ac
        self.T = params['T']
        self.diffusive = (self.boltzmann * self.T) / self.epsilon
        self.brownian_x = self.Brown_xforce(self.epsilon, self.diffusive, f4)
        self.brownian_y = self.Brown_yforce(self.epsilon, self.diffusive, f1)

    def Brown_xforce(self, epsilon, diffusivity, f4):
        mu, sigma = 0, 1.
        Fbt = epsilon * np.sqrt(((2 * diffusivity)/(f4 * self.ts))) * np.random.normal(mu, sigma, (len(f4), len(f4[0])))
        return Fbt

    def Brown_yforce(self, epsilon, diffusivity, f1):
        mu, sigma = 0, 1.
        Fbn = epsilon * np.sqrt(((2 * diffusivity)/(f1 * self.ts))) * np.random.normal(mu, sigma, (len(f1), len(f1[0])))
        return Fbn


class Drag:
    """
    Class to calculate colloidal drag forces from fluid velocity arrays. Based from calculations
    outlined in Gao et, al 2010 and Qui et. al. 2011.

    .. math::
        F_{x}^{D} = \\frac{\\xi}{f_{4}} (f_{3}u_{x} - V_{x})

        F_{y}^{D} = \\xi (f_{2} u_{y} - \\frac{V_{y}}{f_{1}})

    Parameters:
    ----------
    :param np.ndarray ux: fluid velocity in the x-direction
    :param np.ndarray uy: fluid velocity in the y-direction
    :param np.ndarray Vx: colloid velocity in the x-direction
    :param np.ndarray Vy: colloid velocity in the y-direction
    :param np.ndarray f1: Hydrodynamic force correction term [Gao et. al. 2010.]
    :param np.ndarray f2: Hydrodynamic force correction term [Gao et. al. 2010.]
    :param np.ndarray f3: Hydrodynamic force correction term [Gao et. al. 2010.]
    :param np.ndarray f4: Hydrodynamic force correction term [Gao et. al. 2010.]
    :keyword float ac: Colloid radius. Default is 1e-6 m
    :keyword float viscosity: Dynamic fluid viscosity of water. Default 8.9e-4 Pa S
    :keyword float rho_colloid: Colloid particle density. Default :math:`2650 kg/m^3`
    :keyword float rho_water: Water density. Default :math:`997 kg/m^3`

    Returns:
    -------
    :return: drag_x (np.ndarray) non-vectorized drag forces in the x-direction
    :return: drag_y: (np.ndarray) non-vectorized drag forces in the y-direction
    """
    def __init__(self, ux, uy, f1, f2, f3, f4, **kwargs):

        params = {'ac': 1e-6, 'viscosity': 8.9e-4, 'rho_colloid': 2650., 'rho_water': 997.,
                  'T': 298.15, 'ts': 1.}
        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]        
        
        self.ac = params['ac']
        self.viscosity = params['viscosity']
        self.rho_water = params['rho_water']
        self.rho_colloid = params['rho_colloid']
        self.epsilon = 6. * np.pi * self.viscosity * self.ac
        self.Vcol = -((self.rho_colloid - self.rho_water)*((2*self.ac)**2)*9.81)/(18*self.viscosity)
        # todo: update this all to a fortran routine that is called each iteration. Replace VCol with stored value!
        self.drag_x = self.drag_xforce(ux, self.Vcol, self.epsilon, f3, f4)
        self.drag_y = self.drag_yforce(uy, self.Vcol, self.epsilon, f1, f2)

        self.all_physical_params = copy.copy(params)
        
    def drag_xforce(self, ux, Vx, epsilon, f3, f4):
        Fdt = (epsilon / f4) * ((f3 * ux) - Vx)
        return Fdt

    def drag_yforce(self, uy, Vy, epsilon, f1, f2):
        Fdn = epsilon * ((f2 * uy) - (Vy / f1))
        return Fdn

        
class Gap:
    """
    Class that calculates the non-dimensional gap distance between colloid and surface.

    This class also calculates hydrodynamic force correction terms outlined in Gao et. al. 2010.
    Note: Passing a np.nan value into here can return an overflow warning!

    .. math::
        f_{1}(\\bar{h}) = 1.0 - 0.443 exp(-1.299\\bar{h}) - 0.5568 exp(-0.32\\bar{h}^{0.75})

    .. math::

        f_{2}(\\bar{h}) = 1.0 + 1.455 exp(-1.2596\\bar{h}) - 0.7951 exp(-0.56\\bar{h}^{0.50})

    .. math::

        f_{3}(\\bar{h}) = 1.0 - 0.487 exp(-5.423\\bar{h}) - 0.5905 exp(-37.83\\bar{h}^{0.50})

    .. math::

        f_{4}(\\bar{h}) = 1.0 - 0.35 exp(-0.25\\bar{h}) - 0.40 exp(-10\\bar{h})

    Parameters:
    ----------
    :param np.ndarray xarr: Array of x-distances to nearest solid surface
    :param np.ndarray yarr: Array of y-distances to nearest solid surface
    :keyword float ac: Radius of a colloid. Default is 1e-6

    Returns:
    -------
    :return: f1 (np.ndarray) Drag force correction term [Gao et al 2010]
    :return: f2 (np.ndarray) Drag force correction term [Gao et al 2010]
    :return: f3 (np.ndarray) Drag force correction term [Gao et al 2010]
    :return: f4 (np.ndarray) Drag force correction term [Gao et al 2010]
    """
    def __init__(self, xarr, yarr, **kwargs):


        params = {'ac': 1e-6}
        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]
        
        self.ac = params['ac']
        self.yhbar = np.abs(yarr/self.ac)
        self.xhbar = np.abs(xarr/self.ac)
        self.f1 = self.set_f1(self.yhbar)
        self.f2 = self.set_f2(self.yhbar)
        self.f3 = self.set_f3(self.xhbar)
        self.f4 = self.set_f4(self.xhbar)
        
    def set_f1(self, yhbar):
        f1 = 1.0 - 0.443 * np.exp(yhbar * -1.299) - 0.5568 * np.exp((yhbar ** 0.75) * -0.32)
        return f1

    def set_f2(self, yhbar):
        f2 = 1.0 + 1.455 * np.exp(yhbar * -1.259) + 0.7951 * np.exp((yhbar ** 0.50) * -0.56)
        return f2

    def set_f3(self, xhbar):
        f3 = 1.0 - 0.487 * np.exp(xhbar * -5.423) - 0.5905 * np.exp((xhbar ** 0.50) * -37.83)
        return f3

    def set_f4(self, xhbar):
        f4 = 1.0 - 0.35 * np.exp(xhbar * -0.25) - 0.40 * np.exp(xhbar * -10.)
        return f4


class DLVO:
    """
    Class method to calculate vectorized DLVO force arrays for colloid surface interaction using
    methods outlined in Qui et. al. 2011 and Liang et. al. 2008? *Check this later*

    Parameterization of this class is handled primary through the ChemistryDict by **kwargs

    Mathematics used in calcuation of DLVO interaction energies are:

    .. math::

        \\frac{1}{\kappa} = (\\frac{\epsilon_{r} \epsilon_{0} k T}{e^{2} N_{A} I^{*}})^{\\frac{1}{2}}

    .. math::

        \Phi^{EDL} = \pi \epsilon_{0} \epsilon_{r} a_{c}
        (2 \psi_{s} \psi_{c}
        ln(\\frac{1 + exp(-\kappa h)}{1 - exp(-\kappa h)})
        + (\psi_{s}^{2} + \psi_{c}^{2})
        ln(1 - exp(-2 \kappa h)))

    Parameters:
    -------
    :param np.ndarray xarr: Physical distance from solid boundaries in the x direction
    :param np.ndarray yarr: Physical distance from solid boundaries in the y direction
    :keyword dict valence: Valences of all species in solution. (Optional)
    :keyword dict concentration: Concentration of all species in solution (Optional)
    :keyword float zeta_colloid: Measured_zeta potential of colloid (Reccomended).
        Default -40.5e-3 Na-Kaolinite Colloid [Chorom 1995. Eur. Jour. of Soil Science]
    :keyword float zeta_surface: Bulk_zeta potential of porous media (Reccomended).
        Default -60.9e-3 Glass bead media [Ducker 1992, Langmuir V8]
    :keyword float I: Ionic strength of simulated solution (Reccomended). Default 1e-3 M
    :keyword float ac: Colloid radius in meters. Default 1e-6 m.
    :keyword float epsilon_r: Relative dielectric permativity of water. (Optional)
        Default 78.304 @ 298 K [Malmberg and Maryott 1956. Jour. Res. Nat. Beau. Std. V56(1)

    :keyword float sheer_plane: Equivelent to the thickness of one layer of water molecules. (Optional)
        Default 3e-10 m [Interface Science and Technology, 2008. Volume 16 Chapter 3]
    :keyword float T: Temperature of simulation fluid. Default 298.15 k
    :keyword float lvdwst_colloid: Lifshits-van der Waals surface tension component from colloid. (Reccomended)
        Default is 39.9e-3 J/m**2 [Giese et. al. 1996, Jour. Disp. Sci. & Tech. 17(5)]
    :keyword float lvdwst_solid: Lifshits-van der Waals surface tension component from solid. (Reccomended)
        Default is 33.7e-3 J/m**2 [Giese et. al. 1996]
    :keyword float lvdwst_water: Lifshits-van der Waals surface tension component from water. (Reccomended)
        Default is 21.8e-3 J/m**2 [Interface Science and Technology, 2008. V16(2)]
    :keyword float psi+_colloid: Lewis acid base electron acceptor parameter. (Reccomended)
        Default is 0.4e-3 J/m**2 [Giese et. al. 1996]
    :keyword float psi+_solid: Lewis acid base electron acceptor parameter. (Reccomended)
        Default is 1.3e-3 J/m**2 [Giese et. al. 1996]
    :keyword float psi+_water: Lewis acid base electron acceptor parameter. (Reccomended)
        Default is 25.5e-3 J/m**2 [Interface Science and Technology, 2008. V16(2)]
    :keyword float psi-_colloid: Lewis acid base electron donor parameter. (Reccomended)
        Default is 34.3e-3 J/m**2 [Giese et. al. 1996]
    :keyword float psi-_solid: Lewis acid base electron donor parameter. (Reccomended)
        Default is 62.2e-3 J/m**2 [Giese et. al. 1996]
    :keyword float psi-_water: Lewis acid base electron donor parameter. (Reccomended)
        Default is 25.5e-3 J/m**2 [Interface Science and Technology, 2008. V16(2)]
    :keyword np.ndarray xvArr: Array of vector directions.This array is applied to properly represent attractive
        and repulsive forces
    :keyword np.ndarray yvArr: Array of vector directions.This array is applied to properly represent attractive
        and repulsive forces

    Return:
    ------
    :return: EDLx (np.ndarray) vectorized np.array of electric-double-layer force values in the x-direction
    :return: EDLy (np.ndarray) vectorized np.array of electric-double-layer force values in the y-direction
    :return: LVDWx (np.ndarray) vectorized np.array of lifshitz-van-der-walls force values in the x-direction
    :return: LVDWy (np.ndarray) vectorized np.array of lifshitz-van-der-walls force values in the y-direction
    :return: LewisABx (np.ndarray) vectorized np.array of lewis acid base force values in the x-direction
    :return: LewisABy (np.ndarray) vectorized np.array of lewis acid base force values in the y-direction
    """
    def __init__(self, xarr, yarr, **kwargs):

        params = {'concentration': {'Na': 10e-4}, 'adjust_zeta': False, 'I_initial': False, 'I': 10e-4, 'ac': 1e-6,
                  'epsilon_r': 78.304, 'valence': {'Na': 1.}, 'sheer_plane': 3e-10, 'T': 298.15,
                  'lvdwst_water': 21.8e-3, 'lvdwst_colloid': 39.9e-3, 'lvdwst_solid': 33.7e-3, 'zeta_colloid': -40.5e-3,
                  'zeta_solid': -60.9e-3, 'psi+_colloid': 0.4e-3, 'psi-_colloid': 34.3e-3, 'psi+_water': 25.5e-3,
                  'psi-_water': 25.5e-3, 'psi+_solid': 1.3e-3, 'psi-_solid': 62.2e-3, 'rho_colloid': 2650.}

        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]

        self.rho_colloid = params['rho_colloid']
        self.epsilon_0 = 8.85e-12
        self.epsilon_r = params['epsilon_r']
        self.ac = params['ac']
        self.e = 1.6e-19
        self.valence = params['valence']
        self.concentration = params['concentration'] 
        self.boltzmann = 1.38e-23
        self.stern_z = params['sheer_plane'] 
        self.T = params['T']
        self.zeta_colloid = params['zeta_colloid']
        self.zeta_solid = params['zeta_solid']
        self.lvdwst_water = params['lvdwst_water']
        self.lvdwst_colloid = params['lvdwst_colloid']
        self.lvdwst_solid = params['lvdwst_solid']
        self.eplus_water = params['psi+_water']
        self.eplus_colloid = params['psi+_colloid']
        self.eplus_solid = params['psi+_solid']
        self.eneg_water = params['psi-_water']
        self.eneg_colloid = params['psi-_colloid']
        self.eneg_solid = params['psi-_solid']
        self.xvArr = params['xvArr']*-1
        self.yvArr = params['yvArr']*-1
        self.hamaker = None

        self.all_chemical_params = copy.copy(params)

        self.__resolution = params['lbres']/params['gridref']
        
        if params['I']:
            self.ionic_strength = 2 * params['I']  # 2I is what is used in the debye equation
        else:
            self.ionic_strength = self.ionic(params['valence'], params['concentration'])
        
        self.k_debye = self.debye(self.epsilon_0, self.epsilon_r, self.boltzmann, self.T, self.e,
                                  self.ionic_strength)
                                         
        self.colloid_potential = self._colloid_potential(self.zeta_colloid, self.ac, self.k_debye, self.stern_z)
        self.surface_potential = self._surface_potential(self.zeta_solid, self.k_debye, self.stern_z)

        # Calculate the chemical potential
        self.EDLx = self._EDL_energy(self.epsilon_0, self.epsilon_r, self.ac, self.colloid_potential,
                                     self.surface_potential, self.k_debye, xarr)/xarr*self.xvArr

        self.EDLy = self._EDL_energy(self.epsilon_0, self.epsilon_r, self.ac, self.colloid_potential,
                                     self.surface_potential, self.k_debye, yarr)/yarr*self.yvArr
        
        self.LVDWx = self._Lifshitz_van_der_Walls(xarr, self.ac, self.lvdwst_water, self.lvdwst_colloid,
                                                  self.lvdwst_solid)/xarr*self.xvArr

        self.LVDWy = self._Lifshitz_van_der_Walls(xarr, self.ac, self.lvdwst_water, self.lvdwst_colloid,
                                                  self.lvdwst_solid)/yarr*self.yvArr

        self.LewisABx = self._lewis_acid_base(xarr, self.ac, self.eplus_colloid, self.eplus_solid, self.eplus_water,
                                              self.eneg_colloid, self.eneg_solid, self.eneg_water)/xarr*self.xvArr

        self.LewisABy = self._lewis_acid_base(yarr, self.ac, self.eplus_colloid, self.eplus_solid, self.eplus_water,
                                              self.eneg_colloid, self.eneg_solid, self.eneg_water)/yarr*self.yvArr
        self._combined_hamaker_constant()

        self.attractive_x = self._combined_lvdw_lewis_ab(xarr)/xarr * self.xvArr
        self.attractive_y = self._combined_lvdw_lewis_ab(yarr)/yarr * self.yvArr

    def ionic(self, valence, concentration):
        """
        Calculates the 2*I from user supplied valence and concentraitons

        .. math::

            I^{*} = \sum_{i} Z_{i}^{2} M_{i}

        Parameters:
        ----------
        :param dict valence: Dictionary of chemical species, valence
        :param dict concentration: Dictionary of chemical species, concentration

        Returns:
        -------
        :return: I (float) 2*ionic stength
        """
        I = 0
        for key in valence:
            I += (float(concentration[key])*(float(valence[key])**2))
        return I

    def debye(self, epsilon_0, epsilon_r, kb, T, e, ionic_strength):
        """
        Method to calculate Debye length

        Parameters:
        ----------
        :param float epsilon_0: Permativity of a vacuum
        :param float epsilon_r: Relative permativity of water
        :param float kb: Boltzmann constant
        :param float T: fluid temperature in K
        :param float e: electron charge
        :param float ionic_strength: 2*ionic strength

        Returns:
        -------
        :return: Debye length (float)
        """
        NA = 6.02e23
        k_inverse = np.sqrt((epsilon_0*epsilon_r*kb*T)/(e*e*NA*ionic_strength))
        return 1./k_inverse

    def _colloid_potential(self, zeta, ac, kd, z):
        """
        Calculates the surface potential on a colloid

        Parameters:
        ----------
        :param float zeta: Zeta potential of colloid
        :param float ac: Colloid radius
        :param float kd: Debye length
        :param float z: Thickness of the sheer plane (stern layer)

        Returns:
        -------
        :return: (float) colloid surface potential
        """
        potential = zeta*(1.+(z/ac))*np.exp(kd*z)
        return potential

    def _surface_potential(self, zeta, kd, z):
        """
        Calculates the surface potential of the solid phase

        Parameters:
        ----------
        :param float zeta: Zeta potential of Solid phase
        :param float kd: Debye length
        :param float z: Thickness of the sheer plane (stern layer)

        Returns:
        -------
        :return: (float) Solid phase surface potential
        """
        potential = zeta*np.exp(kd*z)
        return potential
    
    def _EDL_energy(self, E0, Er, ac, cp, sp, kd, arr):
        """
        Parameters:
        ----------
        E0: (float) dilectric permativity in a vacuum
        Er: (float) fluid permativity
        ac: (float) colloid radius
        cp: (float) colloid potential
        kd: (float) debye length
        arr: (np.array: np.float) array of distances from surfaces

        Output:
        -------
        EDL: (np.array: np.float) array of EDL energies in relation to porous surfaces

        Note:
        -----
        Mathematical calcualtion is broken in three sections for ease of programming
        """

        edl0 = np.pi*E0*Er*ac
        edl1 = 2.*sp*cp
        edl2 = np.log((1. + np.exp(-kd*np.abs(arr)))/(1. - np.exp(-kd*np.abs(arr))))
        edl3 = sp*sp + cp*cp
        edl4 = np.log(1. - np.exp(-2.*kd*np.abs(arr)))

        edl = edl0*(edl1*edl2 + edl3*edl4)

        return edl

    def _adjust_zeta_colloid(self, potential, ac, kd, z):
        zeta = potential/((1. + (z/ac))*np.exp(kd*z))
        return zeta

    def _adjust_zeta_surface(self, potential, kd, z):
        zeta = potential/(np.exp(kd*z))
        return zeta

    def _combined_hamaker_constant(self):
        """
        Method to calculate the hamaker constant for surface-colloid interactions
        based on Israelachvili 1991
        """
        s_ah = self.surface_potential * (24 * np.pi * 0.165e-9 ** 2)
        c_ah = self.colloid_potential * (24 * np.pi * 0.165e-9 ** 2)

        self.hamaker = np.sqrt(s_ah * c_ah)

    def _combined_lvdw_lewis_ab(self, arr):
        """
        Method to calculate the combined attractive force profile based on liang et. al.
        instead of using vdw and lewis acid base profiles seperately
        Parameters:
        ----------
        :param np.ndarray arr: distance array

        :return: (np.ndarray) attractive force profile for porous media
        """
        lvdw_lab0 = -self.hamaker / 6.
        lvdw_lab1 = (self.ac / arr) + (self.ac / (arr + (2.* self.ac)))
        lvdw_lab2 = np.log(arr / (arr + self.ac))

        return lvdw_lab0 * (lvdw_lab1 + lvdw_lab2)

    # todo: remove attractive force calculations and replace with Hamaker & Liang calcs.
    def _Lifshitz_van_der_Walls(self, arr, ac, vdw_st_water, vdw_st_colloid, vdw_st_solid):
        """
        Parameters:
        ----------
        arr: (np.array, np.float) array of distances from solid surfaces
        ac: (float) colloid radius
        vdw_st_water: (float) vdW surface tension of water
        vdw_st_colloid: (float) vdW surface tension of colloid
        vdw_st_solid: (float) vdW surface tension (bulk) of solid phase

        constant:
        --------
        h0: contact plane between colloid and surface {Interface Science and Technology, 2008. Volume 16. Chapter 3}

        Returns:
        -------
        lvdw: (np.array, np.float) array of lifshitz_vdW interaction energies
        """
        
        h0 = 1.57e-10

        lvdw0 = -4.*np.pi*(h0*h0/arr)*ac
        lvdw1 = np.sqrt(vdw_st_water) - np.sqrt(vdw_st_solid)
        lvdw2 = np.sqrt(vdw_st_water) - np.sqrt(vdw_st_colloid)

        lvdw = lvdw0*lvdw1*lvdw2
        return lvdw

    def _lewis_acid_base(self, arr, ac, eplus_colloid, eplus_solid, eplus_water, eneg_colloid,
                         eneg_solid, eneg_water):
        """
        Parameters:
        ----------
        arr: (np.array, np.float) array of distances from solid surfaces
        e_plus_*: (float) electron acceptor parameter for each specific phase
        e_minus_*: (float) electron donor parameter for each specific phase

        Constants:
        ----------
        h0: contact plane between colloid and surface  {Interface Science and Technology, 2008. Volume 16. Chapter 3}
        chi: water decay length {van Oss 2008}

        Returns:
        -------
        lab: (np.array, np.float) array of lewis acid base interaction energies
        
        """
        h0 = 1.57e-10
        chi = 0.6e-10

        lab0 = -4.*np.pi*h0*ac
        lab1 = np.exp((h0-arr)/chi)
        lab2 = np.sqrt(eplus_water)*(np.sqrt(eneg_colloid) + np.sqrt(eneg_solid) - np.sqrt(eneg_water))
        lab3 = np.sqrt(eneg_water)*(np.sqrt(eplus_colloid) + np.sqrt(eplus_solid) - np.sqrt(eplus_water))
        lab4 = np.sqrt(eplus_colloid*eneg_solid)
        lab5 = np.sqrt(eneg_colloid*eplus_solid)

        lab = lab0*lab1*(lab2+lab3-lab4-lab5)
        return lab


class ColloidColloid(object):
    """
    The ColloidColloid class is used to calculate colloid-colloid interaction forces
    using the formulations presented in Liang 2008, Qui 2012, and Israelichevi 1996.
    Attractive forces are based on the Liang & Israelichevi formulation. Electric
    doulbe layer forces are calculated using Qui et. al. 2012.

    The ColloidColloid object also provides methods to update ColloidColloid force
    array fields during model streaming.

    Colloid colloid interaction energies are calculated via:

    .. math::
        \Phi^{EDL} = 32 \pi \epsilon_{0} \epsilon_{r} a_{c}
        (\\frac{kT}{Ze})^{2} * [tanh(\\frac{Ze\psi_c}{4kT})]^{2}
        * exp(-\kappa h)

    .. math::

        A_{H} = 384 \pi \\frac{\psi_{c}^{2} h k T I^{*}}{\kappa^{2}} exp(- \kappa h)

    .. math::

        \Phi^{A} = - \\frac{A_{H}}{6}[\\frac{2a_{c}^{2}}{h^{2} + 4a_{c}h} +
        \\frac{2a_{c}^{2}}{(h + 2a_{c})^{2}} + ln(1 - \\frac{4a_{c}^{2}}{(h + 2a_{c})^{2}})]

    Parameters:
    ----------
    :param np.ndarray arr: A np.ndarray that represents the shape of the colloid
        domain
    :param float resolution: Colloid model resolution

    :keyword dict valence: Valences of all species in solution. (Optional)
    :keyword dict concentration: Concentration of all species in solution (Optional)
    :keyword float zeta_colloid: Measured_zeta potential of colloid (Reccomended).
        Default -40.5e-3 Na-Kaolinite Colloid [Chorom 1995. Eur. Jour. of Soil Science]
    :keyword float zeta_surface: Bulk_zeta potential of porous media (Reccomended).
        Default -60.9e-3 Glass bead media [Ducker 1992, Langmuir V8]
    :keyword float I: Ionic strength of simulated solution (Reccomended). Default 1e-3 M
    :keyword float ac: Colloid radius in meters. Default 1e-6 m.
    :keyword float epsilon_r: Relative dielectric permativity of water. (Optional)
        Default 78.304 @ 298 K [Malmberg and Maryott 1956. Jour. Res. Nat. Beau. Std. V56(1)
    :keyword float sheer_plane: Equivelent to the thickness of one layer of water molecules. (Optional)
        Default 3e-10 m [Interface Science and Technology, 2008. Volume 16 Chapter 3]
    :keyword float T: Temperature of simulation fluid. Default 298.15 k
    """
    def __init__(self, arr, **kwargs):

        self.__params = {'concentration': False, 'adjust_zeta': False, 'I_initial': False,
                         'I': 10e-4, 'ac': 1e-6, 'epsilon_0': 8.85e-12 , 'epsilon_r': 78.304, 'valence': {'Na': 1.},
                         'sheer_plane': 3e-10, 'T': 298.15, 'lvdwst_water': 21.8e-3, 'lvdwst_colloid': 39.9e-3,
                         'lvdwst_solid': 33.7e-3, 'zeta_colloid': -40.5e-3, 'zeta_solid': -60.9e-3,
                         'psi+_colloid': 0.4e-3, 'psi-_colloid': 34.3e-3, 'psi+_water': 25.5e-3,
                         'psi-_water': 25.5e-3, 'psi+_solid': 1.3e-3, 'psi-_solid': 62.2e-3, 'kb': 1.38e-23,
                         'e': 1.6e-19, 'rho_colloid': 2650.}

        for kwarg, value in kwargs.items():
            self.__params[kwarg] = value

        self.__arr = arr
        self.__xarr = np.zeros(arr.shape)
        self.__yarr = np.zeros(arr.shape)
        self.__debye = False
        self.__colloid_potential = False
        self.__ionic_strength = False
        self.__resolution = copy.copy(self.__params['lbres'])/self.__params['gridref']
        self.__pos = []
        self.__x_distance = False
        self.__y_distance = False
        self.__x = False
        self.__y = False
        self.__center = False

    def __reset(self):
        """
        Resets the calculation arrays
        """
        self.__xarr = np.zeros(self.__arr.shape)
        self.__yarr = np.zeros(self.__arr.shape)
        self.__pos = []
        self.__x = False
        self.__y = False

    def __get_colloid_positions(self):
        """
        Get the specific x, y positions of each colloid in the system

        Parameters:
        -----------
            colloids: (list, <class: Colloids.LB_Colloid.Colloid)

        Returns:
        --------
            pos: (list) list of colloid positions within the model space
        """
        self.__pos = Colloid.positions
        return self.__pos

    def update(self, colloids):
        """
        Updates the colloidal positions and force arrays for the system

        Parameters:
        ----------
        :param list colloids: (list, <class: Colloids.LB_Colloid.Colloid)
        """
        self.__reset()
        # self.positions

    @property
    def x_array(self):
        """
        Property method to generate the full x force array for colloid-colloid interaction
        """
        return self.__get_full_dlvo_array("x")

    @property
    def y_array(self):
        """
        Property method to generate the full y force array for colloid-colloid interaction
        """
        return self.__get_full_dlvo_array("y")

    @property
    def x(self):
        """
        Property method to generate the x force array for colloid-colloid interaction
        """
        if isinstance(self.__x, bool):
            self.__x = self.__dlvo_interaction_energy("x")
        return self.__x

    @property
    def y(self):
        """
        Property method to generate or return the y force array for colloid-colloid interaction
        """
        if isinstance(self.__y, bool):
            self.__y = self.__dlvo_interaction_energy("y")
        return self.__y

    @property
    def x_distance_array(self):
        """
        Generates an angular distance array in the x direction.
        """
        if isinstance(self.__x_distance, bool):
            self.__x_distance = self.__angular_array("x")
        return self.__x_distance

    @property
    def y_distance_array(self):
        """
        Generates an angular distance array in the y direction
        """
        if isinstance(self.__y_distance, bool):
            self.__y_distance = self.__angular_array("y")
        return self.__y_distance

    @property
    def positions(self):
        """
        Property method to generate colloid positions if they are not stored yet
        """
        if not self.__pos:
            self.__get_colloid_positions()
        return self.__pos

    @property
    def ionic_strength(self):
        """
        Property method to calculate ionic_strength on the fly
        """
        if not self.__params['concentration']:
            return self.__params['I']*2
        else:
            I = 0
            for key in self.__params['concentration']:
                I += (float(self.__params['concentration'][key])
                      * (float(self.__params['valence'][key]) ** 2))
            return I

    @property
    def debye(self):
        """
        Property method to calculate the debye length on the fly
        """
        if isinstance(self.__debye, bool):
            na = 6.02e23
            k_inverse = np.sqrt((self.__params['epsilon_0']*self.__params['epsilon_r']
                                *self.__params['kb']*self.__params['T'])/
                                (self.__params['e']*self.__params['e']*na*self.ionic_strength))
            self.__debye = 1./k_inverse
        return self.__debye

    @property
    def colloid_potential(self):
        """
        Property method that generates colloid potential
        """
        if isinstance(self.__colloid_potential, bool):
            self.__colloid_potential = self.__params['zeta_colloid']*(1. +
                                       (self.__params['sheer_plane']/self.__params['ac']))\
                                        *np.exp(self.debye*self.__params['sheer_plane'])
        return self.__colloid_potential

    def __get_full_dlvo_array(self, arr_type):
        """
        Handler definition to call subroutes to generate dvlo_force_array

        Parameters:
            arr_type: (str) x direction or y direction , "x", "y"

        Returns:
            dvlo: (np.ndarray) full array of dlvo interaction forces from colloids
        """
        if arr_type.lower() == "x":
            arr = self.__xarr
            dlvo_colloid = self.x
        elif arr_type.lower() == "y":
            arr = self.__yarr
            dlvo_colloid = self.y
        else:
            raise TypeError("arr_type {} is not valid".format(arr_type))

        dlvo = self.__create_colloid_colloid_array(arr, dlvo_colloid)

        return dlvo

    def __dlvo_interaction_energy(self, arr_type):
        """
        Uses formulation of Israelachvili 1992 Intermolecular surface forces
        to calculate Hamaker constant, followed by the Liang et. al. 2007 to calc
        attractive and repulsive forces

        Parameters:
            arr_type: (str) x direction or y direction , "x", "y"

        Returns:
            dvlo: (np.ndarray) dlvo interaction force from colloids
        """
        kb = 1.31e-23

        if arr_type.lower() == "x":
            c_arr = self.x_distance_array
        elif arr_type.lower() == "y":
            c_arr = self.y_distance_array
        else:
            raise TypeError("arr_type {} is not valid".format(arr_type))

        """
        A = 384. * np.pi * c_arr * kb * self.__params['T']\
            * self.ionic_strength * self.colloid_potential * self.colloid_potential \
            * np.exp(-self.debye * np.abs(c_arr))/ (self.debye * self.debye)
        """
        # use Israelachvili 1991 for hamaker constant
        A = self.colloid_potential * 24 * np.pi * 0.165e-9 ** 2

        lwdv0 = -A / 6.
        lvdw1 = (2. * self.__params['ac'] ** 2.) / (self.__params['ac'] ** 2. + 4. * self.__params['ac'] * c_arr)
        lvdw2 = (2. * self.__params['ac'] ** 2.) / (c_arr + 2. * self.__params['ac']) ** 2.
        lvdw3 = np.log(1. - ((4. * self.__params['ac'] ** 2.) / (c_arr + 2. * self.__params['ac']) ** 2.))

        lewis_vdw = lwdv0 * (lvdw1 + lvdw2 + lvdw3)

        """
        edl0 = 128. * np.pi * self.__params['ac'] * self.__params['ac'] *\
               0.5 * self.ionic_strength * 1.38e-23 * self.__params['T']
        edl1 = (2. * self.__params['ac']) * self.debye ** 2.

        z = 0.
        nz = 0.
        for key, value in self.__params['valence'].items():
            z += float(value)
            nz += 1

        z /= nz # todo: this term may be more correct!

        # z /= 58.44  # todo: look up this term (might be stern length insted!)
        # todo: look into Liang for attractive energy of col-col interaction. Replace for simplification.
        edl2 = np.tanh((z * 1.6e-19 * self.colloid_potential)/(4. * 1.38e-23 * self.__params['T']))
        edl3 = np.exp(-self.debye * c_arr)

        edl = (edl0 / edl1) * (edl2 ** 2.) * edl3
        """

        # original formulation by Derjaguin 1939
        edl0 = (self.__params['epsilon_0'] * self.__params['epsilon_r'] * self.__params['ac'] * self.colloid_potential * self.colloid_potential) / 2.
        edl1 = np.log(1. + np.exp(-self.debye * c_arr))

        edl = edl0 * edl1

        # todo: look more into the dlvo col-col interactions
        dlvo = (edl - lewis_vdw)/c_arr # lewis_vdw + edl)/c_arr

        if arr_type.lower() == "x":
            dlvo[:, :self.__center] *= -1

        elif arr_type.lower() == "y":
            dlvo[self.__center + 1:, :] *= -1

        else:
            raise TypeError("arr_type {} is not valid".format(arr_type))

        dlvo[self.__center, self.__center] = 0.

        return dlvo

    def __angular_array(self, arr_type):
        """
        Calculates the angular proportion of the force a colloid particle
        exerts in grid space, with regard to distance from the colloid.

        Parameters:
            arr_type: (str) delimiter to determine if the array is in the
                x-direction or y-direction

        Return:
            arr (np.ndarray) Array of angular distances adjusted for the proportion
                of force the colloid would be exposed to.
        """

        if 1.01e-6 >= self.__resolution >= 1e-7:
            self.__center = 2
            arr = np.ones((5, 5))
            center = 2

        elif 1e-7 > self.__resolution >= 1e-8:
            self.__center = 25
            arr = np.ones((51, 51))
            center = 25

        elif 1e-8 > self.__resolution >= 1e-9:
            self.__center = 250
            arr = np.ones((501, 501))
            center = 250

        else:
            raise AssertionError("model resolution: {} is out of bounds".format(self.__resolution))

        for i, n in enumerate(arr):
            for j, m in enumerate(n):
                y = float(i - center)
                x = float(j - center)
                if x == 0 and y == 0:
                    arr[i, j] = 0.1

                elif x == 0:
                    arr[i, j] = 1 * np.abs(y)

                elif y == 0:
                    arr[i, j] = 1 * np.abs(x)
                else:
                    arr[i, j] = np.sqrt(x**2 + y**2) + np.abs((m * (np.arctan(y / x) / (np.pi / 2.))))

        if arr_type.lower() == 'x':
            arr = arr.T

        elif arr_type.lower() == 'y':
            pass

        else:
            raise TypeError("arr_type {} is not valid".format(arr_type))

        return arr * self.__resolution  # /1e-6

    def __create_colloid_colloid_array(self, f_arr, c_arr):
        """
        Method to set colloidal forces to a model array.

        Parameters:
        -----------
            f_arr: (np.ndarray) np.zeros array to calculate foces into
            c_arr: (np.ndarray) calculated colloid force array

        Return:
            f_arr: (np.ndarray) an array of colloidal forces in a single primary
                dimension
        """
        center = (c_arr.shape[0] - 1) // 2
        colloids = self.positions

        for colloid in colloids:
            x, y = colloid

            if np.isnan(x) or np.isnan(y):
                pass

            else:
                x -= center
                y -= center

                if x < 0:
                    c_left_x = -x
                    c_right_x = c_arr.shape[1]
                    f_right_x = c_arr.shape[1] + x
                    f_left_x = 0

                elif x + c_arr.shape[1] > f_arr.shape[1]:
                    f_left_x = x
                    f_right_x = f_arr.shape[1]
                    c_left_x = 0
                    c_right_x = -(x - f_arr.shape[1])

                else:
                    c_left_x = 0
                    c_right_x = c_arr.shape[1]
                    f_left_x = x
                    f_right_x = x + c_arr.shape[1]

                if y < 0:
                    c_top_y = -y
                    c_bottom_y = c_arr.shape[0]
                    f_top_y = 0
                    f_bottom_y = c_arr.shape[0] + y

                elif y + c_arr.shape[0] > f_arr.shape[0]:
                    c_top_y = 0
                    c_bottom_y = -(y - f_arr.shape[0])
                    f_top_y = y
                    f_bottom_y = f_arr.shape[0]

                else:
                    c_top_y = 0
                    c_bottom_y = c_arr.shape[0]
                    f_top_y = y
                    f_bottom_y = y + c_arr.shape[0]

                try:
                    f_arr[f_top_y:f_bottom_y, f_left_x:f_right_x] += c_arr[c_top_y:c_bottom_y, c_left_x:c_right_x]
                except ValueError:
                    pass

        return f_arr


# todo: write conversion of force to chemical potential
def force_to_kT(arr, T):
    k = 1.38e-23
    return