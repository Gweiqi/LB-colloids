from LB_Colloid import Colloid
import numpy as np
import sys
import copy


class ForceToVelocity:
    def __init__(self, forces, **kwargs):
        """
        Class that calculates velocity from force

        Inputs:
        -------
        forces: (np.array, np.float) Array of forces felt by a colloid
        ts: (float) Lattice Boltzmann time step value
        rho_colloid: (float) Colloid particle density 
        ac: (float) colloid radius

        Defaults:
        ---------
        rho_colloid: (float) 2650 kg/m^3
        ac: (float) 1e-6

        Returns:
        --------
        velocity (np.array, np.float) Array of velocities calc. from forces
        """
        params = {'rho_colloid': 2650., 'ac': 1e-6, 'ts': 1.}
        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]

        rho_colloid = params['rho_colloid']
        ac = params['ac']
        ts = params['ts']
        self.mass_colloid = (4. / 3.) * np.pi * (ac * ac * ac) * rho_colloid
        self.velocity = (forces * ts) / self.mass_colloid
        
        
class Velocity:
    def __init__(self, LBx, LBy, velocity_factor, **kwargs):
        """
        Class that dimensionalizes LB velocity

        Inputs:
        -------
        LBx: (np.array, np.float) array of Lattice Boltzmann velocities in the x-direction
        LBy: (np.array, np.float) array of Lattice Boltzmann velocities in the y-direction
        ts: (float) time step value, setup is 1. model run should be much less!
        velocity_factor: (float) LB to physical velocity conversion factor

        Returns:
        --------
        xVelocity: (np.array, np.float) array of dimensionalized velocities in the x-direction
        yVelocity: (np.array, np.float) array of dimensionalized velocities in the y-direction
        """
        # todo: add a lb time step for dimensionalization? What is the best way to recover this?
        # todo: Maybe look at reynolds number for dimensionalization

        params = {'lb_timestep': 1e-5, 'ts': 1}

        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]

        # ts = params['lb_timestep']
        ts = params['ts']
        # todo: use the reynolds number calculation and then divide by gridref!
        self.xvelocity = LBx * velocity_factor
        self.yvelocity = LBy * velocity_factor
        

class Gravity:
    def __init__(self, **kwargs):
        """
        Class to estimate the gravitational force experienced by a colloid:

        Inputs:
        -------
        rho_colloid: (float) particle density of a colloid in kg/m*3
        ac: (float) colloid radius in m

        Defaults:
        ---------
        rho_colloid = 2650 kg/m*3 (standard particle density of soil)
        ac = 1e-6 m
        
        Returns:
        --------
        colloid_mass: (float) assumes that colloids are spherical in nature
        gravity: (float) gravitational force that a colloid experiences in vector form
        """

        params = {'rho_colloid': 2650., 'ac': 1e-6}
        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]

        ac = params['ac']
        rho_colloid = params['rho_colloid']
        
        self.colloid_mass = (4./3.)*np.pi*(ac*ac*ac)*rho_colloid
        self.gravity = (self.colloid_mass*-9.81)


class Bouyancy:
    def __init__(self, **kwargs):
        """
        Class to estimate the gravitational force experienced by a colloid:

        Inputs:
        -------
        rho_water: (float) density of water kg/m*3
        rho_colloid: (float) particle density of a colloid in kg/m*3
        ac: (float) colloid radius in m

        Defaults:
        ---------
        rho_colloid = 2650 kg/m*3 (standard particle density of soil)
        ac = 1e-6 m
        
        Returns:
        --------
        colloid_mass: (float) assumes that colloids are spherical in nature
        gravity: (float) gravitational force that a colloid experiences

        Note:
        -----
        acceleration due to gravity is kept positive to maintain the proper vector direction
        """

        params = {'rho_water': 1000., 'rho_colloid': 2650., 'ac': 1e-6}
        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]
        
        rho_water = params['rho_water']
        rho_colloid = params['rho_colloid']
        ac = params['ac']
        
        self.water_mass = (4./3.)*np.pi*(ac*ac*ac)*rho_water
        self.bouyancy = (self.water_mass*rho_water*9.81)/(rho_colloid)


class Brownian:
    def __init__(self, f1, f4, **kwargs):
        """
        Class to estimate brownian forces on colloids

        Inputs:
        -------
        f1: (np.array, np.float) Drag force correction term {Gao et. al. 2010. Computers and Math with App}
        f4: (np.array, np.float) Drag force correction term {Gao et. al. 2010. Computers and Math with App}
        ac: (float) Colloid radius
        viscosity: (float) dynamic viscosity of water
        T = (float) Absolute Temperature in K

        Defaults:
        ---------
        ac: 1e-6 m
        viscosity: 1.002e-3 (dynamic viscosity of water @ 20 C)
        epsion: 6*pi*viscosity*ac {Gao et. al. 2010. Computers and Math with App}
        T = 298.17 K

        Returns:
        --------
        brownian_x: (np.array, np.float) array of browian (random) forces in the x direction {Qiu et. al 2011. VZJ}
        brownian_y: (np.array, np.float) array of browian (random) forces in the x direction {Qiu et. al 2011. VZJ}
        """

        params = {'viscosity': 1.002e-3, 'ac': 1e-6, 'T': 298.17}
        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]
                
        self.ac = params['ac']
        self.viscosity = params['viscosity']
        self.boltzmann = 1.38e-23
        self.epsilon = 6. * np.pi * self.viscosity * self.ac
        self.T = params['T']
        self.diffusive = (self.boltzmann * self.T) / self.epsilon
        self.brownian_x = self.Brown_xforce(self.epsilon, self.diffusive, f4)
        self.brownian_y = self.Brown_yforce(self.epsilon, self.diffusive, f1)

    def Brown_xforce(self, epsilon, diffusivity, f4):
        mu, sigma = 0, 1.
        Fbt = epsilon * np.sqrt(((2 * diffusivity)/(f4 * 1.))) * np.random.normal(mu, sigma, (len(f4), len(f4[0])))
        return Fbt

    def Brown_yforce(self, epsilon, diffusivity, f1):
        mu, sigma = 0, 1.
        Fbn = epsilon * np.sqrt(((2 * diffusivity)/(f1 * 1.))) * np.random.normal(mu, sigma, (len(f1), len(f1[0])))
        return Fbn


class Drag:
    def __init__(self, ux, uy, f1, f2, f3, f4, **kwargs):
        
        """
        Class to calculate colloidal drag forces from fluid velocity arrays
        
        Inputs:
        -------
        ux: (np.array, np.float) fluid velocity in the x-direction
        uy: (np.array, np.float) fluid velocity in the y-direction
        Vx: (np.array, np.float) colloid velocity in the x-direction (assuming equal to ux for calculation)
        Vy: (np.array, np.float) colloid velocity in the y-direction (assuming equal to uy for calculation)
        f1: (np.array, np.float) Drag force correction term {Gao et. al. 2010. Computers and Math with App}
        f2: (np.array, np.float) Drag force correction term {Gao et. al. 2010. Computers and Math with App}
        f3: (np.array, np.float) Drag force correction term {Gao et. al. 2010. Computers and Math with App}
        f4: (np.array, np.float) Drag force correction term {Gao et. al. 2010. Computers and Math with App}
        ac: (float) Colloid radius
        viscosity: (float) dynamic fluid viscosity of water
        rho_colloid: (float) particle density
        rho_water: (float) water density 
        
        Constants:
        ----------
        ac: 1e-6 m
        viscosity: 1/6 (non dimensional viscosity): assumes LB tau == 1
        epsion: 6*pi*viscosity*ac {Gao et. al. 2010. Computers and Math with App}
        viscosity: 1.002e-3 (dynamic viscosity of water @ 20 C)
        rho_colloid: 2650. kg/m**3 (standard particle density of soil)
        rho_water: 1000. kg/m**3 (density of water @ 20C)
        
        Returns:
        --------
        drag_x: (np.array, np.float) vectorized drag forces in the x-direction non-vectorized
        drag_y: (np.array, np.float) vectorized drag forces in the y-direction non-vectorized
        """
        params = {'ac': 1e-6, 'viscosity': 1.002e-3, 'rho_colloid': 2650., 'rho_water': 1000.,
                  'T': 298.17, 'ts': 1.}
        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]        
        
        self.ac = params['ac']
        self.viscosity = params['viscosity']
        self.rho_water = params['rho_water']
        self.rho_colloid = params['rho_colloid']
        self.epsilon = 6. * np.pi * self.viscosity * self.ac
        self.Vcol = -((self.rho_colloid - self.rho_water)*((2*self.ac)**2)*9.81)/(18*self.viscosity)
        self.drag_x = self.drag_xforce(ux, self.Vcol, self.epsilon, f3, f4)  # *xvArr
        self.drag_y = self.drag_yforce(uy, self.Vcol, self.epsilon, f1, f2)  # *yvArr

        self.all_physical_params = copy.copy(params)
        
    def drag_xforce(self, ux, Vx, epsilon, f3, f4):
        Fdt = (epsilon / f4) * ((f3 * ux) - Vx)
        return Fdt

    def drag_yforce(self, uy, Vy, epsilon, f1, f2):
        Fdn = epsilon * ((f2 * uy) - (Vy / f1))
        return Fdn

        
class Gap:
    def __init__(self, xarr, yarr, **kwargs):
        """
        Inputs:
        -------
        xarr: (np.array, np.float) array of x-distances to nearest solid surface
        yarr: (np.array, np.float) array of x-distances to nearest solid surface
        ac: (np.float) radius of a colloid

        Defaults:
        ---------
        ac = 1e-6 m

        Returns:
        -------
        yhbar: (np.array, np.float) Normalized gap distance by colloid radius in y-direction
        xhbar: (np.array, np.float) Normalized gap distance by colloid radius in x-direction
        f1: (np.array, np.float) Drag force correction term {Gao et. al. 2010. Computers and Math with App}
        f2: (np.array, np.float) Drag force correction term {Gao et. al. 2010. Computers and Math with App}
        f3: (np.array, np.float) Drag force correction term {Gao et. al. 2010. Computers and Math with App}
        f4: (np.array, np.float) Drag force correction term {Gao et. al. 2010. Computers and Math with App}

        Note:
        -----
        Passing np.nan can return an overflow warning. 
        """

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
    def __init__(self, xarr, yarr, **kwargs):
        """
        Defaults and parameterization handled by kwargs dictionary feeding into the params dict.
        
        Inputs:
        -------
        xarr: (np.array() float) distances from solid boundaries in the x direction
        yarr: (np.array() float) distances from solid boundaries in the y direction
        valence: (dictionary, int) valences of all species in solution
        concentration: (dictionary, float) optional dictionary of species concentrations for model run Molar
        zeta_colloid: (float) measured_zeta potential of colloid in Volts
        zeta_surface: (float) bulk_zeta potential of porous media in Volts
        adjust_zeta: (bool) boolean flag that adjusts zeta if ionic strength is varied from zeta potential measurement
            ionic strength
        I_initial: (float) required if adjust_zeta is True: Molar ionic strength that zeta potential was measured at
        I: (float) optional [recommended!]: Molar ionic strength of simulated solution 
        ac: (float) colloid radius in meters
        e: (float) electron charge
        epsilon_0: (float) dielectric permativity of a vacuum
        epsilon_r: (float) relative permativity of water
        boltzmann: (float) boltzmann constant
        sheer_plane: (float) equivelent to the thickness of one layer of water molecules
        T: (float) temperature of simulation fluid in degrees Kelvin
        lvdwst_*: (float) * = colloid, water, solid. Lifshits-van der Waals component
        psi+_*: (float) * = colloid, water, solid. Electron acceptor parameter
        psi-_*: (float) * = colloid, water, solid. Electron donor parameter
        xvArr: (np.array() float) np.array of direction vectors. It is inverted in DLVO to represent attractive
            and repulsive forces properly
        yvArr: (np.array() float) np.array of directiion vectors. It is inverted in DLVO to represent attractive
            and repulsive forces properly


        Defaults:
        ---------
        I: 10e-4 Molar
        ac: 1E-6 meters 
        e: 1.6e-19 C (charge of one electron)
        epsilon_0: 8.85e-12 C**2/(J*m) 
        epsilon_r: 78.304 @ 298 K {Malmberg and Maryott 1956. Jour. Res. Nat. Beau. Std. V56(1)}
        valence: {'Na': 1.} [default assumes a Na+ ion solution]
        concentration: {'Na': 10e-4} [default assumes ionic strength of 10e-4]
        boltzmann: 1.38e-23 J/K
        sheer_plane: 3e-10 meters {Interface Science and Technology, 2008. Volume 16 Chapter 3}
        T: 298.17 deg K
        zeta_colloid: -40.5e-3 zeta potential of Na-kaolinite {Chorom 1995. Eur. Jour. of Soil Sci.}
        zeta_solid: -60.9e-3 zeta potential of glass-beads {Ducker 1992. Langmuir V8}
        lvdwst_colloid: 39.9e-3 J/m**2 {Giese et. al. 1996, Jour. Disp. Sci. & Tech. 17(5)}
        lvdwst_solid: 33.7e-3 J/m**2 {Giese et. al. 1996, Jour. Disp. Sci. & Tech. 17(5)}
        lvdwst_water: 21.8e-3 J/m**2 {Interface Science and Technology, 2008. Volume 16. Chapter 2}
        psi+_colloid: 0.4e-3 J/m**2 {Giese et. al. 1996, Jour. Disp. Sci. & Tech. 17(5)}
        psi+_solid: 1.3e-3 J/m**2 {Giese et. al. 1996, Jour. Disp. Sci. & Tech. 17(5)}
        psi+_water: 25.5e-3 J/m**2 {Interface Science and Technology, 2008. Volume 16. Chapter 2}
        psi-_colloid: 34.3e-3 J/m**2 {Giese et. al. 1996, Jour. Disp. Sci. & Tech. 17(5)}
        psi-_solid: 62.2e-3 J/m**2 {Giese et. al. 1996, Jour. Disp. Sci. & Tech. 17(5)}
        psi-_water: 25.5e-3 J/m**2 {Interface Science and Technology, 2008. Volume 16. Chapter 2}

        Output:
        ------
        DLVO forces
        EDLx: (np.array, float) vectorized np.array of electric-double-layer force values in the x-direction
        EDLy: (np.array, float) vectorized np.array of electric-double-layer  force values in the y-direction
        LVDWx: (np.array, float) vectorized np.array of lifshitz-van-der-walls force values in the x-direction
        LVDWy: (np.array, float) vectorized np.array of lifshitz-van-der-walls force values in the y-direction
        LewisABx: (np.array, float) vectorized np.array of lewis acid base force values in the x-direction
        LewisABy: (np.array, float) vectorized np.array of lewis acid base force values in the y-direction
        """

        params = {'concentration': {'Na': 10e-4}, 'adjust_zeta': False, 'I_initial': False, 'I': 10e-4, 'ac': 1e-6,
                  'epsilon_r': 78.304, 'valence': {'Na': 1.}, 'sheer_plane': 3e-10, 'T': 298.17,
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

    def ionic(self, valence, concentration):
        """
        Input is dictionaries of valence and concentration values (Molar)

        Output is float 2*Ionic strength as Z**2 * M
        """
        I = 0
        for key in valence:
            I += (float(concentration[key])*(float(valence[key])**2))
        return I

    def debye(self, epsilon_0, epsilon_r, kb, T, e, ionic_strength):
        NA = 6.02e23
        k_inverse = np.sqrt((epsilon_0*epsilon_r*kb*T)/(e*e*NA*ionic_strength))
        return 1./k_inverse

    def _colloid_potential(self, zeta, ac, kd, z):
        potential = zeta*(1.+(z/ac))*np.exp(kd*z)
        return potential

    def _surface_potential(self, zeta, kd, z):
        potential = zeta*np.exp(kd*z)
        return potential
    
    def _EDL_energy(self, E0, Er, ac, cp, sp, kd, arr):
        """
        Inputs:
        -------
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

    def _Lifshitz_van_der_Walls(self, arr, ac, vdw_st_water, vdw_st_colloid, vdw_st_solid):
        """
        Inputs:
        -------
        arr: (np.array, np.float) array of distances from solid surfaces
        ac: (float) colloid radius
        vdw_st_water: (float) vdW surface tension of water
        vdw_st_colloid: (float) vdW surface tension of colloid
        vdw_st_solid: (float) vdW surface tension (bulk) of solid phase

        constant:
        --------
        h0: contact plane between colloid and surface {Interface Science and Technology, 2008. Volume 16. Chapter 3}

        Output:
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
        Inputs:
        -------
        arr: (np.array, np.float) array of distances from solid surfaces
        e_plus_*: (float) electron acceptor parameter for each specific phase
        e_minus_*: (float) electron donor parameter for each specific phase

        Constants:
        ----------
        h0: contact plane between colloid and surface  {Interface Science and Technology, 2008. Volume 16. Chapter 3}
        chi: water decay length {van Oss 2008}

        Output:
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


# todo: colloid-colloid interaction forces
class ColloidColloid(object):
    """
    Class to include colloid-colloid interaction forces via DLVO chemical
    potential forces.

    Parameters:
    ----------
        arr: (np.ndarray) Any nd.array that represents the shape of the colloid
            domain
        resolution (float) Colloid model resolution
        **kwargs: Please see DLVO class for documentation on kwarg options

    Properties:
    -----------
        x: (np.ndarray) full model colloidal force array for the x direction
        y: (np.ndarray) full model colloidal force array for the y direction
        x_distance_array: (np.ndarray) Calculated angular distance array for
            a colloid in x direction
        y_distance_array: (np.ndarray) Calculated angular distance array for
            a colloid in y direction
        positions: (list) Nx2 list of all colloid positions in model space
        ionic_strength: (float) fluid ionic strength calculated as 2I
        debye: (float) the inverse debye length of the system
        colloid_potential: (float) the calculated colloid surface potential

    Methods:
    -------
        update(colloids):
    """
    def __init__(self, arr, **kwargs):

        self.__params = {'concentration': False, 'adjust_zeta': False, 'I_initial': False,
                         'I': 10e-4, 'ac': 1e-6, 'epsilon_0': 8.85e-12 , 'epsilon_r': 78.304, 'valence': {'Na': 1.},
                         'sheer_plane': 3e-10, 'T': 298.17, 'lvdwst_water': 21.8e-3, 'lvdwst_colloid': 39.9e-3,
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
        self.__resolution = self.__params['lbres']/self.__params['gridref']
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
            colloids: (list, <class: Colloids.LB_Colloid.Colloid)
        """
        self.__reset()
        self.positions

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
        Property method to calculate the inverse debye length on the fly
        """
        if isinstance(self.__debye, bool):
            na = 6.02e23
            k_inverse = np.sqrt((self.__params['epsilon_r']*self.__params['epsilon_r']
                                *self.__params['kb']*self.__params['T'])/
                                (self.__params['e']*self.__params['e']*na*self.ionic_strength))
            self.__debye = 1./k_inverse
        return self.__debye

    @property
    def colloid_potential(self):
        if isinstance(self.__colloid_potential, bool):
            self.__colloid_potential = self.__params['zeta_colloid']*(1. +
                                       (self.__params['sheer_plane']/self.__params['ac']))\
                                        *np.exp(self.debye*self.__params['zeta_colloid'])
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
        if arr_type.lower() == "x":
            c_arr = self.x_distance_array
        elif arr_type.lower() == "y":
            c_arr = self.y_distance_array
        else:
            raise TypeError("arr_type {} is not valid".format(arr_type))

        A = 384. * np.pi * c_arr * self.debye * self.__params['T']\
            * self.ionic_strength * self.colloid_potential * self.colloid_potential\
            * np.exp(-self.debye * np.abs(c_arr))

        lwdv0 = -A / 6.
        lvdw1 = (2. * self.__params['ac'] ** 2.) / (self.__params['ac'] ** 2. + 4. * self.__params['ac'] * c_arr)
        lvdw2 = (2. * self.__params['ac'] ** 2.)/ (c_arr + 2. * self.__params['ac']) ** 2.
        lvdw3 = np.log(1. - ((4. * self.__params['ac'] ** 2.) / (c_arr + 2. * self.__params['ac']) ** 2.))

        lewis_vdw = lwdv0 * (lvdw1 + lvdw2 + lvdw3)

        edl0 = 128. * np.pi * self.__params['ac'] * self.__params['ac'] *\
               self.ionic_strength * 1.38e-23 * self.__params['T']
        edl1 = (2. * self.__params['ac']) * self.debye ** 2.

        z = 0.
        for key, value in self.__params['valence'].items():
            z += float(value)

        z = z / 58.44

        edl2 = np.tanh((z * 1.6e-19 * self.colloid_potential)/(4. * 1.38e-23 * self.__params['T']))
        edl3 = np.exp(-self.debye * c_arr)

        edl = (edl0 / edl1) * (edl2 ** 2.) * edl3

        dlvo = (lewis_vdw + edl)

        if arr_type.lower() == "x":
            dlvo[:, :self.__center] *= -1

        elif arr_type.lower() == "y":
            dlvo[self.__center + 1:, :] *= -1

        else:
            raise TypeError("arr_type {} is not valid".format(arr_type))

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

        if 1e-6 > self.__resolution >= 1e-7:
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
            raise AssertionError("model resolution is out of bounds")


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

        return arr * self.__resolution / 1e-6

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

                f_arr[f_top_y:f_bottom_y, f_left_x:f_right_x] += c_arr[c_top_y:c_bottom_y, c_left_x:c_right_x]

        return f_arr


# todo: write conversion of force to chemical potential
def force_to_kT(arr, T):
    k = 1.38e-23
    return