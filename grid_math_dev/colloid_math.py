import numpy as np
import matplotlib.pyplot as plt
import sys

class Brownian:
    def __init__(self, xarr, yarr, f1, f4, ac=1e-6, viscosity=1./6., T=298.17):
        '''
        Class to estimate brownian forces on colloids

        Inputs:
        -------
        xarr: (np.array, np.float) array of x-distances to nearest solid surface
        yarr: (np.array, np.float) array of x-distances to nearest solid surface
        f1: (np.array, np.float) Drag force correction term {Gao et. al. 2010. Computers and Math with App}
        f4: (np.array, np.float) Drag force correction term {Gao et. al. 2010. Computers and Math with App}
        ac: (float) Colloid radius
        viscosity: (float) macroscopic fluid viscosity term
        T = (float) Absolute Temperature in K

        Defaults:
        ---------
        ac: 1e-6 m
        viscosity: 1/6 (non dimensional viscosity): assumes LB tau == 1
        epsion: 6*pi*viscosity*ac {Gao et. al. 2010. Computers and Math with App}
        T = 298.17 K

        Returns:
        --------
        brownian_x: (np.array, np.float) array of browian (random) forces in the x direction {Qiu et. al 2011. VZJ}
        brownian_y: (np.array, np.float) array of browian (random) forces in the x direction {Qiu et. al 2011. VZJ}
        '''
        self.ac = ac
        self.viscosity = viscosity
        self.boltzmann = 1.38e-23
        self.epsilon = 6.* np.pi * self.viscosity * self.ac
        self.T = T
        self.diffusive = (self.boltzmann * self.T) / self.epsilon
        self.brownian_x = self.Brown_xforce(self.epsilon, self.diffusive, f4)
        self.brownian_y = self.Brown_yforce(self.epsilon, self.diffusive, f1)

    def Brown_xforce(self, epsilon, diffusivity, f4):
        mu, sigma = 0, 1.
        Fbt = epsilon * np.sqrt(((2 * diffusivity)/(f4 * 1.))) * np.random.normal(mu, sigma, (len(f4),len(f4[0])))
        return Fbt

    def Brown_yforce(self, epsilon, diffusivity, f1):
        mu, sigma = 0, 1.
        Fbn = epsilon * np.sqrt(((2 * diffusivity)/(f1 * 1.))) * np.random.normal(mu, sigma, (len(f1),len(f1[0])))
        return Fbn


class Drag:
    def __init__(self, ux, uy, Vx, Vy, f1, f2, f3, f4, ac=1e-6, viscosity=1./6.):
        
        '''
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
        viscosity: (float) macroscopic fluid viscosity term
        
        Constants:
        ----------
        ac: 1e-6 m
        viscosity: 1/6 (non dimensional viscosity): assumes LB tau == 1
        epsion: 6*pi*viscosity*ac {Gao et. al. 2010. Computers and Math with App}

        Returns:
        --------
        drag_x: (np.array, np.float) drag forces in the x-direction
        drag_y: (np.array, np.float) drag forces in the y-direction
        '''
        self.ac = ac
        self.viscosity = viscosity
        self.epsilon = 6. * np.pi * self.viscosity * self.ac
        self.drag_x = self.drag_xforce(ux, Vx, self.epsilon, f3, f4)
        self.drag_y = self.drag_yforce(uy, Vy, self.epsilon, f1, f2)
                                       
    def drag_xforce(self, ux, Vx, epsilon, f3, f4):
        Fdt = (epsilon / f4) * ((f3 * ux) - Vx)
        return Fdt

    def drag_yforce(self, uy, Vy, epsilon, f1, f2):
        Fdn = epsilon * ((f2 * uy) - (Vy / f1))
        return Fdn

        
class Gap:
    def __init__(self, xarr, yarr, ac=1e-6):
        ### add docstrings and **kwargs support
        '''
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
        Passing np.nan to these returns an overflow warning. 
        '''
        self.ac = ac
        self.yhbar = np.abs(yarr/self.ac)
        self.xhbar = np.abs(xarr/self.ac)
        self.f1 = self.set_f1(self.yhbar)
        self.f2 = self.set_f2(self.yhbar)
        self.f3 = self.set_f3(self.xhbar)
        self.f4 = self.set_f4(self.xhbar)
        
    def set_f1(self, yhbar):
         f1 = 1.0 - 0.443 * np.exp(yhbar * -1.299) - 0.5568 * np.exp((yhbar ** 0.75) * -0.32)
         print(f1)
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
        '''
        Defaults and parameterization handled by kwargs dictionary feeding into the params dict.
        
        Inputs:
        -------
        xarr: (np.array() float) distances from solid boundaries in the x direction
        yarr: (np.array() float) distances from solid boundaries in the y direction
        valence: (dictionary, int) valences of all species in solution
        concentration: (dictionary, float) optional dictionary of species concentrations for model run Molar
        zeta_colloid: (float) measured_zeta potential of colloid in Volts
        zeta_surface: (float) bulk_zeta potential of porous media in Volts
        adjust_zeta: (bool) boolean flag that adjusts zeta if ionic strength is varied from zeta potential measurement ionic strength
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
        
        
        '''

        params = {'concentration': {'Na': 10e-4}, 'adjust_zeta': False, 'I_initial': None, 'I': 10e-4, 'ac': 1e-6,
                  'e': 1.6e-19, 'epsilon_0': 8.85e-12, 'epsilon_r': 78.304, 'valence': {'Na': 1.}, 'boltzmann': 1.38e-23,
                  'sheer_plane': 3e-10, 'T': 298.17, 'lvdwst_water': 21.8e-3, 'lvdwst_colloid': 39.9e-3,
                  'lvdwst_solid': 33.7e-3, 'zeta_colloid': -40.5e-3, 'zeta_solid': -60.9e-3, 'psi+_colloid': 0.4e-3,
                  'psi-_colloid': 34.3e-3, 'psi+_water': 25.5e-3, 'psi-_water': 25.5e-3, 'psi+_solid': 1.3e-3,
                  'psi-_solid': 62.2e-3}

        for kwarg in kwargs:
            params[kwarg] = kwargs[kwarg]
                  
        self.epsilon_0 = params['epsilon_0']
        self.epsilon_r = params['epsilon_r']
        self.ac = params['ac']
        self.e = params['e']
        self.valence = params['valence']
        self.concentration = params['concentration'] 
        self.boltzmann = params['boltzmann']
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

        if params['adjust_zeta'] is not False:
            if parmas['I_initial'] is not None:
                self.k_debye_init(self.epsilon_0, self.epsilon_r, self.boltzmann, self.T, self.e, params['I_initial'])
                
                self.colloid_potential_init = self._colloid_potential(self.zeta_colloid, self.ac, self.k_debye,
                                                                      self.stern_z)
                self.surface_potential_init = self._surface_potential(self.zeta_solid, self.k_debye, self.stern_z)

                if params['I'] is not None:
                    self.ionic_strength = 2*params['I'] #2I is what is used in the debye equation
                else:
                    self.ionic_strength = self.ionic(parmas['valence'], params['concentration'])

                self.k_debye = self.debye(self.epsilon_0, self.epsilon_r, self.boltzmann, self.T, self.e,
                                          self.ionic_strength)

                self.zeta_colloid = self._adjust_zeta_colloid(self.colloid_potential_init, self.ac, self.k_debye,
                                                              self.stern_z)

                self.zeta_solid = self._adjust_zeta_surface(self.surface_potential_init, self.k_debye,
                                                              self.stern_z)
            else:
                print('Please suppy initial ionic strength')
                sys.exit(-1)
                

        else: 
            pass
        
        if params['I'] is not None:
            self.ionic_strength = 2*params['I'] #2I is what is used in the debye equation
        else:
            self.ionic_strength = self.ionic(params['valence'], params['concentration'])
        
        self.k_debye = self.debye(self.epsilon_0, self.epsilon_r, self.boltzmann, self.T, self.e,
                                  self.ionic_strength)
                                         
        self.colloid_potential = self._colloid_potential(self.zeta_colloid, self.ac, self.k_debye, self.stern_z)
        self.surface_potential = self._surface_potential(self.zeta_solid, self.k_debye, self.stern_z)

        # returns values in energy. need forces, so divide by array
        
        self.EDLx = self._EDL_energy(self.epsilon_0, self.epsilon_r, self.ac, self.colloid_potential,
                                     self.surface_potential, self.k_debye, xarr)/xarr

        self.EDLy = self._EDL_energy(self.epsilon_0, self.epsilon_r, self.ac, self.colloid_potential,
                                     self.surface_potential, self.k_debye, yarr)/yarr
        # EDL forces: TODO: we need to re-add direction vectors after calculation
        
        self.LVDWx = self._Lifshitz_van_der_Walls(xarr, self.ac, self.lvdwst_water, self.lvdwst_colloid,
                                                  self.lvdwst_solid)/xarr

        self.LVDWy = self._Lifshitz_van_der_Walls(xarr, self.ac, self.lvdwst_water, self.lvdwst_colloid,
                                                  self.lvdwst_solid)/yarr

        self.LewisABx = self._lewis_acid_base(xarr, self.ac, self.eplus_colloid, self.eplus_solid, self.eplus_water,
                                              self.eneg_colloid, self.eneg_solid, self.eneg_water)/xarr

        self.LewisABy = self._lewis_acid_base(yarr, self.ac, self.eplus_colloid, self.eplus_solid, self.eplus_water,
                                              self.eneg_colloid, self.eneg_solid, self.eneg_water)/yarr
        

        
    def ionic(self, valence, concentration):
        '''
        Input is dictionaries of valence and concentration values (Molar)

        Output is float 2*Ionic strength as Z**2 * M
        '''
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
        '''
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
        '''

        EDL0 = np.pi*E0*Er*ac
        EDL1 = 2.*sp*cp
        EDL2 = np.log((1. + np.exp(-kd*np.abs(arr)))/(1. - np.exp(-kd*np.abs(arr))))
        EDL3 = sp*sp + cp*cp
        EDL4 = np.log(1. - np.exp(-2.*kd*np.abs(arr)))

        EDL = EDL0*(EDL1*EDL2 + EDL3*EDL4)

        return EDL

    def _adjust_zeta_colloid(self, potential, ac, kd, z):
        zeta = potential/((1. + (z/ac))*np.exp(kd*z))
        return zeta

    def _adjust_zeta_surface(self, potential, kd, z):
        zeta = potential/(np.exp(kd*z))
        return zeta

    def _Lifshitz_van_der_Walls(self, arr, ac, vdw_st_water, vdw_st_colloid, vdw_st_solid):
        '''
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
        LVDW: (np.array, np.float) array of lifshitz_vdW interaction energies
        '''
        
        h0= 1.57e-10 

        LVDW0 = -4.*np.pi*(h0/arr)
        LVDW1 = np.sqrt(vdw_st_water) - np.sqrt(vdw_st_solid)
        LVDW2 = np.sqrt(vdw_st_water) - np.sqrt(vdw_st_colloid)

        LVDW = LVDW0*LVDW1*LVDW2
        return LVDW

    def _lewis_acid_base(self, arr, ac, eplus_colloid, eplus_solid, eplus_water, eneg_colloid,
                         eneg_solid, eneg_water):
        '''
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
        LAB: (np.array, np.float) array of lewis acid base interaction energies
        
        '''
        h0 = 1.57e-10
        chi = 0.6e-10

        LAB0 =  -4.*np.pi*h0*ac
        LAB1 = np.exp((h0-arr)/chi)
        LAB2 = np.sqrt(eplus_water)*(np.sqrt(eneg_colloid) + np.sqrt(eneg_solid) - np.sqrt(eneg_water))
        LAB3 = np.sqrt(eneg_water)*(np.sqrt(eplus_colloid) + np.sqrt(eplus_solid) - np.sqrt(eplus_water))
        LAB4 = np.sqrt(eplus_colloid*eneg_solid)
        LAB5 = np.sqrt(eneg_colloid*eplus_solid)

        LAB = LAB0*LAB1*(LAB2+LAB3-LAB4-LAB5)
        return LAB
        
