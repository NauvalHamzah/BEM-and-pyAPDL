# -*- coding: utf-8 -*-

"""
Module for storing rotor properties and calculation of induction factors and forces for the airfoil sections.
"""
import pandas as pd
import numpy as np
from configparser import NoOptionError
from math import radians, degrees, sqrt, cos, sin, atan2, atan, pi, acos, exp
from .airfoil import load_airfoil

class Rotor: 
    """
    Holds rotor properties and a list of all airfoil sections.

    :param configparser.SafeConfigParser cfg: Configuration object
    :param string name: Name of rotor
    :param string mode: Solver mode
    """
    def __init__(self, cfg, name, mode, input_name):
        self.n_blades = cfg.getint(name, 'nblades')
        self.diameter = cfg.getfloat(name, 'diameter')
        self.radius_hub = cfg.getfloat(name,'radius_hub')
        s,c,r,self.theta=self.SecToTip(cfg, name, mode, input_name)

        n_elements = len(s)

        self.sections = []
        for i in range(n_elements): 
            sec = Section(load_airfoil(s[i]), float(r[i]), radians(self.theta[i]), float(c[i]), self, mode)
            self.sections.append(sec)
        

        self.precalc(twist=0.0)
    
    def SecToTip(self, cfg, name, mode, input_name):
        s_input = cfg.get(name, 'section').split()
        c_input = [float(p) for p in cfg.get(name, 'chord').split()]
        r_input = [float(p) for p in cfg.get(name, 'radius').split()]
        theta_input = [float(p) for p in cfg.get(name, 'pitch').split()]
        n_sections = len(s_input)
        n_elements = n_sections-1

        blade_name=input_name
        f = open("output/"+blade_name+"-Airfoil.txt", "w+")
        for i in range (len(s_input)):    
            f.write(s_input[i]+"\n")
        f.close()

        s=s_input
        c=c_input
        r=r_input
        theta=theta_input
        
        return s,c,r,theta


    def precalc(self, twist):
        """
        Calculation of properties before each solver run, to ensure all parameters are correct for parameter sweeps.

        :return: None
        """
        self.blade_radius = 0.5*self.diameter
        self.area = pi*self.blade_radius**2

        # Apply twist
        for i,sec in enumerate(self.sections):
            sec.pitch = radians(self.theta[i] + twist)


    def sections_dataframe(self):
        """
        Creates a pandas DataFrame with all calculated section properties.

        :return: DataFrame with section properties
        :rtype: pd.DataFrame
        """

        columns = ['radius','chord','pitch','Cl','Cd','Cl_Cd','Cn','Ct','Pn','Pt','dT','dQ','F','a','ap','Re','v','v_theta','v_rel','alpha','phi']
        data = {}
        for param in columns:
            array = [getattr(sec, param) for sec in self.sections]
            data[param] = array
            
        
        return pd.DataFrame(data)
    

 

class Section: 
    """
    Class for calculating induction factors and forces according to the BEM theory for a single airfoil section.

    :param Airfoil airfoil: Airfoil of the section
    :param float radius: Distance from center to mid of section
    :param float width: Width of section
    :param float pitch: Pitch angle in radians
    :param float chord: Chord length of section
    :param Rotor rotor: Rotor that section belongs to
    :param string mode: Solver mode
    """
    def __init__(self, airfoil, radius, pitch, chord, rotor, mode):
        self.airfoil = airfoil
        self.radius = radius
        self.width = 0
        self.pitch = pitch
        self.chord = chord
        self.rotor = rotor
        
        if mode == 'turbine':
            self.C = -1
        else:
            self.C = 1
        
        self.v = 0.0
        self.v_theta = 0.0
        self.v_rel = 0.0
        self.a=0.0
        self.ap=0.0
        self.Re = 0.0
        self.alpha = 0.0
        self.Pn = 0.0
        self.Pt = 0.0
        self.F = 0.0
        self.Cl = 0.0
        self.Cd = 0.0
        self.Cl_Cd = 0.0
        self.Cn = 0.0
        self.Ct = 0.0
        self.phi = 0.0
        self.F_t = 0.0
        self.F_q = 0.0
        self.sigma = self.rotor.n_blades*self.chord/(2*pi*self.radius)
        self.dT = 0
        self.dQ = 0

    def iterate(self,v_inf,omega):
        self.sigma
        gamma = (omega*self.radius)/v_inf
        a_init = 0.25*(2+(pi*gamma*self.sigma)-(sqrt(4-4*pi*gamma*self.sigma+pi*gamma**2*self.sigma*(8*self.pitch+pi*self.sigma))))
        ap_init = 0
        #a_init=0.3
        #ap_init=0.3
        n_iter = 0
        conv = 1e-5
        error = 1
        n=0
        a = a_init
        ap = ap_init
        while n<10 and error>conv:
            v = (1-a)*v_inf
            vp = (1+ap)*omega*self.radius
            phi = degrees(np.arctan(v/vp))
            check = np.isnan([phi])

            if check==True:
                a=0.3
                ap=0.3
                n_iter=0
                v = (1-a)*v_inf
                vp = (1+ap)*omega*self.radius
                phi = degrees(np.arctan(v/vp))

            alpha = phi-degrees(self.pitch)
            Cl = self.airfoil.Cl(radians(alpha))
            Cd = self.airfoil.Cd(radians(alpha))
            Cn = Cl*cos(radians(phi)) + Cd*sin(radians(phi))
            Ct = Cl*sin(radians(phi)) - Cd*cos(radians(phi))
            F = self.tip_loss(radians(phi))
            K = (4*F*(sin(radians(phi)))**2)/(self.sigma*Cn)
            a_temp = 1/(K+1)

            if a_temp<0.2:
                a_temp = 1/(K+1)
                a_new = a_temp
            else:
                bracket = ((K*(1-0.4)+2)**2)+4*(K*0.2**2-1)
                a_corr = 0.5*(2+K*(1-0.4)-sqrt(bracket))
                a_new = a_corr

            Kp = 4*F*(sin(radians(phi)))*(cos(radians(phi)))/(self.sigma*Ct)
            ap_new = 1/(Kp-1)

            a_error = abs(a_new-a)
            ap_error = abs(ap_new-ap)
            error = max(a_error,ap_error)
            n_iter = n_iter + 1
            n=n+1

            a = a_new
            ap= ap_new
            
        return radians(phi)
        

    def tip_loss(self, phi):
        """
        Prandtl tip loss factor, defined as

        .. math::
            F = \\frac{2}{\\pi}\\cos^{-1}e^{-f} \\\\
            f = \\frac{B}{2}\\frac{R-r}{r\\sin\\phi}

        A hub loss is also caluclated in the same manner.

        :param float phi: Inflow angle
        :return: Combined tip and hub loss factor
        :rtype: float
        """
        def prandtl(dr, r, phi):
            f = self.rotor.n_blades*dr/(2*r*(sin(phi)))
            if (-f > 500): # exp can overflow for very large numbers
                F = 1.0
            else:
                F = 2*acos(min(1.0, exp(-f)))/pi
                
            return F
        
        if phi == 0:
            F = 1.0
        else:    
            r = self.radius
            Ftip = prandtl(self.rotor.blade_radius - r, r, phi)
            Fhub = prandtl(r - self.rotor.radius_hub, r, phi)
            F = Ftip*Fhub
        #print("\t\t7. Masuk tip-loss")    
        self.F = F
        return F
 
                    
    def airfoil_forces(self, phi):
        """
        Force coefficients on an airfoil, decomposed in axial and tangential directions:

        .. math::
            C_T = C_l\\cos{\\phi} - CC_d\\sin{\\phi} \\\\
            C_Q = C_l\\sin{\\phi} + CC_d\\cos{\\phi} \\\\

        where drag and lift coefficients come from
        airfoil tables.

        :param float phi: Inflow angle
        :return: Axial and tangential force coefficients
        :rtype: tuple
        """
        #print("\t\t\t8. Masuk airfoil_forces")
        C = self.C

        alpha = C*(self.pitch - phi)
                
        Cl = self.airfoil.Cl(alpha)
        Cd = self.airfoil.Cd(alpha)
                
        Cn = Cl*cos(phi) - C*Cd*sin(phi)
        Ct = Cl*sin(phi) + C*Cd*cos(phi)

        self.Cl = Cl
        self.Cd = Cd
        self.Cl_Cd = Cl/Cd
        self.Cn = Cn
        self.Ct = Ct
        self.alpha = alpha
        
        return Cn, Ct
    
    def induction_factors(self, phi):
        """
        Calculation of axial and tangential induction factors,

        .. math::
            a = \\frac{1}{\\kappa - C} \\\\
            a\' = \\frac{1}{\\kappa\' + C} \\\\
            \\kappa = \\frac{4F\\sin^2{\\phi}}{\\sigma C_T} \\\\
            \\kappa\' = \\frac{4F\\sin{\\phi}\\cos{\\phi}}{\\sigma C_Q} \\\\
            
        :param float phi: Inflow angle
        :return: Axial and tangential induction factors
        :rtype: tuple
        """
        #print("\t6. Masuk induction factor")
        C = self.C
        
        F = self.tip_loss(phi)
        
        Cn, Ct = self.airfoil_forces(phi)
        
        kappa = 4*F*sin(phi)**2/(self.sigma*Cn)
        kappap = 4*F*sin(phi)*cos(phi)/(self.sigma*Ct)

        a = 1.0/(kappa - C)
        ap = 1.0/(kappap + C)

        
        #Spera Correction
        bracket = ((kappa*(1-0.4)+2)**2)+4*(kappa*0.2**2-1)
        if a>0.2:
            if bracket > 0:
                a = 0.5*(2+kappa*(1-0.4)-sqrt(bracket))
            else:
                a=a
        '''
        #Modified Glauert Correction
        lim = 1+((1-a)**2*self.sigma*Ct/(sin(phi)**2))
        bracket = lim*(50-36*F)+12*F*(3*F-4)
        if a>0.2:
            if bracket>0:
                a = (18*F-20-3*sqrt(bracket)/(36*F-50))
            else:
                a = a
        '''
        
        return a, ap
        
    def func(self, phi, v_inf, omega):
        """
        Residual function used in root-finding functions to find the inflow angle for the current section.

        .. math::
            \\frac{\\sin\\phi}{1+Ca} - \\frac{V_\\infty\\cos\\phi}{\\Omega R (1 - Ca\')} = 0\\\\

        :param float phi: Estimated inflow angle
        :param float v_inf: Axial inflow velocity
        :param float omega: Tangential rotational velocity
        :return: Residual
        :rtype: float
        """
        # Function to solve for a single blade element
        C = self.C
        #print("phi : "+ str(degrees(phi)))
        a, ap = self.induction_factors(phi)
        
        resid = sin(phi)/(1 + C*a) - v_inf*cos(phi)/(omega*self.radius*(1 - C*ap))

        self.a = a
        self.ap = ap
        self.phi = phi

        return resid
    
    def forces(self, phi, v_inf, omega, fluid):
        """
        Calculation of axial and tangential forces (thrust and torque) on airfoil section. 

        The definition of blade element theory is used,

        .. math::
            \\Delta T = \\sigma\\pi\\rho U^2C_T r\\Delta r \\\\
            \\Delta Q = \\sigma\\pi\\rho U^2C_Q r^2\\Delta r \\\\
            U = \\sqrt{v^2+v\'^2} \\\\
            v = (1 + Ca)V_\\infty \\\\
            v\' = (1 - Ca\')\\Omega R \\\\

        Note that this is equivalent to the momentum theory definition,

        .. math::
            \\Delta T = 4\\pi\\rho r V_\\infty^2(1 + Ca)aF\\Delta r \\\\
            \\Delta Q = 4\\pi\\rho r^3 V_\\infty\\Omega(1 + Ca)a\'F\\Delta r \\\\


        :param float phi: Inflow angle
        :param float v_inf: Axial inflow velocity
        :param float omega: Tangential rotational velocity
        :param Fluid fluid: Fluid 
        :return: Axial and tangential forces
        :rtype: tuple
        """
        #print("5. Masuk forces")
        C = self.C
        r = self.radius
        rho = fluid.rho
        
        a, ap = self.induction_factors(phi)
        #print(a)
        Cn, Ct = self.airfoil_forces(phi)
        
        v = (1 + C*a)*v_inf
        vp = (1 - C*ap)*omega*r
        U = sqrt(v**2 + vp**2)   
        
        self.a = a
        self.ap = ap
        self.phi = phi
        self.Re = rho*U*self.chord/fluid.mu
        self.v = v
        self.v_theta = vp
        self.v_rel = U

        
        # From blade element theory
        self.Pn = 0.5*rho*U**2*self.chord*Cn
        self.Pt = 0.5*rho*U**2*self.chord*Ct
             
        return self.Pn, self.Pt

    def integrate (self, Pn, Pt, i):
        columns=['radius', 'Pn', 'Pt']
        data = {}
        for param in columns:
            array = [getattr(sec, param) for sec in self.rotor.sections]
            data[param] = array
        
        d=data.values()
        d_list=list(d)

        dr = (d_list[0][i]-d_list[0][i-1])
        r = ((d_list[0][i]+d_list[0][i-1])/2)

        self.dT = self.rotor.n_blades*dr*((d_list[1][i]+d_list[1][i-1])/2)
        self.dQ = self.rotor.n_blades*dr*r*((d_list[2][i]+d_list[2][i-1])/2)

        return self.dT, self.dQ
    

        
    

