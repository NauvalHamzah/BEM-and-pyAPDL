# -*- coding: utf-8 -*-
"""
Module for the solver class.
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import optimize
from configparser import SafeConfigParser
from math import radians, degrees, sqrt, pi, floor
from .fluid import Fluid
from .rotor import Rotor
from .rotor import Section


class Solver: 
    """
    The Solver object loads the config file and contains functions for running a single simulation,
    parameter sweeps and optimization.

    :param string config_path: Path to config file
    """
    def __init__(self, config_path, input_name):

        # Read configuration file
        cfg = SafeConfigParser()
        cfg.read(config_path)
        
        # Case
        self.v_inf = cfg.getfloat('case', 'v_inf')
        self.tsr = cfg.getfloat('case', 'tsr')
        

        if cfg.has_option('case', 'twist'):
            self.twist = cfg.getfloat('case', 'twist')
        else:
            self.twist = 0.0
        if cfg.has_option('case', 'coaxial'):
            self.coaxial = cfg.getboolean('case', 'coaxial')
        else:
            self.coaxial = False
        
        # Rotor
        if cfg.has_section('turbine'):
            self.mode = 'turbine'
            self.rotor = Rotor(cfg, 'turbine', self.mode, input_name)
        else:
            self.mode = 'rotor'
            self.rotor = Rotor(cfg, 'rotor', self.mode)
        
        # Fluid
        self.fluid = Fluid(cfg)
        
        # Output
        self.T = 0 # Thrust
        self.Q = 0 # Torque
        self.P = 0 # Power
        
        # Coaxial
        if self.coaxial:
            self.rpm2 = cfg.getfloat('case','rpm2')
            if cfg.has_option('case', 'twist2'):
                self.twist2 = cfg.getfloat('case', 'twist2')
            else:
                self.twist2 = 0.0
            self.rotor2 = Rotor(cfg, 'rotor2', self.mode)
            self.zD = cfg.getfloat('case','dz')/self.rotor.diameter
            self.T2 = 0
            self.Q2 = 0
            self.P2 = 0

        # Solver
        self.solver = 'bisect'
        self.Cs = 0.625
        if cfg.has_section('solver'):
            self.solver = cfg.get('solver','solver')
            if cfg.has_option('solver', 'Cs'):
                self.Cs = cfg.getfloat('solver','Cs')
       
    def rotor_coeffs(self, T, Q, P):
        """ 
        Dimensionless coefficients for a rotor. 

        .. math::
            \\text{J} = \\frac{V_\\infty}{nD} \\\\
            C_T = \\frac{T}{\\rho n^2 D^4} \\\\
            C_Q = \\frac{Q}{\\rho n^2 D^5} \\\\
            C_P = 2\\pi C_Q \\\\
            \\eta = \\frac{C_T}{C_P}J \\\\

        :param float T: Thrust
        :param float Q: Torque
        :param float P: Power
        :return: Advance ratio, thrust coefficient, torque coefficient, power coefficient and efficiency
        :rtype: tuple
        """

        D = self.rotor.diameter
        R = 0.5*D
        rho = self.fluid.rho
        self.rpm = (self.tsr*self.v_inf/R)*(60/(2*np.pi))
        n = self.rpm/60.0
        J = self.v_inf/(n*D)
        omega = self.rpm*2*pi/60.0
 
        CT = T/(rho*n**2*D**4)
        CQ = Q/(rho*n**2*D**5)
        CP = 2*pi*CQ

        if J==0.0:
            eta = (CT/CP)
        else:
            eta = (CT/CP)*J

        return J, CT, CQ, CP, eta

    def turbine_coeffs(self, T, Q, P):
        """
        Dimensionless coefficients for a turbine.

        .. math::
            \\text{TSR} = \\frac{\\Omega R}{V_\\infty} \\\\
            C_T = \\frac{2T}{\\rho A V_\\infty^2} \\\\
            C_P = \\frac{2P}{\\rho A V_\\infty^3} \\\\

        :param float T: Thrust
        :param float Q: Torque
        :param float P: Power
        :return: Tip-speed ratio, power coefficient and thrust coefficient
        :rtype: tuple
        """
        #print("11. Turbine coeff")
        rho = self.fluid.rho
        V = self.v_inf
        omega = self.rpm*2*pi/60.0
        rpm = omega*60/(2*np.pi)
        TSR = omega*self.rotor.blade_radius/V
        CT = T/(0.5*rho*self.rotor.area*V**2)
        CP = P/(0.5*rho*self.rotor.area*V**3)

        return rpm, CP, CT

        
    def run_sweep(self, parameter, n, low, high):
        """
        Utility function to run a sweep of a single parameter.

        :param string parameter: Parameter to sweep, must be a member of the Solver class.
        :param int n: Number of runs
        :param float low: Minimum parameter value
        :param float high: Maximum parameter value

        :return: DataFrame of results and list of sections for each run
        :rtype: tuple
        """

        if self.mode == 'turbine':
            df = pd.DataFrame(columns = [parameter, 'T', 'Q', 'P', 'RPM', 'CT', 'CP'], index=range(n))
        else:
            if self.coaxial:
                cols = [parameter, 'T', 'Q', 'P', 'J', 'CT', 'CQ', 'CP', 'eta', 
                        'CT2', 'CQ2', 'CP2', 'eta2']
            else:
                cols = [parameter, 'T', 'Q', 'P', 'J', 'CT', 'CQ', 'CP', 'eta']

            df = pd.DataFrame(columns = cols, index=range(n))
        #print("1. Masuk run_sweep")
        sections = []
        trial = []
        for i,p in enumerate((np.linspace(low, high, n))):
            setattr(self, parameter, p)
            if self.mode == 'turbine':
                T,Q,P,sec_df = self.run()
                rpm,CP,CT = self.turbine_coeffs(T, Q, P)
                df.iloc[i] = [p, T, Q, P, rpm, CT, CP]
                
            else:
                if self.coaxial:
                    T,Q,P,sec_df,T2,Q2,P2,sec_df2 = self.run()
                    J,CT,CQ,CP,eta = self.rotor_coeffs(T, Q, P)
                    J,CT2,CQ2,CP2,eta = self.rotor_coeffs(T2, Q2, P2)
                    df.iloc[i] = [p, T, Q, P, T2, Q2, P2, J, CT, CQ, CP, eta, CT2, CP2, eta2]
                else:
                    T,Q,P,sec_df = self.run()
                    J,CT,CQ,CP,eta = self.rotor_coeffs(T, Q, P)
                    df.iloc[i] = [p, T, Q, P, J, CT, CQ, CP, eta]

            sections.append(sec_df)          
           
        
        return df, sections
    
    def solve(self, rotor, twist, tsr, v_inflow, D):
        """
        Find inflow angle and calculate forces for a single rotor given rotational speed, inflow velocity and radius.

        :param Rotor rotor: Rotor to solve for
        :param float twist: Angle to adjust rotor pitch
        :param float tsr: Tip Speed Ratio
        :param float v_inflow: Inflow velocity
        :param float r_inflow: Inflow radius (equal to blade radius for single rotors)
        :return: Calculated thrust, torque and power for the rotor
        :rtype: tuple
        """
        #print("3. Masuk solve")
        rotor.precalc(twist)
        
        rpm = (tsr*v_inflow/(D/2))*(60/(2*np.pi))
        omega = rpm*2*pi/60.0
        # Axial momentum (thrust)
        T = 0.0
        # Angular momentum
        Q = 0.0
        elem=0
        for sec in rotor.sections:
            v_inf = v_inflow
            #phi = sec.iterate(v_inf,omega)
            
            if self.solver == 'brute':
                phi = self.brute_solve(sec, v_inf, omega)
            else:
                try:
                    phi = optimize.bisect(sec.func, 0.001*pi, 0.99*pi, args=(v_inf, omega))
                except ValueError as e:
                    print(e)
                    print('Bisect failed, switching to brute solver')
                    phi = self.brute_solve(sec, v_inf, omega)
            
            Pn, Pt = sec.forces(phi, v_inf, omega, self.fluid)
            
            if elem == 0:
                elem=elem+1
                continue
            else:
                dT, dQ = sec.integrate(Pn, Pt, elem)
            
            elem=elem+1
            
            T = T+dT
            Q = Q+dQ
        # Power
        P = Q*omega  

        return T, Q, P

    def slipstream(self):
        """
        For coaxial calculations. Calculates slipstream radius and velocity for the upper rotor according to
        momentum theory. Currently only the static case is included.

        .. math::
            r_s = \\frac{R}{\\sqrt{2}} \\\\
            v_s = C_s\\sqrt{\\frac{2 T}{\\rho A}} \\\\

        :return: Radius and velocity of the slipstream
        :rtype: tuple
        """

        r_s = self.rotor.blade_radius/sqrt(2.0)
        v_s = self.Cs*sqrt(2*self.T/(self.fluid.rho*self.rotor.area))

        return r_s, v_s

 
    def run(self):
        """
        Runs the solver, i.e. finds the forces for each rotor.

        :return: Calculated thrust, torque, power and DataFrame with properties for all sections.
        :rtype: tuple
        """
        #print("2. Masuk run()")
        self.rpm = (self.tsr*self.v_inf/(self.rotor.diameter/2))*(60/(2*np.pi))
        self.T, self.Q, self.P = self.solve(self.rotor, self.twist, self.tsr, self.v_inf, self.rotor.diameter)
       
        print('--- Results ---')
        print('rpm:\t',self.rpm)
        print('Trust (N):\t',self.T)
        print('Torque (Nm):\t',self.Q)
        print('Power (W):\t',self.P)

        # Coaxial calculaction
        if self.coaxial:
            self.r_s, self.v_s = self.slipstream()
           
            self.T2, self.Q2, self.P2 = self.solve(self.rotor2, self.twist2, self.rpm2, self.v_s, self.r_s)

            print('Trust 2 (N):\t',self.T2)
            print('Torque 2 (Nm):\t',self.Q2)
            print('Power 2 (W):\t',self.P2)

            return self.T, self.Q, self.P, self.rotor.sections_dataframe(), self.T2, self.Q2, self.P2, self.rotor2.sections_dataframe()
        
        else:
            return self.T, self.Q, self.P, self.rotor.sections_dataframe()


    def brute_solve(self, sec, v, omega, n=3600):
        """ 
        Solve by a simple brute force procedure, iterating through all
        possible angles and selecting the one with lowest residual.

        :param Section sec: Section to solve for
        :param float v: Axial inflow velocity
        :param float omega: Tangential rotational velocity
        :param int n: Number of angles to test for, optional
        :return: Inflow angle with lowest residual
        :rtype: float
        """
        resid = np.zeros(n)
        #phis = np.linspace(-0.9*np.pi,0.9*np.pi,n)
        phis = np.linspace(0,0.5*np.pi,n)
        for i,phi in enumerate(phis):
            res = sec.func(phi, v, omega)
            if not np.isnan(res):
                resid[i] = res
            else:
                resid[i] = 1e30
        i = np.argmin(abs(resid))
        return phis[i]

    def optimize_pitch(self):
        """
        Optimize rotor pitch for either maximum thrust (propeller) or maximum power (turbine)
        using a genetic evolution algorithm.

        This is intended as an example of how optimization can be done using the scipy.optimize
        package. The overall procedure can be readily modified to optimize for other parameters,
        e.g. a parametrized function for the pitch, or for a parameter sweep instead of 
        a single parameter set.

        return: Array of optimized pitches
        """

        def run_bemt(x):
            print('Current iteration:',x)
            for sec,pitch in zip(self.rotor.sections, x):
                sec.pitch = np.radians(pitch)

            T,Q,P,df = self.run()
            if self.mode == 'turbine':
                return -P
            else:
                return -T
 
        x = [sec.pitch for sec in self.rotor.sections]
        bounds = [(0,30)]*len(x)

        result = optimize.differential_evolution(run_bemt, bounds, tol=1e-1)

        return result 
    
    def OutputBlade(self,df,section_df,input_name,airfoil):       
        tsr=self.tsr
        twist=self.twist
        v_inf=self.v_inf
        v_rel=section_df[0]["v_rel"].tolist()
        radius=section_df[0]["radius"].tolist()
        chord=section_df[0]["chord"].tolist()
        dT=section_df[0]["dT"].tolist()
        dQ=section_df[0]["dQ"].tolist()
        
        file1 = open(airfoil,'r')
        lines = file1.readlines()
        airfoil_name=[None]*(len(lines))
        alpha=section_df[0]["alpha"].tolist()
        Cm=np.zeros(len(alpha))
        Pm=np.zeros(len(alpha))        

        Axial_force=np.zeros(len(dT)-1)
        Tangent_force=np.zeros(len(dQ)-1)
        Moment=np.zeros(len(Cm)-1)
        

        for i in range(len(lines)):
            a=lines[i].strip()
            airfoil_name[i]=a
            airfoil_name[i]="pybemt/airfoils/Cm/"+airfoil_name[i]+"_Cm.txt"
            file_alpha, file_Cm = np.loadtxt(airfoil_name[i], unpack=True)
            alpha[i]=degrees(alpha[i])
            for j in range (len(file_alpha)):
                temp=abs(file_alpha-alpha[i])
                pos=np.where(temp==np.amin(temp))
                index=int(pos[0])
                if file_alpha[index]>alpha[i]:
                    index = index-1
                x1=file_alpha[index]
                x2=file_alpha[index+1]
                y1=file_Cm[index]
                y2=file_Cm[index+1]
                x=alpha[i]
                y=y1+(x-x1)*(y2-y1)/(x2-x1)
            Cm[i]=y
            Pm[i]=0.5*1.225*(v_rel[i]**2)*(chord[i]**2)*Cm[i]
        
        for i in range(len(Axial_force)):
            dr = radius[i+1] - radius[i]
            ave_Pm = (Pm[i]+Pm[i+1])/2
            Axial_force[i]=dT[i+1]/3
            Tangent_force[i]=dQ[i+1]/(3*((radius[i]+radius[i+1])/2))
            Moment[i]=ave_Pm*dr


        blade_name=input_name
        f = open("output/"+blade_name+"-Load.txt", "w+")
        f.write(blade_name+'\n')
        f.write("TSR: "+str(tsr)+"\nPitch: "+str(twist)+"\nV_inf: "+str(v_inf)+"\n")
        f.write("Radius\tAxial Force\tTangential Force\tMoment\n")
        for i in range (len(Axial_force)):    
            f.write(str(round((radius[i]+radius[i+1])/2,3))+'\t'+str(round(Axial_force[i],3))+\
                '\t'+str(round(Tangent_force[i],3))+'\t'+str(round(Moment[i],3))+'\n')
        f.close()
    

      
