# -*- coding: utf-8 -*-

#  __authors__ = Elena Khramenkova
#  __institution__ = TU Delft
#
import numpy as np
import math
np.set_printoptions(suppress=True)
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import seaborn as sns
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
from math import exp, log, factorial
import os

# The energies of the most stable systems (recombination energies)
ref = {
"CH4"     :     -24.06553,
"R"       :      8.3144598, #J/K*mol
"P_0"     :      0.986923169314,
 "P_1"    :      101325, #Paskal (1atm = 101325 Pa) kg*m^(-1)*s^(-2)
"k_b"     :      1.38064852E-23, # kg*m^2/(s^2*K)
"H2/MOR" :     -1733.253350871367047,
"H2O"    :     -17.220233011652674,
"O2"    :       -31.928779841938734,
"CuO"    :  np.array([-256.5142,-256.6045,-256.6052,-256.6060,-256.6064])/4,
"AlOH"    :     [-1750.8553204661, 1, 1, 1], # AlOH confined in the side posket of mordenite from Harshini's paper
"Cu3O3"   :     [-1924.612981, 3, 0, 3, 0],  # the lowest model from the JPC letters CumAlnOpHq
"Cu2AlO4H"   :    [-1895.341404, 2, 1, 4, 1], # Cu2AlO4H from 
}

 #### FORMATION OF Cu2AlO4H in the presence of CuO and EFAl
 #### 2x*CuO + (x-1)H2/MOR + AlxOyHz + (5x + Z -2y + 2)/4 *O2 -> x*Cu2AlO4H + (x - 2 + z)/2 *H2O

class AITD:     
    #Let's assume that Sconf and n depend on the concentration of CuO and EFAL
    ### n (number of defects) ~ 0.5*c(Al) - c(Cu) and S_conf ~ c 0.5*c(Al) - c(Cu)
    ##### BAS (mikromol/g) 1130 of H-MOR-B)
    ##### LAS_8MR_bottom (mikromol/g) 255 of H-MOR-B)
    ##### LAS_total 490 (miromol/g)
    cl_conc_Al = 260e-6 #mikromol/g
    cl_conc_Cu = np.arange(1e-6, 500e-6, 10e-6) #concentration in mol/gramm
    # n is the number of defects
    # n = v*Na # where v is the value form the range [0; 180]mikromol/gramm

    # n = [i*ref["Na"] for i in cl_conc] #the number of defects depending on the concentration of cluster

    # N = mol*Na = m/M * Na = 1gramm/M(Al2Si46O96) * Na
    N = 1*ref["Na"]/2881.84 #dimensionless

    # range of water chemical potentials
    mu_H2O = np.arange(-2.0, -0.1, 0.038)
    
    # range of oxygen chemical potentials
    mu_O2 = np.arange(-2.0, -0.1, 0.038)
    # upload the tabulated data for H2O
    list_H2O = np.loadtxt("list_H2O.txt")
    H_H2O_0K, H_H2O_700, S_H2O_700 = list_H2O[0,2], list_H2O[8,2], list_H2O[8,1]
    
    # upload the tabulated data for O2
    list_O2 = np.loadtxt("list_O2.txt")
    H_O2_0K, H_O2_700, S_O2_700 = list_O2[0,2], list_O2[8,2], list_O2[8,1]
    
    def __init__(self,  reference, cluster):
        self.reference = reference
        self.cluster = cluster
    
    def calculate_E(self): # CumAlnOpHq
        #The equation for EFAl AlOH/zeol
        # (2*p - q + 2 - 3*n - 2*m) / 4 * O2 + (q + n - 2) / 2 * H2O + n * AlOH/zeol + m * CuO + (1 - n) H2/zeol ->
        # -> CumAlnOpHq/zeol
        # convert a.u. to eV
        DE = 27.211 * ( self.cluster[0] - ( ( ( 2*self.cluster[3] - self.cluster[4] + 2 - 3*self.cluster[2] - 2*self.cluster[1] ) / 4 ) * self.reference['O2'] ) \
            - ( ( ( self.cluster[4] + self.cluster[2] - 2 ) / 2 ) * self.reference['H2O'] ) - self.cluster[2] * self.reference['AlOH'][0] \
                - ( 1 - self.cluster[2] ) * self.reference['H2/MOR'] - self.cluster[1] * self.reference['CuO'][4] )

        return DE
              
    def calculate_3d_G(self, DE):
        # X = oxygen chemical potential change
        # Y = concentrational change of Cu content
        X, Y = np.meshgrid(AITD.mu_O2, AITD.cl_conc_Cu)

        ### GIBBS FREE ENERGY CHANGE: delta_G = delta_H * n - n * S_vib - T * delta_S_conf can be approximated ->
        # -> DG = ( n * DH - T * DS ) / N
        # 700 K and 1 atm for X definition (water chemical potential)
        mu_H2O_700K = ( ( ( AITD.H_H2O_700 - AITD.H_H2O_0K ) * 1000 ) - ( 700 * AITD.S_H2O_700 ) + \
                       self.reference["R"] * 700 * math.log( 1 / self.reference["P_0"] ) ) / 96485 # eV of water chemical potentail
            
        # Add the change in the oxygen chemical potentials and constant 
        DH =  DE - ( ( 2*self.cluster[3] - self.cluster[4] + 2 - 3*self.cluster[2] - 2*self.cluster[1] ) / 2 ) * AITD.mu_O2 - \
                     ( ( self.cluster[4] + self.cluster[2] - 2 ) / 2 ) * mu_H2O_700K

        # Number of defects depends on the concetration of Ga and Ca.      
        # If there is Al then do the following: Al_concentration - Cu _concentration
        if self.cluster[2] > 0:
            n = ( ( AITD.cl_conc_Al - Y ) * self.reference["Na"] )
            # configurational entropy
            # Sconf = kb * ln(N + n)/(N!n!)
            # S_conf = [ref["k_b"]*N*(math.log(1+(l/N)) + (l/N) * math.log(1+N/l)) for l in n] # eV/K #DFT and Thermodynamics
            DS = self.reference["k_b"] * AITD.N * np.log( ( AITD.N + self.reference["Na"] * ( AITD.cl_conc_Al - Y ) ) / AITD.N ) + \
            ( ( self.reference["Na"] * ( AITD.cl_conc_Al - Y ) ) / AITD.N ) * \
            np.log( ( AITD.N + self.reference["Na"] * ( AITD.cl_conc_Al - Y ) )/( self.reference["Na"] * ( AITD.cl_conc_Al - Y ) ) ) # eV/K #DFT and Thermodynamics
        else:
        # If there is no Ca then nly consider the concentration of Ga
            n  = ( Y * self.reference["Na"] )
            DS = self.reference["k_b"] * AITD.N * np.log( ( AITD.N + self.reference["Na"] * Y ) / AITD.N ) + \
                ( ( self.reference["Na"] * Y ) / AITD.N ) * np.log( ( AITD.N + self.reference["Na"] * Y ) / ( self.reference["Na"] * Y)) # eV/K #DFT and Thermodynamics
        DG = (  n * DH - 700 * DS ) / AITD.N

        return DG
        
    def print_3d_diagram(self, list_DG, list_names, list_colors):
        # Create the figure
        fig = plt.figure()
        # Add the axes
        ax = fig.gca(projection='3d')
        X, Y = np.meshgrid(AITD.mu_O2, AITD.cl_conc_Cu)
        # Add the axes
        # Plotting the planes using 2d mesh and 
        # Iterating trough the list of structures to keep them together in one figure
        for DG, name, color in zip(list_DG, list_names, list_colors):
            surf = ax.plot_surface(X, Y, DG, color=color, alpha = 1.0, rstride=20, cstride=20, label = name)
            surf._facecolors2d = surf._facecolor3d
            surf._edgecolors2d = surf._edgecolor3d
            ax.set_xlabel(r'water chem. pot., eV')        
            ax.set_ylabel(r'c(Cu), mol')
            ax.set_zlabel('$ΔG, eV$')
            ax.set_title(r'$Phase$ $diagram$')
            plt.rcParams['svg.fonttype'] = 'none'
            lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0, prop={'size': 10})
            ax.view_init(0, 0)
        # Save the graph
        plt.savefig("Phase_diagram_Cu2AlO4H_vs_Cu3O3.svg", dpi=300)
        plt.show()

    def calculate_2d_G(self, DE):
        # 700 K and 1 atm for water chemical potential
        mu_H2O = ( ( ( AITD.H_H2O_700 - AITD.H_H2O_0K ) * 1000 ) - ( 700 * AITD.S_H2O_700 ) + \
                       self.reference["R"] * 700 *\
                math.log( 1 / self.reference["P_0"] ) ) / 96485 # eV of water chemical potentail
        # 700 K and 1 atm for water chemical potential
        mu_O2 = ( ( ( AITD.H_O2_700 - AITD.H_O2_0K ) * 1000 ) - ( 700 * AITD.S_O2_700 ) + \
                       self.reference["R"] * 700 *\
                math.log( 1 / self.reference["P_0"] ) ) / 96485 # eV of oxygen chemical potentail
        
        ### GIBBS FREE ENERGY CHANGE: delta_G = delta_H * n - n * S_vib - T * delta_S_conf -> 
        # DG = ( n * DH - T * DS ) / N
        
        # add the change in the water and oxygen chemical potentials
        DH =  DE - ( ( ( 2*self.cluster[3] - self.cluster[4] + 2 - 3*self.cluster[2] - 2*self.cluster[1] ) / 2 ) * mu_O2 - \
                     ( ( self.cluster[4] + self.cluster[2] - 2 ) / 2 ) * mu_H2O ) 

        # Number of defects depends on the concetration of Ga and Ca.      
        # If there is Ca then do the following: Ga_concentration - Ca_concentration
        Y = AITD.cl_conc_Cu
        if self.cluster[2] > 0:
            n = ( ( AITD.cl_conc_Al - Y ) * self.reference["Na"] )
            #configurational entropy
            #Sconf = kb * ln(N + n)/(N!n!)
            #S_conf = [ref["k_b"]*N*(math.log(1+(l/N)) + (l/N) * math.log(1+N/l)) for l in n] # eV/K #DFT and Thermodynamics
            DS = self.reference["k_b"] * AITD.N * np.log( ( AITD.N + self.reference["Na"] * ( AITD.cl_conc_Al - Y ) ) / AITD.N ) + \
            ( ( self.reference["Na"] * ( AITD.cl_conc_Al - Y) ) / AITD.N ) * np.log( ( AITD.N + self.reference["Na"] * \
             ( AITD.cl_conc_Al - Y ) ) / ( self.reference["Na"] * ( AITD.cl_conc_Al - Y ) ) ) # eV/K #DFT and Thermodynamics
        else:
        # If there is no Al then only consider the concentration of Ga
            n  = ( Y * self.reference["Na"] )
            DS = self.reference["k_b"] * AITD.N * np.log( ( AITD.N + self.reference["Na"] * Y ) / AITD.N ) + \
            ( ( self.reference["Na"] * Y ) / AITD.N ) * np.log( ( AITD.N + self.reference["Na"] * Y ) / ( self.reference["Na"] * Y) ) # eV/K #DFT and Thermodynamics
        DG = (  n * DH - 700 * DS ) / AITD.N

        return DG
        
    def print_2d_diagram(self, list_DG, list_names, list_colors):
        os.getcwd()
        sns.set_style('whitegrid')
        for DG, name, color in zip(list_DG, list_names, list_colors):
            plt.plot(AITD.cl_conc_Cu, DG, color=color, label=name)
            plt.xlabel(r'c(Cu), mol')        
            plt.ylabel(r'Gibbs free energy change, eV')
            plt.title(r'Phase diagram CuAl vs Cu')
            plt.xlim(0, 0.0005)
            plt.rcParams['svg.fonttype'] = 'none'
            plt.legend()

        #plt.savefig("\concentration_dependencies_CuAl_vs_Cu_2d.svg", dpi=300)
        plt.savefig("concentr_dep_CuAl_vs_Cu_2d.svg", dpi=300)
        plt.show()
        

if __name__ == '__main__':
    ### create the objects
    aitd_Al = AITD(ref, ref['Cu2AlO4H'])
    aitd = AITD(ref, ref['Cu3O3'])
    ### Cu/Al-containing plane (3D) ###
    ##calculate the Gibbs free enegry change (DG) in eV
    free_energies_Al = aitd_Al.calculate_3d_G(aitd_Al.calculate_E())
    
    ### only Cu-containing plane ###
    free_energies = aitd.calculate_3d_G(aitd.calculate_E())

    ### Plot all the planes together ###
    aitd_Al.print_3d_diagram([free_energies_Al, free_energies],\
                             ['CuAlO4H', 'Cu3O3'], ['lawngreen', 'red'])
        
    ### 2D representation
    ##calculate the Gibbs free enegry change (DG) in eV
    free_energies_Al_2d = aitd_Al.calculate_2d_G(aitd_Al.calculate_E())
    
    ### only Ga-containing plane ###
    free_energies_2d = aitd.calculate_2d_G(aitd.calculate_E())
    
    ### Plot all the planes together ###
    aitd.print_2d_diagram([free_energies_Al_2d, free_energies_2d],\
                             ['CuAlO4H', 'Cu3O3'], ['lawngreen', 'red'])
