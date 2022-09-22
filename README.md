# CONCENTRATION DEPENDENCIES

This script can be executed if the respective enregies of the metal-modified zeolites, respective metal bulk and other components of the reaction.
This approach is based on the ab initio thermodynamics analysis (aiTA) and defects formation in zeolites.
Hence, the Cu-oxo and CuAl-oxo clusters in mordenite are treated as defects.

## Installation

```
https://github.com/ekhramen/concentration_dependencies.git

```

* [example]: This folder contains working script with inserted value corresponding to the formation Cu2AlO4H and Cu3O3 complexes confined in mordenite.
Formation of the CuAl-oxo complexes in mordenite in the presence of CuO and extraframework aluminium species in mordenite.
Additionally, water and oxygen are present. The corresponding equilibrium is following:
```
(2p - q + 2 - 3n - 2m)/4 * O2 + (q + n - 2)/2 * H2O + n * AlOH/zeol + m * CuO + (1 - n) H2/zeol -> CumAlnOpHq/zeol
```
The GIBBS FREE ENERGY CHANGE:
```
delta_G = delta_G_1def * n - T * delta_S_conf
```
where n is the number of defects per 1 gramm of catalyst, delta_G_1def is the Gibbs free energy formation of a single defect.
S_vib is the vibrational entropy, T is the temperature, delta_S_conf is configurational entropy.
```
Sconf = kb * ln(N + n)/(N!n!)
```
The Gibbs free energy of asingle defect formation can be defined as following:
```
delta_G_1def =  CumAlnOpHq/zeol - (2p - q + 2 - 3n - 2m)/4*O2 - (q + n - 2)/2*H2O - n*AlOH/zeol - m*CuO - (1 - n) H2/zeol - (2p - q + 2 - 3n - 2m)/2*delta_mu(O) - (q + n - 2)/2*delta_mu(H2O)
```

* [scripts]: This folder contains the scripts used to run the aiTA.

## Instructions

This script contains the following steps with regard to the present example:
1. Insert the energies of the all components of the aiTA equations in the dictionary.
2. In class AITD, define N (the number of the framework units), concetration of the Cu,
n (the number of defects as a function of the cluster), concetration of Al.
3. For the gaseous species (molecular oxygen, water), the chemical potentials have to be defined. 
The thermochemical data on other gaseous species could be found here:
        ```
        https://webbook.nist.gov/chemistry/fluid/
        ```
4. Using the equation where AlOH/zeol is used as a EFAl, the reaction energy of the CuAl-oxo and Cu-oxo clusters formation.
Function calcualte_E gives the reaction energy.
5. Using the tabulated data for the gaseous species (molecular oxygen and water) and temperature of 700 K, the Gibbs free enregy of a single defect formation is calculated.
Then the number of the defects is itroduced through the concentration dependecies and configurationa; entropy.
Function calculate_3d_G calculates the Gibbs free energy as a function of Cu concentration.
6. Function print_3d_diagram prints a 3D diagram, where Gibbs free energy is a funcion of Cu concentration and oxygen chemical potantial.
7. Function calculate_2d_G the Gibbs free energy as a function of Cu concentration, while the chemical potentials are defined at 700 K.
8. Function print_2d_diagram prints a 2D diagram, where Gibbs free energy is a function of Cu concentration.
