from ase.io import read
from ase import atoms
import numpy as np
import eam_calculator as EAM_calc
from scipy.optimize import least_squares
from ase.build import bulk
from ase.eos import EquationOfState

def opt_func(params, atoms, weights, forces, e_ref, a_ref):

    #Set calculator using given program and parameters
    calc = EAM_calc.get_calc(params)

    #Configurational testing for given data
    F_f = []
    for atom,force in zip(atoms, forces):
        atom.set_calculator(calc)
        F_f = np.array(weights[1]*(atom.get_forces() - force) / (force))
    
    #Aluminiun bulk initialization using for fitting
    Al = bulk('Al', 'fcc', a = a_ref)
    Al = set_calculator(calc)
    e = []
    v = []
    

    


atoms_90 = read(filename="snapshots_with_forces_xyz/res_POSCAR_0.9.xyz")
atoms_100 = read(filename="snapshots_with_forces_xyz/res_POSCAR_1.0.xyz")
atoms_110 = read(filename="snapshots_with_forces_xyz/res_POSCAR_1.1.xyz")

#Lists of forces and atoms
atoms = [atoms_90, atoms_100, atoms_110]
forces = [atoms_90.get_forces(), atoms_100.get_forces(), atoms_110.get_forces()]

#Initial parameters for the fit
params = [1000, 3, 5, 1] #eV, AA-1, AA, AA-1
a_0 = 4.032
E_0 = -3.36
w = [1, 1, 1] #Task 3

