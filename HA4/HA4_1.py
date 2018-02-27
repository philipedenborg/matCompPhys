from ase.io import read
from ase import atoms
import numpy as np
import eam_calculator as EAM_calc
from scipy.optimize import least_squares
from ase.build import bulk
from ase.eos import EquationOfState

def opt_func(params, atoms, weights, forces, E_0, a_0):

    #Set calculator using given program and parameters
    calc = EAM_calc.get_calc(params)
    
    #Calculate force fit for optimization function
    F_f = []
    for atom,force in zip(atoms,forces):
        atom.set_calculator(calc)
        #F_i = np.array(np.sqrt(weights[0])*(atom.get_forces()-force)).flatten()
        F_i = np.sqrt(weights[0])*(atom.get_forces()-force)
        F_f.extend(F_i)
    
    
    #Initialize Al for energy and lattice fit    
    al = bulk('Al', 'fcc', a = a_0)
    al.set_calculator(calc)
    e = []
    v = []
    
    cell = al.get_cell()
    for x in np.linspace(.5, 1.5, 10):
        al.set_cell(cell*x, scale_atoms=True)
        e.append(al.get_potential_energy())
        v.append(al.get_volume())

    eos = EquationOfState(v, e)
    V, E, B = eos.fit()

  
    #Create lattice parameter fit for optimization function
    a = np.cbrt(V*4.0)
    F_l = np.sqrt(weights[1])*(a - a_0)

    #Create energy fit for optimization function
    F_e  = np.sqrt(weights[2])*(E-E_0)

    #Full optimization function
    F = np.append(F_f,[F_l, F_e])

    return F

def check(params, forces, E_0, a_0):
    
    calc = EAM_calc.get_calc(params)

    al = bulk('Al', 'fcc', a = a_0)
    al.set_calculator(calc)
    e = []
    v = []
    
    cell = al.get_cell()
    for x in np.linspace(.5, 1.5, 10):
        al.set_cell(cell*x, scale_atoms=True)
        e.append(al.get_potential_energy())
        v.append(al.get_volume())

    eos = EquationOfState(v, e)
    V, E, B = eos.fit()

    a = np.cbrt(V*4.0)

    return a, E
    
atoms_90 = read(filename="snapshots_with_forces_xyz/res_POSCAR_0.9.xyz")
atoms_100 = read(filename="snapshots_with_forces_xyz/res_POSCAR_1.0.xyz")
atoms_110 = read(filename="snapshots_with_forces_xyz/res_POSCAR_1.1.xyz")

#Lists of atoms and forces
atoms = [atoms_90, atoms_100, atoms_110]
forces = [atoms_90.get_forces(), atoms_100.get_forces(), atoms_110.get_forces()]

#Initial parameters for the fit
params = [1000, 3, 5, 1] #eV, AA-1, AA, AA-1
a_0 = 4.032
E_0 = -3.36
weights = [1., 100., 200.] #Task 3
#'''
fit = least_squares(fun = opt_func, x0 = params, args = (atoms, weights, forces, E_0, a_0)) #Activated when fit is needed

file = open("fit_params.txt", "w") 
for i in fit.x: 
    file.write("%f\n" % i) 
file.close()
#print(fit)
#'''
params = fit.x #[3436.87, 3.490, 4.387265, 0.86178] #fit if fitting is done #

a, E = check(params, forces, E_0, a_0)

print(a)
print(E)

