from ase.io import read
from ase import atoms
import numpy as np
import eam_calculator as EAM_calc
from scipy.optimize import least_squares
from ase.build import bulk
from ase.eos import EquationOfState
import matplotlib.pyplot as plt
plt.rc('text',usetex=True)
plt.rc('font', family='serif')
from ase.visualize import view

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

def check(params, E_0, a_0):
    
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

def RMSD(predictions,observed_data):
    Y_p = np.array(predictions)
    Y_o = np.array(observed_data)
    n = float(np.size(Y_o))
    MSE = 1/n*np.sum((Y_o-Y_p)**2)
    RMSD = np.sqrt(MSE)
    return RMSD

def check_errors(params,atoms,obs_forces):
    ons_forces = np.array(obs_forces)
    calc = EAM_calc.get_calc(params)
    atoms.set_calculator(calc)
    pred_forces = np.array(atoms.get_forces())
    rmsd_errors = [RMSD(pred_forces[:,i],obs_forces[:,i]) for i in xrange(0,3)]
    return rmsd_errors

def plot_comparison(params,atoms,obs_forces,err):
    plt.figure()
    plt.hold(True)
    x = np.linspace(-3,3)
    plt.plot(x,x,'black')
    calc = EAM_calc.get_calc(params)
    colours = ['g.','r.','b.']
    for j in [0,1,2]:
        F_obs = np.array(obs_forces[j])
        al = atoms[j]
        al.set_calculator(calc)
        pred_forces = np.array(al.get_forces())
        for i in [0,1,2]:
            plt.errorbar(pred_forces[:,i],F_obs[:,i],yerr=err[j][i],fmt=colours[j],ecolor='black',capthick=2)
    plt.xlabel('Predicted forces [Seriously do not know the unit]',fontsize=16)
    plt.ylabel('DFT data [UNIT UNKNOWN]',fontsize=16) 
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    plt.show()

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
weights = [1., 300., 200.] #Task 3

# Fit the parameters
if "fit" == "fitnot":
    fit = least_squares(fun = opt_func, x0 = params, args = (atoms, weights, forces, E_0, a_0)) 
    file = open("fit_params.txt", "w") 
    for i in fit.x: 
        file.write("%f\n" % i) 
        file.close()

params = [1464.484362, 3.007881, 8.221412, 1.19702] # fit if fitting is done #

# Check parameters by calculating E_co and a_0
if "check" == "checknot":
    a, E = check(params,E_0, a_0)
    print(a)
    print(E)

# Calculated RMSD for forces
if "rmsd" == "rmsdnot":
    f = open("RMSD.txt","w")
    for i in [0,1,2]:
        for x in check_errors(params,atoms[i],forces[i]):
            f.write("%.3f " % x)
            f.write("\n")

# Plot calculated forces against DFT result
if "plot" == "plot":        
    err = []
    f = open("RMSD.txt","r")
    for line in f:
        string = line.strip().split(" ")
        col = [float(s) for s in string]
        err.append(col)      
    plot_comparison(params,atoms,forces,err)

