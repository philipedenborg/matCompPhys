# This code is done by Philip Edenborg and Magnus Fant
# Import libs
import numpy as np
from ase.calculators.lj import LennardJones
from ase import Atoms
from ase.visualize import view
import matplotlib.pyplot as plt
from scipy import optimize

# Lennard Jones potential
def lennard_jones(epsilon, sigma, r1, r2):
    R = 0
    for i,j in zip(r1,r2): # Calculate distance between atoms
        R += (i-j)**2
    R = np.sqrt(R)
    V = 4*epsilon*((sigma/R)**12 - (sigma/R)**6)
    return V

# Lennard Jones potential as function of distance
def lennard_jones2(R):
    epsilon = 0.0104
    sigma = 3.40
    V = 4*epsilon*((sigma/R)**12 - (sigma/R)**6)
    return V

# Calc potential energy of system of atoms using the Lennard Jones potential
def calc_E(positions, epsilon, sigma):
    E = 0
    N = np.size(positions,0) # Number of atoms in system
    for i in range(0,N-1):
        for j in range(i+1,N): # Count interaction between atoms only once
            E += lennard_jones(epsilon, sigma, positions[i], positions[j])
    return E


# Read file with positions to numpy array
pos =[]
with open("positions.txt") as f:
    for line in f:
        string = line.strip().split(',')
        col = [float(s) for s in string]
        pos.append(col)
pos = np.array(pos)

# Physical parameters
epsilon = 0.0104
sigma = 3.40

# Initiate atoms and calculator using ASE
atoms = Atoms('ArArArAr', positions = pos) # Four Argon atoms
calc = LennardJones(epsilon = 0.0104, sigma = 3.4, rc=100)
atoms.set_calculator(calc)  # attach calc to Atoms object

# Calculate energy using calc_E and ASE
E1 = calc_E(pos,epsilon,sigma)
E2 = atoms.get_potential_energy()  # calculate energy

# Code for finding equilibrium distance between two Argon atoms in Lennard Jones potential
N = 2000    # Number of data points
V_LJ = np.array([])     # Value for potential for different distances
R = np.linspace(3.5,4,N)    # Distance between atoms, set small range for higher accuracy in minimum

for x in range(0,N):
    V_LJ = np.append(V_LJ,lennard_jones(epsilon, sigma,[0,0,0], [R[x],0,0]))

# Find minimum in energy using two different methods
ind = np.argmin(V_LJ)   # Find index corresponding to lowest energy
minimum = optimize.fmin(lennard_jones2,1)    # Find minimum using numpy optimization method

# Print results
print("The equilibrium distance (when using argmin) is: R_0 = {0:.5f} Å".format(R[ind]))
print("The equilibrium distance (when using fmin) is: R_0 = {0:.5f} Å".format(minimum[0]))
print("The energy (as calculated by you) for the four Argon atoms is E = {0:.9f} eV".format(E1))
print("The energy (as calculated by ASE) for the four Argon atoms is E = {0:.9f} eV".format(E2))

# Plot Lennard Jones potential and visualize atomic system
R = np.linspace(3.1,7,100)
plt.figure(figsize=(12,8))
plt.plot(R,lennard_jones2(R),'g-')
plt.xlabel('x-axis',fontsize=20)
plt.ylabel('y-axis',fontsize=20)
plt.show()
#view(atoms)
