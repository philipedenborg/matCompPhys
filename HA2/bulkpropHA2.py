#!/usr/bin/env python
from ase import Atoms
from gpaw import GPAW, PW
from ase.build import bulk
from ase.parallel import rank
import numpy as np
from ase.eos import EquationOfState
from ase.io.trajectory import Trajectory
from ase.io import read
from ase.units import kJ


a0 = 4.0  # approximate lattice constant
al = bulk('Al','fcc',a0)

k = 10

# create calculator			
calc = GPAW(mode = PW(500),	# Basis with Cut-off			
		parallel={'band':1},	
		kpts = (k,k,k),		# Number of k-points
		xc = 'PBE',
		basis = 'dzp',		# xc functional
		txt = 'Al-%02d.txt' % k,	# Name of output file
		eigensolver = 'rmm-diis')

al.set_calculator(calc) # set calculator

cell = al.get_cell()
traj = Trajectory('Al.traj', 'w')
for x in np.linspace(0.95, 1.05, 5):
    al.set_cell(cell * x, scale_atoms=True)
    al.get_potential_energy()
    traj.write(al)

configs = read('Ag.traj@0:5')  # read 5 configurations
# Extract volumes and energies:
volumes = [al.get_volume() for al in configs]
energies = [al.get_potential_energy() for al in configs]
eos = EquationOfState(volumes, energies)
v0, e0, B = eos.fit()
print(B / kJ * 1.0e24, 'GPa')
eos.plot('Ag-eos.png')



