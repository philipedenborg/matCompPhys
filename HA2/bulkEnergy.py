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


a0 = 4.044  # lattice constant from EOS

al = bulk('Al','fcc',a0)

k = 16

# create calculator			
calc = GPAW(mode = PW(650),	# Basis with Cut-off			
		parallel={'band':1},	
		kpts = (k,k,k),		# Number of k-points
		xc = 'PBE',
		basis = 'dzp',		# xc functional
		txt = 'Al-bulkEnergy.txt',	# Name of output file
		eigensolver = 'rmm-diis')

al.set_calculator(calc) # set calculator

E_bulk = al.get_potential_energy()

if rank == 0:
    print E_bulk
