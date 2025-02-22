#!/usr/bin/env python
from ase import Atoms
from gpaw import GPAW, PW
from ase.build import bulk
from ase.parallel import rank
import numpy as np

a0 = 4.05
atoms = bulk('Al','fcc',a=a0)

for k in range(4,30,2):
	# create calculator			
	calc = GPAW(mode = PW(500),	# Basis with Cut-off			
			parallel={'band':1},	
			kpts = (k,k,k),		# Number of k-points
			xc = 'PBE',
			basis = 'dzp',		# xc functional
			txt = 'Al-%02d.txt' % k,	# Name of output file
			eigensolver = 'rmm-diis')
	atoms.set_calculator(calc)
	E = atoms.get_potential_energy()	
	if rank == 0:	
		print k,E	



