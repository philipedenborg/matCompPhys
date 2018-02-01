#!/usr/bin/env python
from ase import Atoms
from gpaw import GPAW, PW
from ase.build import bulk
from ase.parallel import rank
import numpy as np

a0 = 4.05
atoms = bulk('Al','fcc',a=a0)
p = 0.10

for k in range(4,30,2):
	# create calculator			
	calc = GPAW(mode = PW(500),	# Basis with Cut-off			
			h = 0.18,	
			kpts = (k,k,k),		# Number of k-points
			xc = 'PBE',		# xc functional
			txt = 'Al-%02d.txt' % k,	# Name of output file
			eigensolver = 'rmm-diis')
	atoms.set_calculator(calc)
	if rank == 0:
		print k,atom.get.potential_energy()			



