#!/usr/bin/env python
from ase import Atoms
from gpaw import GPAW, PW
from ase.build import bulk
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
			txt = 'k_data.txt',	# Name of output file
			eigensolver = 'rmm-diis')
	atoms.set_calculator(calc)		
	print k,atom.get.potential_energy()			
	


