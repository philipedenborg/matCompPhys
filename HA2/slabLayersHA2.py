#!/usr/bin/env python
from ase import Atoms
from gpaw import GPAW, PW
from ase.build import fcc111
from ase.build import fcc100
from ase.parallel import rank
import numpy as np

a0 = 4.04
k = 16 

for N in range(1,23,2):
        fcc = fcc100('Al', (1, 1, N), a=a0, vacuum=7.5)
        fcc.center(axis=2)
	# create calculator			
	calc = GPAW(mode = PW(650),	# Basis with Cut-off			
			parallel={'band':1},	
			kpts = (k,k,1),		# Number of k-points
			xc = 'PBE',
			basis = 'dzp',		# xc functional
			txt = 'slab-%d.txt' %N,	# Name of output file
			eigensolver = 'rmm-diis')
	fcc.set_calculator(calc)
	E = fcc.get_potential_energy()	
	


	


