#!/usr/bin/env python
from ase import Atoms
from gpaw import GPAW, PW
from ase.build import fcc111
from ase.build import fcc100
from ase.parallel import rank
import numpy as np

a0 = 4.04
slab1 = fcc111('Al',size = (20,20,3), a = a0)

for k in range(4,30,2):
	# create calculator			
	calc = GPAW(mode = PW(500),	# Basis with Cut-off			
			parallel={'band':1},	
			kpts = (k,k,k),		# Number of k-points
			xc = 'PBE',
			basis = 'dzp',		# xc functional
			txt = 'Al-slab1%02d.txt' % k,	# Name of output file
			eigensolver = 'rmm-diis')
	slab1.set_calculator(calc)
	E = slab1.get_potential_energy()	
	if rank == 0:	
		print k,E	

slab2 = fcc100('Al',size = (20,20,3), a = a0)

for k in range(4,30,2):
	# create calculator			
	calc = GPAW(mode = PW(500),	# Basis with Cut-off			
			parallel={'band':1},	
			kpts = (k,k,k),		# Number of k-points
			xc = 'PBE',
			basis = 'dzp',		# xc functional
			txt = 'Al-slab2%02d.txt' % k,	# Name of output file
			eigensolver = 'rmm-diis')
	slab2.set_calculator(calc)
	E = slab2.get_potential_energy()	
	if rank == 0:	
		print k,E	


