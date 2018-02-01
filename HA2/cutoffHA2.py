#!/usr/bin/env python
from ase import Atoms
from gpaw import GPAW, PW
from ase.build import bulk
from ase.parallel import rank
import numpy as np
from ase.units import kJ


a0 = 4.05  # approximate lattice constant

al = bulk('Al','fcc',a0)

k = 16

# create calculator
for cf in np.linspace(400,1000,6):
    calc = GPAW(mode = PW(cf),	# Basis with Cut-off			
		parallel={'band':1},	
		kpts = (k,k,k),		# Number of k-points
		xc = 'PBE',
		basis = 'dzp',		# xc functional
		txt = 'Al-cf%02d.txt' % cf,	# Name of output file
		eigensolver = 'rmm-diis')
    al.set_calculator(calc) # set calculator
    E =  al.get_potential_energy()
    if rank == 0:
        print cf,E
