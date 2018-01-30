#!/usr/bin/env python
from ase import Atoms
from gpaw import GPAW
from ase.build import bulk

atoms = bulk('Al')

calc = GPAW(mode = PW(500),
		h = 0.18,
		kpts = (4,4,4),
		xc = 'PBE',
		txt = 'out.txt',
		eigensolver = 'rmm-diis')  # create calculator

atoms.set_calculator(calc)
atom.get.potential_energy()

