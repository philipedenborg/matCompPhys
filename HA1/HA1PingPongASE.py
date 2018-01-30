import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.visualize import view
from ase.build import bulk
from ase.build import surface
from ase.io import read,trajectory


d = 1.1
co = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)])

print co.get_chemical_symbols()

#view(co)

mol = molecule('CH3CH2OH')

#view(mol)

struct = bulk('Fe','bcc',a=2.87,cubic=True)

#view(struct)

slab = surface('Au', (2, 1, 1), 9)
slab.center(vacuum=10, axis=2)

#view(slab)

from ase import Atoms
from ase.io import read,write

co = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.1)])
write('co.xyz',co)

write('traj.traj',co)
