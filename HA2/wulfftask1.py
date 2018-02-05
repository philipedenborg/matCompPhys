from ase.cluster.wulff import wulff_construction
from ase.io import write
from ase.visualize import view

atoms = wulff_construction('Al',
                           surfaces=[(1,0,0),
                                     (1,1,1),
                                     (1,1,0)],
                           energies=[0.8128, 0.8128, 0.8128],
                           size = 10000,
                           structure = 'fcc',
                           roundering='below')
atoms.center(vacuum=10)
view(atoms)
