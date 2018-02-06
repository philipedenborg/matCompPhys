from ase.build import fcc111, add_adsorbate
from ase.visualize import view
slab = fcc111('Al', size=(1,1,1))
add_adsorbate(slab, 'CO')
slab.center(vacuum=10.0, axis=2)

view(slab)
