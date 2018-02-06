from ase.build import fcc111, fcc100, add_adsorbate
from ase.visualize import view
from ase import Atoms
from gpaw import GPAW, PW

co = Atoms('CO',[(0.,0.,0.), (0, 0, 1.16)])
N = 1
k = 16

sites = ["ontop","bridge","hollow"]

f = open('site_energy.txt', 'w')

for x in sites:
    slab = fcc100('Al', size=(1,1,N))
    add_adsorbate(slab, co,1.7,x)
    slab.center(vacuum=10.0, axis=2)
    calc = GPAW(mode = PW(650),	# Basis with Cut-off			
		parallel={'band':1},	
		kpts = (k,k,1),		# Number of k-points
		xc = 'PBE',
		basis = 'dzp',		# xc functional
		txt = 'slab-%d.txt' %N,	# Name of output file
		eigensolver = 'rmm-diis')
    slab.set_calculator(calc)
    f.write(slab.get_potential_energy()+"\n")

f.close()

view(slab)
