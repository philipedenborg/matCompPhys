from ase import Atoms
from gpaw import GPAW, PW
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111,fcc100, add_adsorbate
from ase.parallel import rank

f100=open("slab-19-100.txt")
f111=open("slab-19-111.txt")

f=open("E_ads","w")

E_slab100 = f100.get_potential_energy()
E_slab111 = f111.get_potential_energy()

N = 19
slab111 = fcc111('Al', size=(1, 1, N), vacuum=10.0)
slab100 = fcc100('Al', size=(1, 1, N), vacuum=10.0)

co = Atoms('CO')

k = 16
calc = GPAW(mode = PW(650),	# Basis with Cut-off			
		parallel={'band':1},	
		kpts = (k,k,1),		# Number of k-points
		xc = 'PBE',
		basis = 'dzp',		# xc functional
		txt = 'slab-%d.txt' %N,	# Name of output file
		eigensolver = 'rmm-diis')
co.set_calculator(calc)

E_co = molecule.get_potential_energy()
sites = ["ontop","hollow"]
surf = ["100","111"]
slabs = ["100":slab100,"111":slab111]
energies = ["100":E_slab100,"111":E_slab111]
offset = 2
for s1,s2 in zip(sites,surf):
    add_adsorbate(slabs[s2], molecule, offset, s1)
    slab.set_calculator(calc)
    constraint = FixAtoms(mask=[a.symbol == 'Al' for a in slabs[s2]])
    slabs[s2].set_constraint(constraint)
    dyn = QuasiNewton(slabs[s2], trajectory="COAl%s%s.traj" % s1,s2)
    dyn.run(fmax=0.05)
    E_slab_ads = slabs[s2].get_potential_energy()
    E = energies(s2) + E_co - E_slab_ads
    
    if rank == 0:
        f.write("%s, %s %f" % s1, s2, E)




