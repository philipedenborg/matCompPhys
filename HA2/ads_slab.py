from gpaw import GPAW, PW
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111,fcc100, add_adsorbate
from ase.parallel import rank

h = 1.85
d = 1.16
k = 16
N = 1

surface = 100
site = "ontop"
pos = [(0., 0., 0.), (0., 0., d)] # Use to set orientation

if surface == 100:
	slab = fcc111('Al', size=(1, 1, N), vacuum=10.0)

calc1 = GPAW(mode = PW(650),	# Basis with Cut-off			
			parallel={'band':1},	
			kpts = (k,k,1),		# Number of k-points
			xc = 'PBE',
			basis = 'dzp',		# xc functional
			txt = 'ads-%d-%s-%d.txt' %(surface,site,N),	# Name of output file
			eigensolver = 'rmm-diis')

slab.set_calculator(calc1)
e_slab = slab.get_potential_energy()

calc2 = GPAW(mode='lcao',        # use the LCAO basis mode
            basis='dzp',        # use double zeta polarized basis set
            h=0.18,             # grid spacing
            xc='PBE',           # XC-functional
            #kpts=(16,16,16),    # k-point grid
            txt='out.txt')      # name of GPAW output text file

molecule = Atoms('CO', positions = pos )
a = 8
molecule.set_cell((a,a,a))
molecule.center()

molecule.set_calculator(calc2)
e_CO = molecule.get_potential_energy()


add_adsorbate(slab, molecule, h, site)
constraint = FixAtoms(mask=[a.symbol == 'Al' for a in slab])
slab.set_constraint(constraint)
dyn = QuasiNewton(slab, trajectory='AlCO.traj')
dyn.run(fmax=0.05)

if rank==0:
    print('Adsorption energy:', e_slab + e_CO - slab.get_potential_energy())
