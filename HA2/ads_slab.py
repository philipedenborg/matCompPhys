from gpaw import GPAW, PW
from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111,fcc100, add_adsorbate
from ase.parallel import rank
from ase.visualize import view

h = 1.85	# Off-set of adsorbate
d = 1.16	# Bond length of adsorbate
k = 5		# Number of k-points
N = 1		# Number of cell layers in slab
vac = 10.0	# Vacuum 

surface = 100	# 100 or 111
site = "ontop"  # ontop, bridge or hollow for 100 and ontop or bridge for 111 (fcc or hcp)
pos = [(0., 0., 0.), (0., 0., d)] # Use to set orientation (000, 000d) = carbon facing slab
orientation = "C"  # Atom facing the slab, C or O.

if surface == 100:
	slab = fcc100('Al', size=(1, 1, N), vacuum=vac)
if surface == 111:
	slab = fcc111('Al', size=(1, 1, N), vacuum=vac)

calc1 = GPAW(mode = PW(650),	# Basis with Cut-off			
			parallel={'band':1},	
			kpts = (k,k,1),		# Number of k-points
			xc = 'PBE',
			basis = 'dzp',		# xc functional
			txt = 'ads-%d-%s-%s-%d.txt' %(surface,site,orientation,N),	# Name of output file
			eigensolver = 'rmm-diis')

slab.set_calculator(calc1)
e_slab = slab.get_potential_energy()

calc2 = GPAW(mode='lcao',        # use the LCAO basis mode
            basis='dzp',        # use double zeta polarized basis set
            h=0.18,             # grid spacing
            xc='PBE',           # XC-functional
            #kpts=(16,16,16),    # k-point grid
            txt='CO_gpaw.txt')      # name of GPAW output text file

molecule = Atoms('CO', positions = pos )
a = 8				# Size of molecule cell
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
    print(surface,site,N,'Adsorption energy:', slab.get_potential_energy() - e_slab - e_CO )
