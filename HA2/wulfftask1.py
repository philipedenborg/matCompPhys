from ase.cluster.wulff import wulff_construction
from ase.io import write
from ase.visualize import view

case = 1	# 1 for clean surfaces, 2 for adsorbation 
N = 3000	# Number of atoms 
vac = 10.0	# Vacuum
if case == 1:
	#f100=open("surfaceenergy100.txt","r")
	#f111=open("surfaceenergy111.txt","r")
	E100 = 0.95	#float(f100.read())
	E111 = 0.82	#float(f111.read())
if case == 2:
	#f100=open("surfaceenergy100ads.txt","r")
	#f111=open("surfaceenergy111ads.txt","r")
	E100 = 1.665769	#float(f100.read())
	E111 = 1.708100 	#float(f111.read())


atoms = wulff_construction('Al',
                           surfaces=[(1,0,0),
                                     (1,1,1)],
                           energies=[E100, E111],
                           size = N,
                           structure = 'fcc',
                           rounding='below')
atoms.center(vacuum=vac)
view(atoms)
