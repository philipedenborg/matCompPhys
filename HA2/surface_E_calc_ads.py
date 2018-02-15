import numpy as np
from ase.io import read
from ase.units import J, m

case = 111

if case == 100:
	f = open('surfaceenergy100adsN5ot.txt', 'w')
	E_ads = 1.46
	fcc = read('slab-19-100.txt')
if case == 111:
	f = open('surfaceenergy111adsN5ot.txt', 'w')
	E_ads = 1.57
	fcc = read('slab-19-111.txt')

E_bulk = -3.737236
N = 19
cell = fcc.get_cell()
area = np.linalg.norm(np.cross(cell[0],cell[1]))
print area
sigma = 1/(2*area)*(fcc.get_potential_energy()-N*E_bulk)
theta = 1.0/4.0
sigma += theta*E_ads/area
f.write(str(sigma/(J/m**2)) + '\n')
f.close()
