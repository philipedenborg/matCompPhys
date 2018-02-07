import numpy as np
from ase.io import read
from ase.units import J, m

f = open('surfaceenergy111.txt', 'w')
E = -3.737236
for N in range(1,17,2):
    fcc = read('slab-'+str(N)+'.txt')
    cell = fcc.get_cell()
    area = np.linalg.norm(np.cross(cell[0],cell[1]))
    sigma = 1/(2*area)*(fcc.get_potential_energy()-N*E)
    f.write(str(N) + ' ' + str(sigma/(J/m**2)) + '\n')
f.close()
