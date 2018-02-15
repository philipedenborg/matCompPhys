import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import scipy.sparse as sparse
from HA3_2 import create_second_derivative_matrix_1D
def main():

    ''' Finite difference geometry (1D) '''
    r_max = 8 # Maximum radius of position grid in Hartree units
    r_min = 0 # Minimum radius of position grid in Hartree units
    n_r = 1000 # Number of elements in position grid
    r = np.linspace(r_min, r_max, n_r+1) # Position grid in Hartree units
    r = r[1:] # Remove singularity in r=0

    ''' Differentiation '''
    h = r[1]-r[0] # Step size
    A_dd = create_second_derivative_matrix_1D(r, h)

    ''' Physical constants '''
    Z_helium = 2 # Charge of helium nucleus in hartree units
    Z_hydrogen = 1 # Charge of hydrogen nucleus in hartree units

    ''' Boundary conditions '''
    U_0 = 0
    U_inf = 0

    ''' Analytical solutions '''
    # V_hartree = 1/r - (1 + 1/r)*np.exp(-2*r) # The hartree potential for hydrogen
    phi_s_H = (1/np.sqrt(np.pi))*np.exp(-r) # Wave function for hydrogen in hartree
    # f = np.sqrt(4*np.pi)*phi_s_H*r
    n_s_H = phi_s_H*phi_s_H # Electron density hydrogen ground state

    ''' Initial conditions for U '''
    # R_U = np.zeros(np.size(r))
    # R_U[0] = U_0
    # R_U[-1] = U_inf

    ''' Initial condition for f '''
    # R_f = np.zeros(np.size(r))
    # R_f[0] = 0
    # R_f[-1] = 0

    ''' Computation '''
    Z = Z_hydrogen
    I = np.identity(n_r)
    V_H = 0
    V_x = 0
    V_c = 0
    A = (-1/2.0)*A_dd - I*(Z/r) - V_H - V_x - V_c
    (energy, wave_functions) = sparse.linalg.eigs(A, which='SM') # SM = smallest
    # magnitude of the eigenvectors
    e_min_ind = np.argmin(energy) # index of lowest energy in the energy vector
    energy_min = energy[e_min_ind]
    phi_min = wave_functions[:, e_min_ind] # Find the wave function set corresponding
    # to the lowest energy
    phi_min = phi_min/np.sqrt(np.trapz(phi_min**2, r)) # normalization
    print(energy_min)

    ''' Plotting '''
    fig_1 = plt.figure()
    ax_potential = fig_1.add_subplot(111)
    label1 = 'Calculated radial probability distribution. \n Energy: ' \
           + str(np.round(np.real(energy_min),6)) + ' atomic units'
    label2 = 'Theoretical radial probability distribution. \n Energy: ' \
           + str(-0.5) + ' atomic units'
    ax_potential.plot(r, phi_min**2, label=label1)
    ax_potential.plot(r, 4*np.pi*r**2*phi_s_H**2, '--', label=label2)
    ax_potential.set_xlabel('Radius [atomic units]')
    ax_potential.set_ylabel('Probability [a.u]')
    ax_potential.legend(loc=1)

    plt.savefig('RadProb.eps')
    plt.savefig('RadProb.png')
    plt.show()

    with open("calculation_outputs.txt", "w") as textfile:
        textfile.write("Minimum energy: " + str(np.real(energy_min)) + "\n")
        textfile.write("Normalization: " + str(np.real(trapz(phi_min**2, r))) + "\n")

if __name__ == '__main__':
    main()
