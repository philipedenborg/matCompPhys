import numpy as np
import numpy.linalg as npl
import matplotlib.pyplot as plt
from scipy.integrate import trapz
def main():
    ''' Finite difference geometry (1D) '''
    r_max = 5 # Maximum radius of position grid in Hartree units
    r_min = 0.01 # Minimum radius of position grid in Hartree units
    n_r = 1000 # Number of elements in position grid
    r = np.linspace(r_min, r_max, n_r) # Position grid in Hartree units

    ''' Differentiation '''
    h = r[1]-r[0] # Step size
    A = create_second_derivative_matrix_1D(r, h)

    ''' Physical constants '''
    Z_helium = 2 # Charge of helium nucleus in hartree units
    Z_hydrogen = 1 # Charge of hydrogen nucleus in hartree units

    ''' Boundary conditions for poisson equation '''
    # U_0 = 0
    # U_inf = 0

    ''' Analytical solutions '''
    V_hartree = 1/r - (1 + 1/r)*np.exp(-2*r) # The hartree potential for hydrogen
    phi_s_H = (1/np.sqrt(np.pi))*np.exp(-r) # Wave function for hydrogen

    ''' Computation '''
    U = compute_VsH_and_U(A, r, phi_s_H)[1]

    ''' Plotting '''
    fig_1 = plt.figure()
    ax_potential = fig_1.add_subplot(111)
    ax_potential.plot(r, U, label='Calculated r*V$_{H}$')
    ax_potential.plot(r, V_hartree*r, '--', label='Theoretical r*V$_{H}$')
    ax_potential.set_xlabel('Radius [a.u.]')
    ax_potential.legend(loc=4)

    plt.show()

    ''' Functions '''
def compute_VsH_and_U(A, r, phi_s_H):
    n_s_H = phi_s_H*phi_s_H # Electron density hydrogen ground state
    # print(trapz(n_s_H, r)) # WHY NOT NORMALIZED???
    U_0 = npl.solve(A, -4*np.pi*r*n_s_H)
    U = U_0 + r/r[-1]
    V_s_H = U/r
    return V_s_H, U

def create_second_derivative_matrix_1D(r, h):
    # http://www.cs.cornell.edu/~bindel/class/cs6210-f12/notes/lec32.pdf
    n = np.size(r)
    A = np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            if j == i:
                A[i, j] = -2
                if i != 0:
                    A[i, j-1] = 1
                if i != n-1:
                    A[i, j+1] = 1
    return A/(h**2)

if __name__ == '__main__':
    main()
