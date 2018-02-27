import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from ase.io import read

res = ''
plt.rc('text',usetex=True)
plt.rc('font', family='serif')

while res != '0':
    print "K-points = kp, Cut off  = cf, Equation of state = eos, Surface energies = surf, E_ads convergence = adsconv. Write 0 to exit"
    res = raw_input("What are you looking at today?\n")

    if res == 'cf':
        f = open('cutoff_file','r')
        cf = []
        for line in f:
            string = line.split(' ')
            col = [float(s) for s in string]
            cf.append(col)
        cf = np.array(cf)
        plt.figure()
        plt.plot(cf[:,0],cf[:,1],'-o')
        plt.rc('text',usetex=True)
        plt.rc('font', family='serif') 
        plt.xlabel('Cut-off energy [eV]',fontsize=16)
        plt.ylabel('Energy [eV]',fontsize=16) 
        ax = plt.gca()
	ax.ticklabel_format(useOffset=False)
        plt.grid(True)		
        plt.show()

    elif res == 'eos':
        V = 16.529
        a = (4*V)**(1/3.0)
        B = 77.708
        print "Lattice constant: %.5f Å" % a
        print "Bulkmodulus: B %.3f GPa" % B  
        img=mpimg.imread('Al-eos.png')
        #imgplot = plt.imshow(img)
        #plt.show()
        
        configs = read('Al.traj@0:5')  # read 5 configurations
	# Extract volumes and energies:
	volumes = [al.get_volume() for al in configs]
	energies = [al.get_potential_energy() for al in configs]
        v = np.array(volumes)
	E = np.array(energies)
	degree = 4 # Degree of polynomial fit
	p = np.polyfit(v,E,degree)
	x = np.linspace(v[0],v[4],200)
	y = np.polyval(p,x)

	plt.figure()
	plt.plot(x,y)
	plt.hold(True)
	plt.plot(v,E,'bo')
	plt.xlabel('Volume [\AA]',fontsize=16)
        plt.ylabel('Energy [eV]',fontsize=16)
        plt.grid(True)
	plt.show()
    elif res == 'kp':
        f = open('k_file','r')
        data = []
        for line in f:
            string = line.split(' ')
            col = [float(s) for s in string]
            data.append(col)
        data = np.array(data)
        k = data[:,0]
        E = data[:,1]
        plt.figure()
        plt.plot(k,E,'-o')
        plt.rc('text',usetex=True)
        plt.rc('font', family='serif')
        plt.xlabel('Number of k-points',fontsize=16)
        plt.ylabel('Energy [eV]',fontsize=16)
        plt.grid(True)
        plt.show()
    elif res == 'adsconv':    
        f = open('conv100_OT_O_energy','r')
        data = []
        for line in f:
            string = line.split(' ')
            col = [float(s) for s in string]
            data.append(col)
        data = np.array(data)
        N = data[:,0]
        E = data[:,1]
        plt.figure(1)
        plt.plot(N,E,'-o')
        plt.rc('text',usetex=True)
        plt.rc('font', family='serif')
        plt.xlabel('Number of cell layers',fontsize=16)
        plt.ylabel('Adsorption energy [eV]',fontsize=16)
        ax = plt.gca()
	ax.ticklabel_format(useOffset=False)
	plt.grid(True)
        plt.show()
        print E[2]
    elif res == 'surf':



