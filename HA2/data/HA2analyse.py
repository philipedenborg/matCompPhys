import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np

res = ''
plt.rc('text',usetex=True)
plt.rc('font', family='serif')

while res != '0':
    print "K-points = kp, Cut off  = cf, Equation of state = eos, Surface energies = surf. Write 0 to exit"
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
        plt.plot(cf[:,0],cf[:,1])
        plt.rc('text',usetex=True)
        plt.rc('font', family='serif') 
        plt.xlabel('Cut-off energy [eV]',fontsize=10)
        plt.ylabel('Energy [eV]',fontsize=10) 
        plt.show()
    elif res == 'eos':
        V = 16.529
        a = (4*V)**(1/3.0)
        B = 77.708
        print "Lattice constant: %.3f Å" % a
        print "Bulkmodulus: B %.3f GPa" % B  
        img=mpimg.imread('Al-eos.png')
        imgplot = plt.imshow(img)
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
        plt.show()
    elif res == 'surf':
        
        
        # fcc(100)
        A = 1
        E_slab = 1
        N_slab = 1
        N_bulk = 1
        E_bulk = 1

        sigma100 = 1/(2*A)*(E_slab - N_slab/N_bulk * E_bulk)
        # fcc(111)
        A = 1
        E_slab = 1
        N_slab = 1
        N_bulk = 1
        E_bulk = 1
         
        sigma111 = 1/(2*A)*(E_slab - N_slab/N_bulk * E_bulk)

        
        f = open('k_file','r')
        data = []
        for line in f:
            string = line.split(' ')
            col = [float(s) for s in string]
            data.append(col)
        data = np.array(data)
        k = data[:,0]
        E = data[:,1]
        plt.figure(1)
        plt.plot(k,E,'-o')
        plt.rc('text',usetex=True)
        plt.rc('font', family='serif')
        plt.xlabel('Number of k-points',fontsize=16)
        plt.ylabel('Energy [eV]',fontsize=16)
        plt.show()
        


