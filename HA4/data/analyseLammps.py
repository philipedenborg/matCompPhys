import matplotlib.pyplot as plt
import numpy as np

f = open("thermal_exp2.txt","r")
data = []
for line in f:
    string = line.strip().split('\t')
    col = [float(s) for s in string]
    data.append(col)
f.close()    
data = np.array(data)

T = data[:,0]

alpha = data[:,1]
T_plot = T[6::10]
alpha_plot = alpha[6::10]

T_exp = np.array([0, 100, 200, 300, 400, 500, 600, 650])+273
alpha_exp = np.array([22, 25.4, 26.5, 27.8, 29.9, 32.5, 35.5, 37.2])*3


h1 = plt.plot(T,alpha*1e6,'b-',label="Simulation data")
#plt.ylim([55,110])
plt.hold(True)
h2 = plt.plot(T_exp,alpha_exp,'r-s',label="Experimental data")
plt.legend(labels=["Simulation","Experiment"],loc=0)
plt.plot(T_plot,alpha_plot*1e6,'bo')
plt.rc('text',usetex=True)
plt.rc('font', family='serif')
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
plt.xlabel('Target temperature [K]',fontsize=16)
plt.ylabel('Thermal expansion [$10^{-6}$ K$^{-1}$]',fontsize=16) 
plt.grid(True)


plt.show()

