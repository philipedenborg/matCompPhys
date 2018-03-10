from subprocess import call
import numpy as np

DT=10
temperatures = range(200,900+2*DT,DT)

V = []
T = []
#DT = []#int(temperatures[1])-int(temperatures[0])
alpha = []

files = [str(temps)+"log.lammps" for temps in temperatures]
former_temp = temperatures[0]
for temps in temperatures:
	timesteps = []
	f = open("datafiles/log_files/"+str(temps)+"log1D.lammps","r")
	out = open("datafiles/log_files/"+str(temps)+"lammps.txt","w")	
	data = f.readlines()
	save = False
	for line in data:
		if "Loop time" in line:	
			break	
		if save:		
			string = line.split( )			
			row = [float(s) for s in string]
			timesteps.append(row)
		if "Step Temp" in line:
 				save = True	
				
	for lines in timesteps:
		for elements in lines:
			out.write(str(elements))
			out.write("\t")
		out.write("\n")
	
	timesteps = np.array(timesteps)
	Lx = timesteps[:,8]
	Ly = timesteps[:,9]
	Lz = timesteps[:,10]
	V_t = Lx#*Ly*Lz
	T_t = timesteps[:,1]
	U = timesteps[:,2]
	P = timesteps[:,4]
	
	H = U + P * V_t
	mean_T = np.mean(T_t)
	mean_V = np.mean(V_t)
	mean_H = np.mean(H)
	mean_HV = np.mean(H*V_t)
	
	#k = 8.6173303e-5 # eV/K
	#beta = k*mean_T
	#alpha.append(- k * beta**2 * (mean_HV - mean_H * mean_V)/mean_V)

	V.append(mean_V)

	

	#DT.append(mean_T - former_temp)
	T.append(mean_T)
	#frmer_temp =  mean_T
	f.close()
	out.close()


#V = np.array(V)
#DT = np.array(DT)

output = open("thermal_exp.txt","w")
for i in range(1,len(V)-1):
	alpha = 1/V[i]*(V[i+1] - V[i-1])/(2.0*DT)
	output.write(str(int(T[i]))+"\t"+str(alpha)+"\n")

output.close()

