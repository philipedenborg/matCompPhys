from subprocess import call
import numpy as np

DT = 10
temperatures = range(200,900+2*DT,DT)


for temp in temperatures:
    f = open("equilibrate_solid_minus_temp.in", "r")
    in_file = open("datafiles/in_files/"+str(temp)+".in","w")
    old = f.read() # read everything in the file
    in_file.write("variable\tT equal "+str(temp)+" \n" + old) # write the new line before 
    
    
 
    f2 = open("submitfile_minus_file","r")
    sub_file = open("submitfile","r+")
    old_sub = f2.read()
    sub_file.seek(0) # rewind
    sub_file.write("#!/usr/bin/env bash\n")
    sub_file.write("#SBATCH -o "+str(temp)+"log1D.lammps\n")
    sub_file.write(old_sub + "\n")
    sub_file.write("mpirun ./lmp_mpi -in datafiles/in_files/"+str(temp)+".in") 
    call(["sbatch", "submitfile"])
    
    f.close()
    f2.close()
    sub_file.close()
    in_file.close()
       
