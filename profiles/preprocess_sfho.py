import numpy as np
from math import pi
from matplotlib import pyplot as plt
#filename = "sfho.dat"
filename = "DD2_M13641364_M0_LK_SR_R04.dat"

with open(filename) as original_file:
    quantity_info = original_file.readline()
    data = original_file.readlines()

Mass = np.zeros(len(data))
Entropy = np.zeros(len(data))
Rho = np.zeros(len(data))
Temperature = np.zeros(len(data))
Ye = np.zeros(len(data))
Velocity = np.zeros(len(data))
Tau=np.zeros(len(data))

Radius = np.zeros(len(data))



for i in range(len(data)):
    line = data[i].split()
    Mass[len(data)-1-i] = float(line[0])
    Entropy[len(data)-1-i] = float(line[1])
    Rho[len(data)-1-i] = float(line[2])
    Temperature[len(data)-1-i] = float(line[3])
    Ye[len(data)-1-i] = float(line[4])
    Velocity[len(data)-1-i] = float(line[5])
    Tau[len(data)-1-i] = float(line[6])

r_inner = 2.95e7  #inner radius = 295 km, the point to generate ejecta profiles in whiskyTHC simulations

Radius[0] = r_inner

#print("max entropy: ", max(Entropy))
#print("min entropy: ", min(Entropy))



for i in range(1,len(data)):
    
    dmass = Mass[i]-Mass[i-1]
    r_center = (Radius[i-1]**3 + dmass/2/(4*pi/3*Rho[i-1]))**(1/3)
    Radius[i] = (r_center**3 + dmass/2/(4*pi/3*Rho[i]))**(1/3)
    

#with open("modified_sfho.dat",'w') as outfile:
with open("DD2.dat",'w') as outfile:
    outfile.write("index Mass[g] Radius[cm] Temperature[K] Rho[g/cm^3] Velocity[cm/s] Ye[-] Entropy[kb/baryon] Tau[s]\n")
    outfile.write(str(len(Mass))+'\n')
    for i in range(len(Mass)):
        s0=str(i+1).rjust(6," ")
        s1="     %.9e" %Mass[i]
        s2="     %.9e" %Radius[i] 
        s3="     %.9e" %Temperature[i]
        s4="     %.9e" %Rho[i]
        s5="     %.9e" %Velocity[i]
        s6="     %.9e" %Ye[i]
        s7="     %.9e" %Entropy[i]
        s8="     %.9e" %Tau[i]        
        outfile.write(s0+s1+s2+s3+s4+s5+s6+s7+s8+'\n')

#plt.plot(Mass,Velocity)
#plt.show()











