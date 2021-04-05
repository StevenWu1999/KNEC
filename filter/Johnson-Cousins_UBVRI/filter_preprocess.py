from math import *
import numpy as np 
import matplotlib.pyplot as plt

# https://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves

def preprocess(filename):
    
    with open(filename) as filterfile:
        filterfile.readline()
        data = filterfile.readlines()
        
    
    nm = 1.0e-7 #cm
    clight = 29979245800.0 #cm/s
    
    N = 1801

    Freq = np.zeros(N)
    Transmission = np.zeros(N)
    
    for i in range(len(data)):
        if not data[i]=='\n':            
            line = data[i].split()
            wavelength = float(line[0])*nm
            Transmission[i] = float(line[1])/100
            Freq[i] = clight/wavelength
    
    
    return Freq, Transmission
    



if __name__ == "__main__":
    
    Filelist = ["Bessel_U-1.txt","Bessel_B-1.txt",
                "Bessel_V-1.txt","Bessel_R-1.txt","Bessel_I-1.txt"]
    
    
    Transmission_list=[]
    
    N = 1801
    for file in Filelist:
        Freq, Transmission = preprocess(file)
        Transmission_list.append(Transmission)
    
    Transmission_list[4][0]=0.0
    with open("UBVRI.dat",'w') as outfile:
        outfile.write('Johnson-Cousins UBVRI transmission as a function of frequency\n')
        outfile.write(str(N)+'\n')
        for i in range(N):
            f = "     %.9e" %Freq[i]
            U = "     %.9e" %Transmission_list[0][i]
            B = "     %.9e" %Transmission_list[1][i]
            V = "     %.9e" %Transmission_list[2][i]
            R = "     %.9e" %Transmission_list[3][i]
            I = "     %.9e" %Transmission_list[4][i]   
            outfile.write(f+U+B+V+R+I+'\n')
            
    plt.figure()
    plt.fill(Freq, Transmission_list[0],color='purple',alpha=0.3,label='U')
    plt.fill(Freq, Transmission_list[1],color='blue',alpha=0.3,label='B')  
    plt.fill(Freq, Transmission_list[2],color='green',alpha=0.3,label='V')
    plt.fill(Freq, Transmission_list[3],color='red',alpha=0.3,label='R')  
    plt.fill(Freq, Transmission_list[4],color='lightcoral',alpha=0.3,label='I')  
    plt.xlabel('Frequency')
    plt.ylabel('Transmission')
    plt.legend()
    plt.savefig('Johnson-Cousins_UBVRI.png')


    
