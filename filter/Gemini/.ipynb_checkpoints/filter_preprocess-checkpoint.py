from math import *
import numpy as np 
import matplotlib.pyplot as plt
from scipy import interpolate,integrate

angstrom = 1.0e-8 #cm
nm = 1.0e-7 #cm
clight = 29979245800.0 #cm/s

# https://www.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves


def readdata(filename):
    with open(filename) as filterfile:
        data = filterfile.readlines()
    wavelength = [];transmission = []
    for line in data:
        line = line.split()
        wavelength.append(float(line[0])*angstrom)
        transmission.append(float(line[1]))
    return wavelength,transmission



def preprocess(filename_list):
    
    N_file = len(filename_list)
    Wavelength_list = []
    Transmission_list = []

    for i in range(N_file):
        filename = filename_list[i]
        wavelength, transmission = readdata(filename)
        Wavelength_list.append(wavelength)
        Transmission_list.append(transmission)
    
   
    
    len_wavelength = [len(wavelength) for wavelength in Wavelength_list]
    min_wavelength = [min(wavelength) for wavelength in Wavelength_list]
    max_wavelength = [max(wavelength) for wavelength in Wavelength_list]

    maxlen = max(len_wavelength)
    
    for i in range(N_file):
        wavelength = Wavelength_list[i]
        transmission = Transmission_list[i]
        if len(wavelength) != maxlen:
            dlambda = wavelength[-1] - wavelength[-2]
            for j in range(maxlen - len(wavelength)):
                wavelength.append(wavelength[-1]+dlambda)
                transmission.append(0.0e0)
            Wavelength_list[i] = np.array(wavelength)
            Transmission_list[i] = np.array(transmission)
      
    Wavelength_list = np.array(Wavelength_list) 
    Transmission_list = np.array(Transmission_list)      
    return Wavelength_list,Transmission_list,min_wavelength,max_wavelength        
            
    




if __name__ == "__main__":
    
    website = "http://svo2.cab.inta-csic.es/svo/theory/fps/index.php?mode=browse&gname=SLOAN&asttype="
    
    #GMOS
    
    GMOSlist = ["Gemini_GMOS-N.u.dat","Gemini_GMOS-N.g.dat",
                "Gemini_GMOS-N.r.dat","Gemini_GMOS-N.i.dat","Gemini_GMOS-N.z.dat"]

    Wavelength_list, Transmission_list,min_wavelength,max_wavelength = preprocess(GMOSlist)
    
    with open("Gemini_GMOS-N_ugriz.dat",'w') as outfile1:
        outfile1.write(website+'\n')
        outfile1.write("Gemini_GMOS-N.u  Gemini_GMOS-N.g  Gemini_GMOS-N.r  Gemini_GMOS-N.i  Gemini_GMOS-N.z\n")
        outfile1.write("min(wavelength)/cm and max(wavelength)/cm:\n" )
        s = ""
        for i in range(len(min_wavelength)):
            s += "  %.9e" %min_wavelength[i]
        s += '\n'
        outfile1.write(s)
        
        s = ""
        for i in range(len(max_wavelength)):
            s += "  %.9e" %max_wavelength[i]
        s += '\n'
        outfile1.write(s)
        outfile1.write("(wavelength and transmission)*5 bands:\n")
        outfile1.write(str(len(Wavelength_list[0]))+'\n')
        
        for j in range(len(Wavelength_list[0])):
            s = ""
            for i in range(len(Wavelength_list)):
                s += "  %.9e" %Wavelength_list[i][j]
                s += "  %.9e" %Transmission_list[i][j]
            s += '\n'
            outfile1.write(s)
                
    
              
    plt.figure()
    plt.plot(Wavelength_list[0]/angstrom, Transmission_list[0],color='blue',alpha=0.8,label='Gemini_GMOS-N.u')
    plt.plot(Wavelength_list[1]/angstrom, Transmission_list[1],color='green',alpha=0.8,label='Gemini_GMOS-N.g')  
    plt.plot(Wavelength_list[2]/angstrom, Transmission_list[2],color='orange',alpha=0.8,label='Gemini_GMOS-N.r')
    plt.plot(Wavelength_list[3]/angstrom, Transmission_list[3],color='coral',alpha=0.8,label='Gemini_GMOS-N.i')  
    plt.plot(Wavelength_list[4]/angstrom, Transmission_list[4],color='indianred',alpha=0.8,label='Gemini_GMOS-N.z')  
    plt.xlabel('Wavelength[angstrom]')
    plt.ylabel('Transmission')
    plt.xlim([0,12000])
    plt.legend()
    #plt.savefig('Gemini_GMOS-N_ugriz.png')  

   
    #test
    filter_z = interpolate.interp1d(Wavelength_list[4],Transmission_list[4],'linear')
    options={'limit':100}
    def func(lambda_x):
        global filter_z
        return 3631*1e-23/(6.626e-27)/lambda_x*filter_z(lambda_x)
    
    options={'limit':1000}
    
    I=integrate.nquad(func,[[8.2e-5,1.1e-4]],opts=[options])
    
    print(I[0])
    
    
 
        
    
    
    
    
    
    
    
    
    '''  
    
     #Flamingos2list       
            
    Flamingos2list = ["Gemini_Flamingos2.J.dat","Gemini_Flamingos2.H.dat","Gemini_Flamingos2.Ks.dat"]
    Wavelength_list, Transmission_list,min_wavelength,max_wavelength = preprocess(Flamingos2list)
    
    with open("Gemini_Flamingos2_JHKs.dat",'w') as outfile2:
        outfile2.write(website+'\n')
        outfile2.write("Gemini_Flamingos2.J  Gemini_Flamingos2.H  Gemini_Flamingos2.Ks\n")
        outfile2.write("min(wavelength)/cm and max(wavelength)/cm:\n" )
        s = ""
        for i in range(len(min_wavelength)):
            s += "  %.9e" %min_wavelength[i]
        s += '\n'
        outfile2.write(s)
        
        s = ""
        for i in range(len(max_wavelength)):
            s += "  %.9e" %max_wavelength[i]
        s += '\n'
        outfile2.write(s)
        outfile2.write("(wavelength and transmission)*3 bands:\n")
        outfile2.write(str(len(Wavelength_list[0]))+'\n')
        
        for j in range(len(Wavelength_list[0])):
            s = ""
            for i in range(len(Wavelength_list)):
                s += "  %.9e" %Wavelength_list[i][j]
                s += "  %.9e" %Transmission_list[i][j]
            s += '\n'
            outfile2.write(s)
    
    
    plt.figure()
    plt.plot(Wavelength_list[0]/angstrom, Transmission_list[0],color='gold',alpha=0.8,label='Gemini_Flamingos2.J')
    plt.plot(Wavelength_list[1]/angstrom, Transmission_list[1],color='darkorange',alpha=0.8,label='Gemini_Flamingos2.H')  
    plt.plot(Wavelength_list[2]/angstrom, Transmission_list[2],color='red',alpha=0.8,label='Gemini_Flamingos2.Ks')
    plt.xlabel('Wavelength[angstrom]')
    plt.ylabel('Transmission')
    plt.xlim([10000,30000])
    plt.legend()
    #plt.savefig('Gemini_Flamingos2_JHKs.png') 

    '''