from math import *
import numpy as np 
import matplotlib.pyplot as plt

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
    
    #CTIO_Andicam BVRIJHK
    
    Andicamlist = ["CTIO_ANDICAM.B_YALO.dat","CTIO_ANDICAM.V_KPNO.dat",
                "CTIO_ANDICAM.R_KPNO.dat","CTIO_ANDICAM.I_KPNO.dat",
                "CTIO_ANDICAM.J.dat","CTIO_ANDICAM.H.dat","CTIO_ANDICAM.K.dat"]

    Wavelength_list, Transmission_list,min_wavelength,max_wavelength = preprocess(Andicamlist)
    
    with open("CTIO_Andicam_BVRIJHK.dat",'w') as outfile1:
        outfile1.write(website+'\n')
        
        s = ""
        for band in Andicamlist:
            s += band+'  '
        s += '\n'
        outfile1.write(s)
        
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
        outfile1.write("(wavelength and transmission)*7 bands:\n")
        outfile1.write(str(len(Wavelength_list[0]))+'\n')
        
        for j in range(len(Wavelength_list[0])):
            s = ""
            for i in range(len(Wavelength_list)):
                s += "  %.9e" %Wavelength_list[i][j]
                s += "  %.9e" %Transmission_list[i][j]
            s += '\n'
            outfile1.write(s)
                
    
              
    plt.figure()
    plt.plot(Wavelength_list[0]/angstrom, Transmission_list[0],color='blue',alpha=0.8,label='CTIO_ANDICAM.B_YALO')
    plt.plot(Wavelength_list[1]/angstrom, Transmission_list[1],color='green',alpha=0.8,label='CTIO_ANDICAM.V_KPNO')  
    plt.plot(Wavelength_list[2]/angstrom, Transmission_list[2],color='orange',alpha=0.8,label='CTIO_ANDICAM.R_KPNO')
    plt.plot(Wavelength_list[3]/angstrom, Transmission_list[3],color='darkorange',alpha=0.8,label='CTIO_ANDICAM.I_KPNO')  
    plt.plot(Wavelength_list[4]/angstrom, Transmission_list[4],color='tomato',alpha=0.8,label='CTIO_ANDICAM.J')  
    plt.plot(Wavelength_list[5]/angstrom, Transmission_list[5],color='red',alpha=0.8,label='CTIO_ANDICAM.H')  
    plt.plot(Wavelength_list[6]/angstrom, Transmission_list[6],color='indianred',alpha=0.8,label='CTIO_ANDICAM.K')  
    plt.xlabel('Wavelength[angstrom]')
    plt.ylabel('Transmission')
    plt.xlim([0,25000])
    plt.legend()
    plt.savefig('CTIO_Andicam_BVRIJHK.png')  



      
    
     #CTIO_DECam     
            
    DECamlist = ["CTIO_DECam.u_filter.dat","CTIO_DECam.g_filter.dat","CTIO_DECam.r_filter.dat",
                 "CTIO_DECam.i_filter.dat","CTIO_DECam.z_filter.dat","CTIO_DECam.Y_filter.dat"]
    Wavelength_list, Transmission_list,min_wavelength,max_wavelength = preprocess(DECamlist)
    
    with open("CTIO_DECam_ugrizY.dat",'w') as outfile2:
        outfile2.write(website+'\n')
      
        s = ""
        for band in DECamlist:
            s += band+'  '
        s += '\n'
        outfile2.write(s)
        
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
        outfile2.write("(wavelength and transmission)*6 bands:\n")
        outfile2.write(str(len(Wavelength_list[0]))+'\n')
        
        for j in range(len(Wavelength_list[0])):
            s = ""
            for i in range(len(Wavelength_list)):
                s += "  %.9e" %Wavelength_list[i][j]
                s += "  %.9e" %Transmission_list[i][j]
            s += '\n'
            outfile2.write(s)
    
    
    plt.figure()
    plt.scatter(Wavelength_list[0]/angstrom, Transmission_list[0],color='blue',alpha=0.8,label='CTIO_DECam.u_filter')
    plt.scatter(Wavelength_list[1]/angstrom, Transmission_list[1],color='green',alpha=0.8,label='CTIO_DECam.g_filter')  
    plt.scatter(Wavelength_list[2]/angstrom, Transmission_list[2],color='orange',alpha=0.8,label='CTIO_DECam.r_filter')
    plt.scatter(Wavelength_list[3]/angstrom, Transmission_list[3],color='darkorange',alpha=0.8,label='CTIO_DECam.i_filter')
    plt.scatter(Wavelength_list[4]/angstrom, Transmission_list[4],color='tomato',alpha=0.8,label='CTIO_DECam.z_filter')  
    plt.scatter(Wavelength_list[5]/angstrom, Transmission_list[5],color='red',alpha=0.8,label='CTIO_DECam.Y_filter')
    
    plt.xlabel('Wavelength[angstrom]')
    plt.ylabel('Transmission')
    plt.xlim([0,12000])
    plt.legend()
   # plt.savefig('CTIO_DECam_ugrizY.png') 

