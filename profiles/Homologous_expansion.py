import numpy as np
from math import pi,log
if __name__ == "__main__":
    
    density_profile = "wind-3"
    T = 1e9 
    Msun=2e33 #g
    Mtotal=0.01*Msun
    clight=3e10 #cm/s
    ye = 0.4
    entropy = 10.0 #[kb/baryon]
    tau = 10.0e-3 #[s] 
    

    #ye 0.4,0.2775,0.05
    #density profile
    if density_profile == "uniform":
        Rmax=1e9 #cm 
        vmax=0.1*clight 
        #vmax=0.01*clight       
        n_slices=3000
        Mass_coordinates=np.linspace(Mtotal/1e5,Mtotal,n_slices)
        R_coordinates=(Mass_coordinates/Mtotal)**(1/3)*Rmax
        
        V_of_R=R_coordinates*vmax/Rmax
        rho=Mtotal/(4/3*pi*Rmax**3)
        Rho=np.ones(len(R_coordinates))*rho
    elif density_profile == "wind-2":
        #rho = k / r^2
        Rmax=1e9 #cm 
        vmax=0.2*clight 
        #vmax=0.01*clight
        
        n_slices=3000
        Mass_coordinates=np.linspace(Mtotal/1e5,Mtotal,n_slices)
        R_coordinates=(Mass_coordinates/Mtotal)**(1/3)*Rmax
        
        V_of_R=R_coordinates*vmax/Rmax
        k = Mtotal/(4*pi*Rmax)
        Rho = k/(R_coordinates**2)

    elif density_profile == "wind-3":
        vmax = 0.1*clight
        vmin = 0.05*clight
        Rmax = 1e9 #cm
        Rmin = Rmax*vmin/vmax

        print("Using wind-3 density profile, Rmin = "+str(Rmin)+" cm")
        #rho = k/r^3
        k = Mtotal/(4*pi*log(Rmax/Rmin))

        n_slices=3000
        V_of_R = np.linspace(vmin,vmax,n_slices)
        R_coordinates = V_of_R/vmax*Rmax
        Mass_coordinates=np.zeros(n_slices)
        for i in range(n_slices-1):
            Mass_coordinates[i+1] = Mass_coordinates[i] + 4*pi*k*log(R_coordinates[i+1]/R_coordinates[i])
        
        Rho = k/(R_coordinates**3)

    elif density_profile == "wind-4":
        pass

    Temperature = T*np.ones(len(R_coordinates))
    

    

    '''
    print('%.9e' % R_coordinates[n_slices-1])
    print(Mass_coordinates[n_slices-1])
    print(V_of_R[n_slices-1])
    '''
    
    with open("wind3_0.01M_0.1c_ye0.4_s10_tau10.dat",'w') as outputfile1:
        outputfile1.write('#      index      mass[g]      radius[cm]      temperature[K]      density[g/cm^3]      velocity[cm/s]      ye(initial)      entropy[kb/baryon]      expansion_timescale[s]\n')
        outputfile1.write(str(len(R_coordinates))+'\n')
        for i in range(len(R_coordinates)):
            s0=str(i+1).rjust(6," ")
            s1="     %.9e" %Mass_coordinates[i]
            s2="     %.9e" %R_coordinates[i]
            s3="     %.9e" %Temperature[i]
            s4="     %.9e" %Rho[i]
            s5="     %.9e" %V_of_R[i]
            s6="     %.9e" %ye
            s7="     %.9e" %entropy
            s8="     %.9e" %tau
            
            outputfile1.write(s0+s1+s2+s3+s4+s5+s6+s7+s8+'\n')
    '''        
    with open("Homologous_expansion_composition"".dat",'w') as outputfile2:
        outputfile2.write(str(len(R_coordinates))+'  '+'1\n')
        outputfile2.write('1.0d0\n1.0d0\n')
        for i in range(len(R_coordinates)):
            s1="     %.9e" %Mass_coordinates[i]
            s2="     %.9e" %R_coordinates[i]
            s3="     %.9e" %1.0
            outputfile2.write(s1+s2+s3+'\n')
        
    '''    
            
            
       
    

    
