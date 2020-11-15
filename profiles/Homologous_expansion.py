import numpy as np
from math import pi
if __name__ == "__main__":
    
    Rmax=1e9 #cm 
    clight=3e10 #cm/s
    vmax=0.1*clight 
    #vmax=0.01*clight
    Msun=2e33 #g
    Mtotal=0.01*Msun

    
    n_slices=3000
    Mass_coordinates=np.linspace(Mtotal/1e5,Mtotal,n_slices)
    R_coordinates=(Mass_coordinates/Mtotal)**(1/3)*Rmax
    
    V_of_R=R_coordinates*vmax/Rmax
    T = 1e9
    density_profile = "wind"
    Temperature_smooth = "nosmooth"
    Temperature_smooth_length = 301
    
    #density profile
    if density_profile == "uniform":
        rho=Mtotal/(4/3*pi*Rmax**3)
        Rho=np.ones(len(R_coordinates))*rho
    elif density_profile == "wind":
        #rho = k / r^2
        k = Mtotal/(4*pi*Rmax)
        Rho = k/(R_coordinates**2)

    #temperature smoothing
    if Temperature_smooth == "smooth":     
        Temperature = T*np.ones(len(R_coordinates))
        for i in range(Temperature_smooth_length):
            Temperature[n_slices-1-i] = T/(Temperature_smooth_length-1)*i
        ratio = 1-Mass_coordinates[n_slices-1-Temperature_smooth_length]/Mass_coordinates[n_slices-1]
        print('temperature is smoothed for '+str(ratio*100)+r'% of the mass')
    else:
        Temperature = T*np.ones(len(R_coordinates))
    

    

    '''
    print('%.9e' % R_coordinates[n_slices-1])
    print(Mass_coordinates[n_slices-1])
    print(V_of_R[n_slices-1])
    '''
    
    with open("Homologous_expansion_wind.dat",'w') as outputfile1:
        outputfile1.write(str(len(R_coordinates))+'\n')
        for i in range(len(R_coordinates)):
            s0=str(i+1).rjust(6," ")
            s1="     %.9e" %Mass_coordinates[i]
            s2="     %.9e" %R_coordinates[i]
            s3="     %.9e" %Temperature[i]
            s4="     %.9e" %Rho[i]
            s5="     %.9e" %V_of_R[i]
            s6="     %.9e" %0.3
            s7="     %.9e" %0.0
            
            outputfile1.write(s0+s1+s2+s3+s4+s5+s6+s7+'\n')
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
            
            
       
    

    
