import numpy as np
from math import pi
if __name__ == "__main__":
    
    Rmax=1e9 #cm, =1e4km 
    clight=3e10 #cm/s
    vmax=0.1*clight 
    #vmax=0.01*clight
    Msun=2e33 #g
    Mtotal=0.01*Msun
    rho=Mtotal/(4/3*pi*Rmax**3)
    T=1e9 #K
    n_slices=3000
    Mass_coordinates=np.linspace(0.0,Mtotal,n_slices)
    R_coordinates=(Mass_coordinates/Mtotal)**(1/3)*Rmax
    
    V_of_R=R_coordinates*vmax/Rmax
    '''
    print('%.9e' % R_coordinates[n_slices-1])
    print(Mass_coordinates[n_slices-1])
    print(V_of_R[n_slices-1])
    '''
    
    with open("Homologous_expansion.dat",'w') as outputfile1:
        outputfile1.write(str(len(R_coordinates))+'\n')
        for i in range(len(R_coordinates)):
            s0=str(i+1).rjust(6," ")
            s1="     %.9e" %Mass_coordinates[i]
            s2="     %.9e" %R_coordinates[i]
            s3="     %.9e" %T
            s4="     %.9e" %rho
            s5="     %.9e" %V_of_R[i]
            s6="     %.9e" %1.0
            s7="     %.9e" %0.0
            
            outputfile1.write(s0+s1+s2+s3+s4+s5+s6+s7+'\n')
            
    with open("Homologous_expansion_composition.dat",'w') as outputfile2:
        outputfile2.write(str(len(R_coordinates))+'  '+'1\n')
        outputfile2.write('1.0d0\n1.0d0\n')
        for i in range(len(R_coordinates)):
            s1="     %.9e" %Mass_coordinates[i]
            s2="     %.9e" %R_coordinates[i]
            s3="     %.9e" %1.0
            outputfile2.write(s1+s2+s3+'\n')
        
        
            
            
       
    

    
