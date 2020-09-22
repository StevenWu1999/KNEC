## KNEC progress



### Sep 21 meeting version: opacity and heating rate from Perego+2017

Opacity_simple.F90 is created to decide the opacity of the ejecta as a function of Ye. Currently the function is discontinuous, with opacity either being 10cm^2/g or 1cm^2/g. 

nuclear_heating_rate.F90 is created to calculate the heating rate from Perego+2017. An exponential factor \alpha=1.3 should be added later to correct the typo. Since the heat is not mainly from gamma rays,  radiation transfer equation is not solved.

### Initial commit--Homologous expansion

The differences between this initial commitment and the original SNEC1.01 files are listed in Homologous expansion--changes to the original files.md. It simulates the homologous expansion of pure hydrogen gas, using the ideal gas EOS. The total mass is 2.0e31 g, the maximum radius is 1.0e12 cm. The initial velocity is proportional to radius, and vmax=0.1c. 
Initial temperature is 10^4 K. These are described in profiles/Homologous_expansion.dat and profiles/Homologous_expansion_composition.dat.  Note that the last two columns of Homologous_expansion.dat do not matter because SNEC doesn't read them. 

There is no bomb energy input. Opacity, luminosity and heating rate are not taken into consideration.  hydro.F90 is used to solve the motion of the fluid.



