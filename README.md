# kNEC 

Zhenyu Wu 

171840687@smail.nju.edu.cn  or  zhenyuwu99@gmail.com


kNEC (KiloNova Explosion Code) is a fortran code that simulates hydrodynamical evolution of BNS merger ejecta and the corresponding kilonova emission. It is based on the SNEC code by Morozova et al., which is a Lagrangian radiation-hydrodynamics code for core-collapse supernova explosion. The current KNEC code is stable and contains the physical processes necessary for kilonova simulation. However, it's not ready for publication because there are some experimental features. Moreover, I keep some features for supernovae now, although they are useless to kilonovae.

If you are not familiar with SNEC, I recommend you to study SNEC first, especially the snec notes.
https://stellarcollapse.org/SNEC.html
https://stellarcollapse.org/codes/snec_notes-1.00.pdf

kNEC code has 2 branches: master and knec-parallel. The master branch, like the original SNEC code, processes only one parameters file at a time, and output part of the information to the screen when running. The knec-parallel branch has a shell script multiknec.sh in it, and src/input_parser.F90 is changed a little bit. With multiknec.sh, this version of the code batches multiple parameters files. Note that kNEC does not involve parallel computing.

## Install and Run

### prerequisite

gfortran  LAPACK

### master branch

* change make.inc according to your machine
* make clean
* make  (this produces an executable *snec* that can be run )
* change parameters file according to the problem. If in the parameters file outdir = "Data", create a folder Data. 
* ./snec

### knec-parallel branch

* change make.inc according to your machine

* mkdir parameter_list,  put all the parameters files in it. Parameters files should be named in accordance with outdir. For example, if there are two parameters files in parameter_list: blh_Apr2_zb, sfho_Apr2_zb, then in blh_Apr2_zb outdir should be "Data_blh_Apr2_zb", and in sfho_Apr2_zb outdir should be "Data_sfho_Apr2_zb". There are some examples in the folder parameter_examples.

* sh multiknec.sh. When the program is finished, the screen_output files will be placed in the Data_ folders.

    

## Main output

bolometric light curves:  Data/lum_observed.dat

multicolor light curves: Data/Magnitude_KNEC.dat

see src/output.F90, src/analysis.F90, src/shock_capture.F90, src/conservation.F90 for more information.





## Parameter Settings

example: blh_Apr2_zb

```txt
#____________LAUNCH_____________
outdir              = "Data_blh_Apr2_zb"      #output data dir

#___________PROFILE_____________
profile_name = "profiles/modified_blh.dat"     # blh profile
read_composition_switch = 0         # We don't calculate compositions in kNEC
comp_profile_name	= ""              # We don't calculate compositions in kNEC

#__________EXPLOSION_____________
initial_data 		= "Thermal_Bomb" 
# "Thermal Bomb" + bomb_mode 2 + final_energy=0 means no explosion (explosion is not needed in kNEC unless you want to inject a shock to study jet-ejecta interactions, in which case bomb or piston parameters matter.)
#Options:
#"Piston_Explosion"
#"Thermal_Bomb"

piston_vel          = 0.6d10      
piston_tstart       = 0.0d0
piston_tend         = 1.0d-2
final_energy        = 0.0d0
bomb_tstart         = 0.0d0
bomb_tend           = 0.0d0
bomb_mass_spread    = 0.0d0 #(in solar mass)
bomb_start_point    = 1
bomb_mode           = 2
#_____________GRID_______________

imax         = 1000              # grid the same as SNEC
gridding = "from_file_by_mass"

#Options:
#"uniform_by_mass"
#"from_file_by_mass"

mass_excision = 0  # SNEC only, ignore this in kNEC
mass_excised = 1.4 #in solar mass, SNEC only, ignore this in kNEC
mass_gravity_switch = 1   # enable gravity
mass_extragravity = 3.0d0 #in solar mass, represents gravitational pull from central engine, 3 Msun by default for kNEC
read_inner_radius_switch = 1   # read inner radius by default for kNEC, in fact this feature is redundant considering the ejecta profile, but we keep it for now.
inner_radius = 2.95d7   # set inner radius to 295km for realistic profiles like blh, sfho, and DD2 because they are extracted from numerical relativistic simulations at 295km. 
#inner_radius = 0.25d9  # set inner radius to 0.25d9 cm for our wind profiles.

#___________EVOLUTION_____________

radiation = 1  # enable radiation 
eoskey = 2   # Paczynski EOS for kNEC
#eoskey Options:
#1 - ideal eos
#2 - Paczynski

continuous_boundary_switch = 0 
#Options for continuous_boundary_switch: 0 means p(imax)=0, or zero boundary condition (zb). 1 means p(imax)=p(imax-1), or continuous boundary condition(cb). You should set continuous_boundary_switch = 0 for kNEC simulations.

Ni_switch = 0    #SNEC only, ignore Ni heating in kNEC, set Ni_switch = 0
Ni_by_hand = 0
Ni_mass = 0.05 			#(in solar mass)
Ni_boundary_mass = 3.0d0    #(in solar mass, should be larger than excised mass)
			    #(attention - smoothing is going to change it, if applied)
Ni_period = 1.0d4

saha_ncomps = 0   # kNEC doesn't solve Saha equations, set saha_ncomps = 0
mu = 100.0d0 # mean molecular weight after nucleosynthesis 
ybar = 2.0d0 # average number of ionized electrons per atom
ye_afternucleosynthesis = 0.4d0
#notice: this ye is ye after nucleosynthesis, used for hydrodynamics; ye pre nucleosynthesis is provided by ejecta profile and used for heating rates, opacities. 
#Through experiments, we find that mu and ybar have little impact on light curves. We set them 100 and 2 respectively by default and the corresponding ye_afternucleosynthesis should be ~ 0.4.

boxcar_smoothing = 0  # In kNEC, we don't need boxcar to simulate mixture of elements due to explosion. 

opacity_floor_envelope = 0.01d0   #SNEC only, ignore this in KNEC
opacity_floor_core     = 0.24d0   #SNEC only, ignore this in KNEC

#nuclear heating rates 

heating_formula = "Apr2"   
# options: (1)"Korobkin" (2)"LR15" (3)"Ricigliano" (4)"Hybrid"(5)"arctan"(6)"Apr2"
# (1) ~ (5) are experimental heating rates. Please use (6) "Apr2", which uses atan + powerlaw with a transition to fit results from SkyNet nuclear reaction network. 

heating_epsilon_th = 0.5d0 # (thermalization efficiency, the proportion of energy deposition. It should be a function of time, but for simplicity we keep it constant in kNEC. 0.5 by default)

# The following parameters are for option(1) "Korobkin" and (4)
heating_epsilon_0 = 2.0d18
heating_sigma = 0.11d0
heating_t_0 = 1.3d0
heating_alpha = 1.3d0

heating_t_cutoff = 10.0d0

#(5) "arctan"
#atan + powerlaw without transition
arctan_t2 = 0.1d0 # day
arctan_eps0 = 2.0d18   #2.0d19
arctan_sigma0 = 0.11d0   #0.03
arctan_t0 = 1.3d0 

#(6)"Apr2"
#atan+powerlaw with transition,  no extra parameters needed here

#___________TIMING_______________  
ntmax               = 10000000000000
tend                = 3.0d6   # simulation ends at ~ 35 day after merger.
dtout               = 1.0d4
dtout_scalar         = 1.0d3
dtout_check          = 1.0d3
#dtout_scalar        = 1.7d3
#dtout_check         = 1.7d3

ntout               = 10000
ntout_scalar        = 10000
ntout_check         = 10000

ntinfo              = 1000

dtmin               = 1.0d-10
dtmax               = 1.0d2

#____________TEST_________________

sedov 		    = 0

```







## From SNEC to kNEC

### 1. Initial Setup

Basic hydrodynamical equations for KNEC are the same as SNEC. The Nickle heating term is replaced by radioactive heating of r-process elements ("simple_heating" in the program). The opacity is replaced by a simple function of initial Ye.

Grid setup is the same as SNEC. We find that grid setup has little impact on final results.

There are 2 types of explosions in SNEC: thermal bomb and piston explosion. However, we don't have explosions in KNEC unless we need to inject a shock. So we set explosion type = thermal bomb and input energy = 0. Boxcar smoothing is turned off in KNEC. 

#### 1.1 Ejecta profiles

Ejecta profiles format:

index     mass[g]      radius[cm]      temperature[K]      density[g/cm^3]      velocity[cm/s]      ye(initial)      entropy[kb/baryon]      expansion_timescale[s]

In fact, the radius data is not used, but recalculated from mass and density in KNEC. The radius of inner boundary is provided in the parameters file.

r-process nucleosynthesis happens in the first seconds. The ye here is initial ye, or ye before nucleosynthesis ("ye_initial" in the code). It is used for opacities and heating rates, playing an important role in determinating kilonova properties. On the other hand, ye after nucleosynthesis ("ye" in the code, "ye_afternucleosynthesis" in parameters file) is provided in the parameters file. It affects EOS and hydrodynamical evoulution. We set this ye to 0.4 and find it has little impacts on final results.

Currently We have 2 kinds of ejecta profiles for KNEC. 

1. wind3 profiles (e.g. wind3_0.01M_0.2c_ye0.05_s10_tau10.dat)

   (see 3. MODELS AND SETUP in Tanaka's paper <https://iopscience.iop.org/article/10.1088/0004-637X/775/2/113/pdf> )

   wind3 means $\rho \propto 1/r^3$.  

   For homologous expansion, $v \propto r$, and velocity of inner boundary is set to 0.05c. 

   Temperature is 1e9K everywhere.

   To reduce the outer boundary effect due to large pressure gradient, and to produce light curves of homologous expansion, we also designed the wind310T6 profile, which means density-radius are 2 powerlaws (powerlaw index 3 and 10), and temperature decreases as a powerlaw near outer boundary (index 6). The profiles are generated by wind.ipynb in KNEC/profiles and you can change the parameters, such as ejecta mass, initial ye, maximum velocity, entropy, tau(expansion timescale), etc.  

   You can also design your own profiles like uniform ejecta. 

   

2. profiles from David Radice's numerical relativity simulations 

   (modified_sfho.dat, modified_blh.dat  and DD2.dat)   

* ignore profiles/ 15Msol_RSG\, stripped_star\, sedov, which are only used in SNEC



#### 1.2 Boundary conditions

We have 2 kinds of boundary conditions ("continuous_boundary_switch" in parameters file):

(1) p(imax) = 0, which is the same as SNEC. You should use this boundary condition in kNEC.

This boundary condition may lead to velocity divergence at outer boundary (sometimes even exceeds speed of light) But we find that this boundary effect has little impact on final results,  because the part whose velocity diverges has so little mass that it is negligible.  

(2) p(imax) = p(imax-1)   not recommended for kNEC

The velocity of ejecta will remain roughly constant, but we find that the condition leads to a large proportion (~ 5%) of pdV term at boundaries in energy conservation. Since we don't model the gas around the ejecta, we hope the influences of boundaries are negligible in energy conservation. Moreover, p(imax) = 0 results match better with homologous expansion. Therefore, we do not recommend this boundary condition for kNEC.

### 2. EOS

kNEC calculates radiation and hydrodynamics of the ejecta according to its key properties like initial ye, ignoring the details of composition. "read_composition_swtich" in parameters file  should be set to 0, and composition profile is not needed.

We use a simplified Paczynski EOS in kNEC. Saha equations are not solved and the correction terms for ionization in Paczynski EOS are ignored. We provide a free parameter "ybar" in the parameters file as the mean degree of ionization. We find this parameter has little impacts on final results. 

 

### 3. Heating rates and opacities

We model opacity as a simple function of initial Ye (opacity_simple.F90): 

$\kappa = 1+ \frac{9}{1+(4Ye)^{12}} \mathrm{[cm^2/g]}$

The maximum opacity is 10 cm^2/g and the minimum is 1 cm^2/g. At Ye = 0.25, opacity = 5.5 cm^2/g. The factor 12 in this formula makes opacity drops steeply near Ye = 0.25.

We have 6 heating rate models in kNEC.  (6) Apr2 models are recommended, the others are experimental. For all of these models, there is an extra parameter for thermalization efficiency, which means the fraction of energy deposited in the ejecta.   ("heating_epsilon_th" in parameters file). For simplicity, we set it a constant 0.5, although it may vary with time in more accurate models. 

(1) Korobkin's 

https://iopscience.iop.org/article/10.3847/2041-8213/aa9ab9   

This formula is only valid for low-Ye low-entropy ejecta.

(2) LR 15  

The heating rate table is obtained from Albino Perego's nuclear network simulations. This model calculates coefficients at every step so it costs much time.

(3) Ricigliano's  (single power-law)

The heating rate table (tables/epsdatafit.dat) is obtained from Albino Perego and Giacomo Ricigliano's nuclear network simulations, providing coefficients of a single powerlaw. But the heating rate goes to infinity at very early times.

(4) Hybrid 

An experimental model combining (1) and (3)

(5) arctan 

Also uses tables/epsdatafit.dat. Combine an arctan function (Korobkin-like) and a single powerlaw to fit Perego's simulations. There are 4 free parameters in this model: t2, eps0, sigma0 and t0. 

Before t2, the heating rate is an arctan function, while after t2 is a single powerlaw. The heating rate is almost flat around eps0 when time is very small. t0 determines the extension of the plateau, sigma0 determines how much heating rate drops after the plateau and before the powerlaw. (see heating_rate_arctan_module.F90 for details)

typical parameters:   

 t2 = 0.1day, eps0 = 2e18 erg/g/s, sigma0 = 0.11, t0 = 1.3 

 t2 = 0.1day, eps0 = 2e19 erg/g/s, sigma0 = 0.03, t0 = 1.3 

(6)  Apr2 (recommended)

(6) = (5) + a transition region between the arctan and powerlaw. Moreover, (6) uses 2 coefficient tables and there are no free parameters in the model, unless you change the code mannually, or change initial Ye/entropy/expansion timescale.

### 4. Output

The main outputs are bolometric luminosities and AB apparent magnitudes in certain bands. 

bolometric luminosity: 

Bolometric luminosity = luminosity at photosphere + radioactive heating above the photosphere, similar to SNEC(luminosity at photosphere + Nickel heating above) . We do not have a temperature floor like SNEC. 

AB apparent magnitudes in certain bands:

We assume black body radiation at the photosphere and black body radiation at each shell above the photosphere. Note that temperatures of these black bodies can be different. The formula is documented in https://www.overleaf.com/read/kwvymfznrkmc. 

We get filters' data from http://svo2.cab.inta-csic.es/svo/theory/fps/index.php?mode=browse&gname=SLOAN&asttype=   The distance between kilonova event and observer is 40 Mpc, the same as GW170817 event.

The very early light curves (<0.01day) may be unreliable because of lack of resolution at the outer boundary. Late time light curves (>10day) may also be unreliable due to deviation from black body models.

















