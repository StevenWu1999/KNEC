## 1. Introduction

## 2. Methods

### 2.1 Breif overview of SNEC

### 2.2 Opacities

* opacity as a function of Ye, compare kNEC and Tanaka et al.

### 2.3 Heating rates

* Heating rate trajectories as obtained by SkyNet
* Heating rate fitted trajectories 
* SkyNet and fitted, 5 different representative sets of thermodynamic variables

### 2.4 Initial and boundary conditions

#### Initial ejecta profile

* SFHo/BLh/DD2: Ye, density, velocity as a function of enclosed mass (2+ figures)
* uniform/wind3/wind310/wind310T6(optimal wind): density, velocity, temperature as a function of radius (2+ figures)
* Question
  1. ignore wind-ex profiles?  
  2. compare hydro and homologous expansion here to explain why we use wind310, or leave this in code validation? 

#### Boundary condition p(imax) = 0 and its problem

* final velocity as a function of enclosed mass, wind310Tx profiles, no heating

* final velocity as a function of enclosed mass, wind310Tx profiles, with heating (thermalization efficiency 0.5)

* final velocity as a function of enclosed mass, blh profile, no heating

* Question 

  1. Do we need to mention that blh (with heating) boundary velocity exceeds speed of light?

     (Comparison with blh-with-modified-vel profile shows that this makes little difference to light curves. For the part of ejecta whose vel > 0.6c, mass and kinetic energy are negligible compared to the whole ejecta.)

### 2.5 Other difference with SNEC

* Lbol: effects of mean molecular weight $\mu$, mean degree of ionization $\bar{y}$

### 2.6 Bolometric luminosities and multicolor luminosities



## 3. Code Validation

### 3.1 Hydrodynamics

* Lbol: compare hydro and homologous expansion, using wind/optimal wind profiles
* Lbol: compare hydro and homologous expansion (v = kr ?) , using blh-with-modified-vel profile 
* velocity fit: blh and blh-with-modified-vel 
* AB Mags: compare hydro and homologous expansion (?) , using blh-with-modified-vel profile
* density of ejecta: compare hydro and homologous expansion

### 3.2 Energy conservation

* E1 = Eejecta = Ekin + Egrav + Eint; E2 = Einitial +  Eheating - Eradiation + boundary pdV work

  E1, E2 as a function of time; optimal wind/blh profiles

### 3.3 Comparison with analytic models

* Lbol: kNEC VS Ricigliano's model, uniform/optimal wind, Ye = 0.1 ~ 0.4
* Photospheric radius Rph, AB mags, Teff
* kNEC VS Rahul's model

## 4. Ab-initio simulations: from mergers to kilonovae

### 4.1 General Features

* ?    Lbol, AB Mags, Energy, Rph, Teff

### 4.2 Impact of uncertainties in heating rates

* Lbol/heating rate: wind/optimal wind, heating x3, x3 (0-10 sec), x0.3, x0.3 (0-10 sec), high Ye 0.4, low Ye0.1

* Lbol/AB Mags: blh profile, heating x3, x3 (0-10 sec), x0.3, x0.3 (0-10 sec)    

  ? combine this with piston shock results

## 5. A first application to AT2017GFO (Combine this with 4?)

### 5.2 Comparison between NR informed models and observations

* blh/sfho extrapolation method: density/mass as a function of time
* AB Mags: blh and extrapolated blh

### 5.3 Impact of shock injection

















###  

















