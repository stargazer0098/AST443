import numpy as np  
import matplotlib.pyplot as plt 

input_file = 'exoplanets_all_list.txt'

''' 
Input file data structure: 

col 0: name of planet 
col 1: planet radius (jupiter radius)
col 2: orbital period of planet (days)  
col 3: semimajor axis (AU) 
col 4: "epoch of periastron" (JD) 
col 5: impact paramater 
col 6: name of host star 
col 7: right ascension (deg) 
col 8: declination (deg)
col 9: V-band magnitude (apparent) 
col 10: star radius (solar radii) 
''' 
planet_name = np.loadtxt(input_file, dtype = str , skiprows = 1, usecols = 0)
star_name = np.loadtxt(input_file, dtype = str, skiprows = 1, usecols = 6)

print(planet_name)
print(star_name)

data = np.loadtxt(input_file, skiprows=1, usecols = (1,2,3,4,5,7,8,9,10))

rad_planet = data[:,0] * 0.10049 # planet radius in [solar radii] 
period     = data[:,1] * 24      # planet's orbital period in [hours] 
semi_major = data[:,2] * 215.032 # semi-major axis in [solar radii] 
epoch_peri = data[:,3]           # "epoch of periastron" in [JD] 
imp_param  = data[:,4]           # impact parameter 
right_asc  = data[:,5]           # [deg]
dec        = data[:,6]           # [deg]
V_mag      = data[:,7]           # magnitude of star in V-band (apparent)
rad_star   = data[:,8]           # star radius in [solar radii] 
print(planet_name[0], star_name[0], rad_planet[0], period[0], semi_major[0], epoch_peri[0], imp_param[0], right_asc[0], dec[0], V_mag[0], rad_star[0])

transit_depth = (rad_planet/rad_star)**2
print(transit_depth[0])

