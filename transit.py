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
print(np.shape(data))

rad_planet = data[:,0] * 0.10049 # planet radius in [solar radii] 
period     = data[:,1] * 24      # planet's orbital period in [hours] 
semi_major = data[:,2] * 215.032 # semi-major axis in [solar radii] 
epoch_peri = data[:,3]           # "epoch of periastron" in [JD] 
imp_param  = data[:,4]           # impact parameter 
right_asc  = data[:,5]           # [deg]
dec        = data[:,6]           # [deg]
V_mag      = data[:,7]           # magnitude of star in V-band (apparent)
rad_star   = data[:,8]           # star radius in [solar radii] 
#print(planet_name[0], star_name[0], rad_planet[0], period[0], semi_major[0], epoch_peri[0], imp_param[0], right_asc[0], dec[0], V_mag[0], rad_star[0])

transit_depth = (rad_planet/rad_star)**2
print(transit_depth[0])

# Range for declination 

dec_min = 21 # [deg]
dec_max = 61 # [deg]

ra_max  = 19 # [hour] i.e. 7PM 
ra_min  = 1 # [hour] i.e. 1AM

# Initial limitations for selecting candidates by magnitude, declination, right ascension 

good = np.where( (V_mag < 12.5) & (dec < dec_max) & (dec > dec_min) & (right_asc > 15) &  (right_asc >ra_max*15) | (right_asc < ra_min * 15) )
data_good = data[good] 

# good data now: 
star_name = star_name[good]
planet_name = planet_name[good] 
rad_planet = rad_planet[good] # [solar radii]
period = period[good]         # [hours] 
semi_major = semi_major[good] # [solar radii]
imp_param = imp_param[good]
right_asc = right_asc[good]   # [deg]
dec = dec[good]               # [deg]
V_mag = V_mag[good]
rad_star = rad_star[good]     # [solar radii] 

print(np.shape(data_good)) # print number of candidates left for selection 
print(np.shape(star_name))
print(np.shape(planet_name))

# New file with these candidates 
with open('./exoplanet_candidates.txt', 'w+') as file:
    for planet_name, rad_planet, period, semi_major, imp_param, star_name,  right_asc, dec, V_mag, rad_star in zip(planet_name, rad_planet, period, semi_major, imp_param, star_name, right_asc, dec, V_mag, rad_star): 
        file.write('{:20}\t{:20}\t{:20}\t{:20}\t{:20}\t{:20}\t{:20}\t{:20}\t{:20}\t{:20}\n'.format(planet_name, rad_planet, period, semi_major, imp_param, star_name, right_asc, dec, V_mag, rad_star))

file.close() 
   

