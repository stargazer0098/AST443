import numpy as np  
import matplotlib.pyplot as plt 
from datetime import datetime, timedelta, timezone 

input_file = 'exoplanets_all_list.txt'

''' 
Input file data structure: 

col 0: name of planet 
col 1: planet radius (jupiter radius)
col 2: orbital period of planet (days)  
col 3: semimajor axis (AU) 
col 4: impact paramater 
col 5: name of host star 
col 6: right ascension (deg) 
col 7: declination (deg)
col 8: V-band magnitude (apparent) 
col 9: star radius (solar radii) 
col 10: primary transit (JD) 
''' 
planet_name = np.loadtxt(input_file, dtype = str , skiprows = 1, usecols = 0)
star_name = np.loadtxt(input_file, dtype = str, skiprows = 1, usecols = 5)

print(planet_name)
print(star_name)

data = np.loadtxt(input_file, skiprows=1, usecols = (1,2,3,4,6,7,8,9,10))
print(np.shape(data))

rad_planet = data[:,0] * 0.10049 # planet radius in [solar radii] 
period     = data[:,1]           # planet's orbital period in [days]
semi_major = data[:,2] * 215.032 # semi-major axis in [solar radii] 
imp_param  = data[:,3]           # impact parameter 
right_asc  = data[:,4]           # [deg]
dec        = data[:,5]           # [deg]
V_mag      = data[:,6]           # magnitude of star in V-band (apparent)
rad_star   = data[:,7]           # star radius in [solar radii] 
prim_trans = data[:,8]           # [JD] 

#print(planet_name[0], star_name[0], rad_planet[0], period[0], semi_major[0], epoch_peri[0], imp_param[0], right_asc[0], dec[0], V_mag[0], rad_star[0])

transit_depth = (rad_planet/rad_star)**2
print(transit_depth[0])

# Range for declination 

dec_min = 21 # [deg]
dec_max = 61 # [deg]

ra_max  = 19 # [hour] i.e. 7PM 
ra_min  = 1 # [hour] i.e. 1AM

# Initial limitations for selecting candidates by magnitude, declination, right ascension 

good = np.where( (V_mag < 12) & (dec < dec_max) & (dec > dec_min) & (right_asc > 15) & ( (right_asc >ra_max*15) | (right_asc < ra_min * 15) ) )
data_good = data[good] 

# good data now: 
star_name = star_name[good]
planet_name = planet_name[good] 
rad_planet = rad_planet[good] # [solar radii]
period = period[good]         # [days]
semi_major = semi_major[good] # [solar radii]
imp_param = imp_param[good]
right_asc = right_asc[good]   # [deg]
dec = dec[good]               # [deg]
V_mag = V_mag[good]
rad_star = rad_star[good]     # [solar radii] 
prim_trans = prim_trans[good]
print(np.shape(data_good)) # print number of candidates left for selection 
print(np.shape(star_name))
print(np.shape(planet_name))

transit_depth = (rad_planet/rad_star)**2
print(prim_trans)
print(V_mag)

# Primary transit in JD represents when the planet is right in front of star (middle of its transit)
# so, if we keep adding periods (days) to the primary transit (JD) we basically would get times for when the planet is at the midpoint of its transit
 
planet = np.where(planet_name == 'HD189733b')
print(prim_trans[planet])

dt_Offset = 2400000.5 # offset between JD and MJD (modern julian date), starts at 11-17-1858

modified_julian_date = prim_trans - dt_Offset
dt_list = [ ]
print(modified_julian_date)

good = np.where(~np.isnan(modified_julian_date))  # only take indices where MDJ is not a nan 
print(good)
print(modified_julian_date[good])

modified_julian_date = modified_julian_date[good] 

# New York is 4 hours behind UTC time 

# today's julian day, i.e. august 31, 2021 is 2459458 JD (Wikipedia for Julian Day) 
today_JD = 2459458 # JD for august 31, 2021 
# 8/31/21 to october 8th is 38 so, iterate T0 + period until 2459458 

oct8_JD = today_JD + 38 

for i in range(len(modified_julian_date)):
   
    dt = datetime(1858, 11, 17, tzinfo = timezone.utc) + timedelta(modified_julian_date[i])
    #print(dt.date()) # get only the dates 

date_obs = prim_trans + period # in JD  
day_to_observe = [ ] 

for i in range(len(date_obs)):
    
    if (date_obs[i] < oct8_JD): 
    
        date_obs[i] = date_obs[i] + period[i] 
        
            date_obs[i] = date_obs[i] + period[i] 
        
    day_to_observe.append(date_obs[i])

print(date_obs)

print(day_to_observe)

#print(dt_list)



# PLOT 
# ---- 

#plt.plot(V_mag, transit_depth, '+')
#plt.show()


# New file with these candidates 

with open('./exoplanet_candidates.txt', 'w+') as file:
    for planet_name, rad_planet, period, semi_major, imp_param, star_name,  right_asc, dec, V_mag, rad_star, prim_trans in zip(planet_name, rad_planet, period, semi_major, imp_param, star_name, right_asc, dec, V_mag, rad_star, prim_trans): 
        file.write('{:15}\t{:15}\t{:15}\t{:15}\t{:15}\t{:15}\t{:15}\t{:15}\t{:15}\t{:15}\n'.format(planet_name, rad_planet, period, semi_major, imp_param, star_name, right_asc, dec, V_mag, rad_star, prim_trans))

file.close() 


 
