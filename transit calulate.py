import math as m

b = float(input('enter impact parameter b '))
P = float(input('enter period (days)'))
a = float(input('enter semi major of the planet (AU) '))
Rs = float(input('enter host star radius (R_sun) '))
Rp = float(input('enter planet radius in (R_jup) '))

# Convert jupiter radius into solar radius
rp = 0.10045*Rp
# Convert days into minute
p = 1440.0*P
# Convert AU in solar radius
A = 215.032*a


T_dur = p/m.pi*m.asin(m.sqrt(((Rs+rp)**2-(b*Rs)**2)/A))

T_dep = (rp/Rs)**2

print('Transit duration is (min):', T_dur)
print('Transit dpeth is:', T_dep)
