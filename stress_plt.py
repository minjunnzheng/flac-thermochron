#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:35:10 2020

@author: summer-tan2
"""

import sys
sys.path.append("/home/summer-tan2/mingjun/geoflac/util")
import os
import flac
import math
import matplotlib.pyplot as plt
fl=flac. Flac()

##read data
x=300
density=fl.read_density(1)
density_1D=density[x,:]
xmesh,zmesh=fl.read_mesh( 1 )
viscosity=fl.read_visc( 1 )
vis_1D=viscosity[x,:]
fric1=30
coh1=4*10**7
fric=(fric1*math.pi)/180
Z_location_1D=[]

z_1D=zmesh[300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300300,:]

for i in range(len(zmesh[0,:])-1):
    
    Z_location_1D.append((zmesh[x,i]+zmesh[x,i+1]+zmesh[x+1,i]+zmesh[x+1,i+1])/4)
##yield_stress
pressure=[]
yield_stress=[]

for i in range(len(Z_location_1D)):
    temp=abs(density_1D[i]*9.8*Z_location_1D[i]*1000)
    pressure.append(temp)
    yi=coh1+math.tan(fric)*pressure[i]    
    yield_stress.append(yi*(10**(-6)))
    
##elast_stress
elast_stress=[]
for i in range(len(vis_1D)):
    elast_stress.append((10**vis_1D[i])*(10**(-15)*(10**(-6))))

##stress
stress=[]
for j in range(len(yield_stress)):
    if (yield_stress[j]<elast_stress[j]):
        stress.append(yield_stress[j])
        
    else:
        stress.append(elast_stress[j])

#Z_location_1D.reverse()
plt.plot(stress,Z_location_1D)
plt.xlabel("stress(MPa)")
plt.ylabel("depth(km)")
plt.show()
plt.plot(vis_1D,Z_location_1D)
plt.xlabel("viscosity(log)")
plt.ylabel("depth(km)")
plt.show()
#plt.invert_yaxis()      
#plt        

temperature=fl.read_temperature(1)
t_1D=temperature[150,:]
plt.plot(t_1D,z_1D)
plt.xlabel("temperature('C)")
plt.ylabel("depth(km)")


plt.figure(figsize = (16, 10), dpi = 80)
#plt.subplots_adjust(left = 0.1, bottom = 0.2, right = 0.9, top = 0.9, wspace = 0.4, hspace = 0.1)
plt.subplot(131)
plt.plot(stress,Z_location_1D)
plt.title("stress vs depth")
plt.xlabel("stress(MPa)")
plt.ylabel("depth(km)")
plt.subplot(132)
plt.title("viscosity vs depth")
plt.plot(vis_1D,Z_location_1D)
plt.xlabel("viscosity(log)")
plt.ylabel("depth(km)")
plt.subplot(133)
plt.title("temperature vs depth")
temperature=fl.read_temperature(1)
t_1D=temperature[150,:]
plt.plot(t_1D,z_1D)
plt.xlabel("temperature('C)")
plt.ylabel("depth(km)")

