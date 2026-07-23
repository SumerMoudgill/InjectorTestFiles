#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 16:31:28 2022

@author: Selina Yang adapted from Andrew Sullivan's code. 
"""

import numpy as np
import csv
from scipy.integrate import quad
from scipy import optimize
import random as rand 

#GENERATE REDSHIFT INFORMATION
#code from andrew sullivan that takes account of the published ligo O3 mass and distance distribution

M=20000 #the number of pairs of masses you want
event_type='BBH' #has the option of BBH, BNS, NSBH 

#From now on it should do everything for BBH, BNS,NSBH case
if event_type=='BBH':
    Max_z=0.74 #for subthreshold maps 
elif event_type=='BNS':
    Max_z=0.087
elif event_type=='NSBH':
    Max_z=0.157 #for subthreshold maps


H0=67.6#hubble constant 
d1=lambda x: ((1+x)**3*0.31+0.69)**-0.5 #lambda function is a quick way to define a mathematical function 
def fr(z):#Function which returns luminosity distance in Mpc given redshift z
    return 3*10**5*(1+z)*quad(d1,0,z)[0]/H0 #returns the luminosity distance in Mpc

'''
astropy has a comoving distance calculator, the input is z and the output is the comoving distance in Mpc
Planck15 is a build-in cosmology model that assumes H0 = 67.7 and Om = 0.307 on a flat universe

from astropy.cosmology import Planck15
co_dist=Planck15.comoving_distance(1.25)
print(co_dist)
'''

Zs=np.empty(M)#create an empty redshift array with length = number of mass
Rs=np.empty(M)#create an empty lum distance array with length = number of mass

Zs[0]=Max_z*np.random.rand()#generate first redshift of the array
rold=fr(Zs[0])#generate first distance info based on that redshift
Rs[0]=rold
probold= rold**2*(1+Zs[0])**2.7/(1+Zs[0])**4
#Perform Metropolis Monte Carlo to get distances
for n in range(0, M-1):
    Ztest=Max_z*np.random.rand()
    rtest=fr(Ztest)
    probtest= rtest**2*(1+Ztest)**2.7/(1+Ztest)**4
    if probtest/probold > np.random.rand():
         Zs[n+1]=Ztest
         Rs[n+1]=rtest
         rold=rtest
         probold=probtest
    else:
         Zs[n+1]=Zs[n]
         Rs[n+1]=rold
         probold=probold

'''
Zs=np.random.uniform(0,Max_z,size=M)
Rs=[fr(z) for z in Zs]
Rs=np.array(Rs)
'''
# Zs is our redshift array
# Rs is our distance array
#print(Zs, Rs)


#GENERATE MASS INFORMATION BASED ON LIGO DISTIRBUTION
#Parameters for the power law plus peak
#andrew sullivan's code that takes O3B mass distribution 

#Parameters for the power law plus peak
max_BH=100#in solar masses
min_BH=5.1
alpha=3.4
beta=1.1
lam=0.039
masspeak=34
peaksig=3.6
delta_m=4.8

Bnorm=(alpha-1)/(min_BH**(1-alpha)-max_BH**(1-alpha))
Gnorm=1/peaksig*(2*np.pi)**(-1/2)
mass1=np.empty(M)
mass2=np.empty(M)
mass3=np.empty(M)
mass1[0]=min_BH+(max_BH-min_BH)*np.random.rand()
if mass1[0] < min_BH+delta_m:
    probold=((1-lam)*Bnorm*mass1[0]**(-alpha)+lam*Gnorm*np.exp(-1/2*((mass1[0]-masspeak)/peaksig)**2))*(1+np.exp(delta_m/(mass1[0]-min_BH)+delta_m/(mass1[0]-min_BH-delta_m)))**(-1)
else:
    probold=(1-lam)*Bnorm*mass1[0]**(-alpha)+lam*Gnorm*np.exp(-1/2*((mass1[0]-masspeak)/peaksig)**2)
integ=lambda x: x**beta*(1+np.exp(delta_m/(mass1[0]*x-min_BH)+delta_m/(mass1[0]*x-delta_m)))**(-1)
if mass1[0] > (min_BH+delta_m):
        B=quad(integ, min_BH/mass1[0], (min_BH+delta_m)/mass1[0])[0]
        C=(1-((min_BH+delta_m)/mass1[0])**(beta+1))/(beta+1)
        A=1/(B+C)
        
        pq2=np.random.rand()
        if pq2 > A*B:
            m2=mass1[0]*(((beta+1)*(pq2-A*B)/A+((min_BH+delta_m)/mass1[0])**(beta+1)))**(1/(beta+1))
        else:
            def f(x):
                bro=pq2-A*quad(integ, min_BH/mass1[0], x)[0]
                return bro
            m2=mass1[0]*optimize.bisect(f, min_BH/mass1[0]+0.0001, (min_BH+delta_m)/mass1[0])
        mass2[0]=m2
        
        pq3=np.random.rand()
        if pq3 > A*B:
            m3=mass1[0]*(((beta+1)*(pq3-A*B)/A+((min_BH+delta_m)/mass1[0])**(beta+1)))**(1/(beta+1))
        else:
            def f(x):
                bro=pq3-A*quad(integ, min_BH/mass1[0], x)[0]
                return bro
            m3=mass1[0]*optimize.bisect(f, min_BH/mass1[0]+0.0001, (min_BH+delta_m)/mass1[0])
        mass3[0]=m3
else:
        B=quad(integ, min_BH/mass1[0], 1)[0]
        
        A=1/(B)
        
        pq2=np.random.rand()
        
        def f(x):
            bro=pq2-A*quad(integ, min_BH/mass1[0], x)[0]
            return bro
        m2=mass1[0]*optimize.bisect(f, min_BH/mass1[0]+0.0001, 1)
        mass2[0]=m2
        
        pq3=np.random.rand()
        
        def f(x):
            bro=pq3-A*quad(integ, min_BH/mass1[0], x)[0]
            return bro
        m3=mass1[0]*optimize.bisect(f, min_BH/mass1[0]+0.0001, 1)
        mass3[0]=m3
        
for n in range(0, M-1):
    newmass=min_BH+(max_BH-min_BH)*np.random.rand()
    if newmass < min_BH+delta_m:
        probnew=((1-lam)*Bnorm*newmass**(-alpha)+lam*Gnorm*np.exp(-1/2*((newmass-masspeak)/peaksig)**2))*(1+np.exp(delta_m/(newmass-min_BH)+delta_m/(newmass-min_BH-delta_m)))**(-1)
    else:
        probnew=(1-lam)*Bnorm*newmass**(-alpha)+lam*Gnorm*np.exp(-1/2*((newmass-masspeak)/peaksig)**2)
    
    if probnew/probold >np.random.rand():
        probold=probnew
        mass1[n+1]=newmass
    else:
        mass1[n+1]=mass1[n]
        
    #Determine the second mass
    
    integ=lambda x: x**beta*(1+np.exp(delta_m/(mass1[n+1]*x-min_BH)+delta_m/(mass1[n+1]*x-delta_m)))**(-1)
    if mass1[n+1] > (min_BH+delta_m):
        B=quad(integ, min_BH/mass1[n+1], (min_BH+delta_m)/mass1[n+1])[0]
        C=(1-((min_BH+delta_m)/mass1[n+1])**(beta+1))/(beta+1)
        A=1/(B+C)
        
        pq=np.random.rand()
        if pq > A*B:
            m2=mass1[n+1]*(((beta+1)*(pq-A*B)/A+((min_BH+delta_m)/mass1[n+1])**(beta+1)))**(1/(beta+1))
        else:
            def f(x):
                bro=pq-A*quad(integ, min_BH/mass1[n+1], x)[0]
                return bro
            m2=mass1[n+1]*optimize.bisect(f, min_BH/mass1[n+1]+0.0001, (min_BH+delta_m)/mass1[n+1])
        mass2[n+1]=m2
    else:
        B=quad(integ, min_BH/mass1[n+1], 1)[0]
        
        A=1/(B)
        
        pq=np.random.rand()
        
        def f(x):
            bro=pq-A*quad(integ, min_BH/mass1[n+1], x)[0]
            return bro
        m2=mass1[n+1]*optimize.bisect(f, min_BH/mass1[n+1]+0.0001, 1)
        mass2[n+1]=m2
        
    #Determine the mass of the randomly generated bh in the second merger
    
    integ=lambda x: x**beta*(1+np.exp(delta_m/(mass1[n+1]*x-min_BH)+delta_m/(mass1[n+1]*x-delta_m)))**(-1)
    if mass1[n+1] > (min_BH+delta_m):
        B=quad(integ, min_BH/mass1[n+1], (min_BH+delta_m)/mass1[n+1])[0]
        C=(1-((min_BH+delta_m)/mass1[n+1])**(beta+1))/(beta+1)
        A=1/(B+C)
        
        pq=np.random.rand()
        if pq > A*B:
            m3=mass1[n+1]*(((beta+1)*(pq-A*B)/A+((min_BH+delta_m)/mass1[n+1])**(beta+1)))**(1/(beta+1))
        else:
            def f(x):
                bro=pq-A*quad(integ, min_BH/mass1[n+1], x)[0]
                return bro
            m3=mass1[n+1]*optimize.bisect(f, min_BH/mass1[n+1]+0.0001, (min_BH+delta_m)/mass1[n+1])
        mass3[n+1]=m3
    else:
        B=quad(integ, min_BH/mass1[n+1], 1)[0]
        
        A=1/(B)
        
        pq=np.random.rand()
        
        def f(x):
            bro=pq-A*quad(integ, min_BH/mass1[n+1], x)[0]
            return bro
        m3=mass1[n+1]*optimize.bisect(f, min_BH/mass1[n+1]+0.0001, 1)
        mass3[n+1]=m3
   
#newmass=np.linspace(min_BH/max_BH, (min_BH+delta_m)/max_BH, 1000)*max_BH
#p=(newmass/max_BH)**beta*(1+np.exp(delta_m/(newmass-min_BH)+delta_m/(newmass-min_BH-delta_m)))**(-1)
#z=np.load('ztriple3.npy')


if event_type=='BBH':
    pass
elif event_type=='BNS':
    mass1=np.full(M,1.4)#we use 1.4 Solar Mass for neutron stars in simulation 
    mass2=np.full(M,1.4)
elif event_type=='NSBH':
    mass2=np.full(M,1.4)


mass1=mass1*(1+Zs) #this is to account for the distance
mass2=mass2*(1+Zs)

#index= list(range(0, len(mass1)))

#GENERATING LOCATION INFO
decl=np.empty(M)
decl1=[]


ra=np.empty(M)
ra1=[]

ran=0
name=np.empty(M)
lum=[]
fudge=[]
for m in range(0, M):
    name[m]=m
    lum.append(' 10')
    fudge.append(' 10')
for m in range(0, M):
    sindeg=2*rand.random()-1 #since declination is uniformly distributed in sin
    deg=np.arcsin(sindeg)*180/np.pi
   
    if deg >= 0:
        decl1.append('+'+ str(int(deg))+':'+str(abs(int((deg-int(deg))*60))))
    else:
        if int(deg)>-1:
            decl1.append('-'+str(int(deg))+':'+str(abs(int((deg-int(deg))*60))))
        else:
            decl1.append(str(int(deg))+':'+str(abs(int((deg-int(deg))*60))))
for m in range(0, M):
    ran=24*rand.random()
    ra1.append('+'+str(int(ran))+':'+(str(abs(int((ran-int(ran))*60)))))

Rs=Rs*1000

with open('o4sub_'+event_type+'_mass.dat', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(mass1, mass2))
with open('o4sub_'+event_type+'_source.dat', 'w') as g:
    writer = csv.writer(g, delimiter='\t')
    writer.writerows(zip(name, ra1, decl1, Rs, lum, fudge))#use zip if things are in array, if not array, no zip
