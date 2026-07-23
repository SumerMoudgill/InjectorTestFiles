"""
Created on Mon Dec 5 16:31:28 2022

@author: Selina Yang adapted from Andrew Sullivan's code. 
email: fy2265@columbia.edu

Use this code for the actual bayestar injection (either submited to habanero or on local computer)
"""

import numpy as np
import csv 
import random 
import os # a package under the python standard utility modules, should work for most of the common operation systems
import subprocess
import os.path
import pandas as pd
from xml.dom import minidom

from astropy.coordinates import SkyCoord
from ligo.skymap.moc import rasterize
from ligo.skymap.io.fits import read_sky_map
from lal import antenna

from astropy.io import fits
from astropy.table import Column, Table, join
from astropy import units as u

import shutil

#load in the mass and source files you generated in the "MonteCarlo_Gen.py" code. 
Mass=np.loadtxt('o4sub_BBH_mass.dat')
Source=np.genfromtxt('o4sub_BBH_source.dat', dtype=str)
#Mass=np.loadtxt('o4sub_BBH_mass.dat')
#Source=np.genfromtxt('o4sub_BBH_source.dat', dtype=str)

M=100#0 #how many good skymaps you want out of this,
#should be way below the number of mass pairs you generated from the Monte Carlo
#due to low sensitivity on the detectors, furthur distance, smaller mass, or other reasons
#especially if you know that your bayestar might not be able to detect some of the injections
oldname = '0.fits'


#reference_psd='o3sub_psd.xml'
reference_psd = 'H1L1_160_PSD.xml'
eventtype='BBH'
#eventtype='BBH'
detector_combo='L1 H1' #please put space in the middle of the two detector names

injectioninfofilename='lhbbh_injectioninfo.dat'
goodinjection_foldername='o4sub_lh_bbh'
#goodinjection_foldername='o4sub_lh_bbh'
badinjection_foldername='o4trash'

good_runs=0 #a tally of how many good runs it has been so far #good run means the injection produced a skymap
n = 0
bad_runs=0
while good_runs < M:
    gpsstart=1238112018
    gpsend = 1382016188 #S231022bd
    interval = gpsend-gpsstart
    #interval=np.random.uniform(86400,31536000) #from a day to 1 years # CHANGED INTERVAL AND GPSEND TO INTEGERS
    #interval = int(np.random.uniform(86400,313536000))
    file_number=1
    #gpsend=gpsstart+(interval*file_number)
    #gpsend = int(gpsstart+(interval*file_number))
    m1=Mass[n, 0]
    m2=Mass[n, 1]
    #this is the first merger
    mass_file1=[m1,m2]

    dfmass = pd.read_csv('o4sub_'+eventtype+'_mass.dat',delimiter='\t',header=None)
    dfsource = pd.read_csv('o4sub_'+eventtype+'_source.dat',delimiter='\t',header=None)

    thismassinfo=list(dfmass.loc[n].to_numpy())
    thissourceinfo=list(dfsource.loc[n].to_numpy())
    timeinfo=[gpsstart,interval,gpsend]


#write the mass and location file for that specific injection 
    with open('thismass.dat', 'w') as f: #M1 stands for Merger 1
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(mass_file1)
        
    with open('thissource.dat', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(Source[n])

#deciding which waveform to use based on the type of event
    if m1+m2<5: #the BNS case
        waveform='TaylorF2twoPointFivePN'
        #Zsuzsa said that TaylorF2 is acceptable for BNS for 190425 paper 
        #Rainer used the TaylorF2threePointFivePN
    elif m2<=3: #the NSBH case
        waveform='IMRPhenomNSBH' #one of the masses has to be 3 solar mass or smaller. 
        #found on https://arxiv.org/pdf/2010.14527.pdf
    else: #The BBH case
        waveform='SEOBNRv4twoPointFivePN'

    waveformlist=[waveform]
    otherlalappsinfo=['uniform','source','source','disable_milkyway','source',20] #t-distr, m-distr, d-distr,milkyway,l-distr,f-lower
    otherrealizecoincsinfo=[8,1,20] #snr-threshold,mintrigger,f-low 
    

    #randomize inclination: 
    inclination=np.rad2deg(np.arccos(np.random.uniform(-1.0,1.0)))


    injectioninforow=thismassinfo+thissourceinfo+timeinfo+waveformlist+otherlalappsinfo+otherrealizecoincsinfo+[inclination]

#doing the actual injection 
    injection=(
        f'lalapps_inspinj -o inj'+str(n)+'.xml --source-file thissource.dat --mass-file thismass.dat --t-distr uniform '
        f'--time-step {interval} --gps-start-time {gpsstart} --gps-end-time {gpsend} '
        f'--m-distr source --d-distr source --disable-milkyway --l-distr source --i-distr fixed --fixed-inc {inclination} '
        f'--disable-spin --f-lower 20 --waveform {waveform}'
    )
    os.system(injection)

    #since we already have an xml made: 

    realizecoincs='bayestar-realize-coincs -o coinc'+str(n)+'.xml inj'+str(n)+'.xml --reference-psd '+str(reference_psd)+' --detector '+str(detector_combo)+' --snr-threshold 4  --net-snr-threshold 8 --min-triggers 2 --keep-subthreshold --waveform '+str(waveform)+' --f-low 20'    #coinc.xml has the data of the output, recovered snr, mass, etc etc
    os.system(realizecoincs)

#since we don't need maps, determine whether or not this is a successful run using realize coincs
    file = minidom.parse('coinc'+str(n)+'.xml')
    element = file.getElementsByTagName('LIGO_LW')
    xmlnamelist=[]
    for i in element: 
        thisname=str(i.getAttribute("Name"))
        xmlnamelist.append(thisname)
    if 'COMPLEX8TimeSeries' in xmlnamelist: 
        successtoken=1
    else: 
        successtoken=0

#if and when it makes a fits file, move it to a subfolder titled "MadeSkymaps" and also move the coinc file that contains the output data with it. 
    if successtoken==1:
        with open(injectioninfofilename, 'a', newline='') as myfile:
            wr = csv.writer(myfile, delimiter='\t')
            wr.writerow(injectioninforow+['success'])
        localizecoincs = (
            f'bayestar-localize-coincs coinc'+str(n)+'.xml --f-low 20 --waveform '+str(waveform)+' --f-low 20')
        os.system(localizecoincs)
        newname = 'fits'+str(good_runs)+'.fits'
        #dfmass['success'].loc[n]='success'
        #dfsource['success'].loc[n]='success'
        if os.path.exists('inj'+str(good_runs)+'.xml'):
            os.system('mv inj'+str(n)+'.xml ./'+goodinjection_foldername+'/inj'+str(good_runs)+'.xml')
        if os.path.exists('coinc'+str(n)+'.xml'):
            os.system('mv coinc'+str(n)+'.xml ./'+goodinjection_foldername+'/coinc'+str(good_runs)+'.xml')
            good_runs+=1
        if os.path.exists(oldname):
            os.rename(oldname, newname)
            os.system(newname+' ./'+goodinjection_foldername+'/'+newname)
        successtoken=2
    elif successtoken==0: 
        with open(injectioninfofilename, 'a', newline='') as myfile:
            wr = csv.writer(myfile, delimiter='\t')
            wr.writerow(injectioninforow+['failed'])

        if os.path.exists('inj'+str(n)+'.xml'):
            os.system('mv inj'+str(n)+'.xml ./'+badinjection_foldername+'/inj'+str(n)+'.xml')
        if os.path.exists('coinc'+str(n)+'.xml'):
            os.system('mv coinc'+str(n)+'.xml ./'+badinjection_foldername+'/coinc'+str(n)+'.xml')
            bad_runs+=1
        successtoken=2
    n += 1


