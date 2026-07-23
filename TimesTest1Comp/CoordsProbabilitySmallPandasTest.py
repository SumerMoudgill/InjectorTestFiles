#!/usr/bin/python3 
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from hpmoc import PartialUniqSkymap
from hpmoc.plot import get_wcs, plot, gridplot
import sys
from CoordsProbabilityCFolder.DetectorMinima import *
from CoordsProbabilityCFolder.RankPoints import *
#from CoordsProbabilityCFolder.O4aSuperevents import *
#from CoordsProbabilityCFolder.O4cSuperevents import *
#from CoordsProbabilityCFolder.O4bSupereventsBilby import *
from CoordsProbabilityCFolder.BBH_test_260720_time1 import *
from CoordsProbabilityCFolder.CoordsProbabilityFunctionsTest import *
from CoordsProbabilityCFolder.SupereventReadTest import *
from CoordsProbabilityCFolder.DetectorMinimaSkymap import *
from CoordsProbabilityCFolder.ConversionF import *
import pandas as pd

#nesting:
#superevent_to_probabilities
#   which_detectors
#   get_probabilities_at_minima
#       get_probability_at_coords
#           match_coords
#           single_pixel_to_probability
#superevent_array_map



#event_dataset = ["event_name", "timeJulian", "detectors", "zero_locations", "detectors_triggered", "sigfigs", "values", "skymap_coord_string"]
#for i in range(len(superevents)):
#    current_superevent=superevents[i]
#    event_name=current_superevent[0]
#    timeJulian=current_superevent[2]
#    detectors=which_detectors(event_array[3])
#    current_superevent_array=superevents[i]
#    print(superevent_to_probabilities(current_superevent_array, file_location="./testing/"))
    #zero_locations from get_probabilities_at_minima, through superevent_to_probabilities
    #detectors_triggered from get_probabilities_at_minima, through superevent_to_probabilities 
    #sigfigs from math_coords, through get_probability_at_coords, get_probabilities_at_minima, and superevent_to_probabilities
    #values from match_coords, through get_probability_at_coords, get_probabilities_at_minima, and superevent_to_probabilities
    #skymap_coord_string from get_probability_at_coords, through get_probabilities_at_minima, and superevent_to_probabilities

#print("hello world")
#fits_a='./detected_fitsmaps/maps1_1.fits'
#fits_b='./detected_fitsmaps/maps2_1.fits'
#fits_c='./detected_fitsmaps/maps4_1.fits'
#a = PartialUniqSkymap.read(fits_a, strategy='ligo')
#b = PartialUniqSkymap.read(fits_b, strategy='ligo')
#c = PartialUniqSkymap.read(fits_c, strategy='ligo')
#print(type(a))
#print(a)
#print(a.coords())
#print(type(a.plot()))
#print(a.plot())

#GW170817_array=["GW170817", 1, GW170817_JDT, "H1L1V1"]
#print(superevent_to_probabilities(GW170817_array, file_location="./testing/"))
#superevent_array_map(GW170817_array, "GW170817", file_location="./testing/")
namescolumn=[]
timescolumn = []
detectorscolumn = []
sf_column = []
p_column = []
f_column = []
cf_column = []
maxpoint_column = []
max_sep_column = []
array_to_use = BBH_test_260720_1_array
print(array_to_use)
folder_d = "./BBH_test_260615_fits/"
for i in range(len(array_to_use)):
    current_superevent_array=array_to_use[i]
    current_superevent_read=superevent_read(current_superevent_array, file_location=folder_d)
    #print("Event name: ", current_superevent_array[0])
    namescolumn.append(current_superevent_array[0])
    #print("Event time: ", current_superevent_array[2])
    timescolumn.append(current_superevent_array[2])
    #print("Detectors list: ", which_detectors(current_superevent_array[3]))
    detectorscolumn.append(which_detectors(current_superevent_array[3]))
    n_detectors=len(which_detectors(current_superevent_array[3]))
    res_p, res_sf, res_c, res_f=superevent_to_probabilities(current_superevent_array, file_location=folder_d)
    #print("Probabilities list: ", res_p)
    p_column.append(res_p)
    #print("Significant figures list: ", res_sf)
    sf_column.append(res_sf)
    #print("Coordinates list: ", res_c)
    #c_column.append(res_c)
    #print("Flags list: ", res_f)
    f_column.append(res_f)
    cf_list=get_all_null_coords(current_superevent_array[2], which_detectors(current_superevent_array[3]))
    res_maxpoint = max_point(current_superevent_read)
    #print("Max point: ", res_maxpoint)
    maxpoint_column.append(res_maxpoint)
    #print("all null points: ", cf_list)
    cf_column.append(cf_list)
    res_sep_min_s=np.zeros((3, 4))
    #print(res_f)
    #print(type(res_f))
    for i in range(len(res_f)):
        for j in range(len(res_f[i])):
            #print("j: ", j)
            if len(cf_list[i][j]) != 0:
                res_sep_min_s[i][j]=angular_separation_deg(res_maxpoint, cf_list[i][j])
    max_sep_column.append(res_sep_min_s.tolist())
    print("RES.OUT")
    #print("90% probability area: ", get_area("./O4A_fits/", current_superevent_array, area_percent=90))
    superevent_array_map(current_superevent_array, current_superevent_array[0], file_location=folder_d) 
results_d = {'Event name':namescolumn, 'Event time':timescolumn, 'Detectors':detectorscolumn, 'Probabilities':p_column, 'Significant figures':sf_column, 'Flags':f_column, 'All null points':cf_column, 'Max point':maxpoint_column, 'Separations from maximum point':max_sep_column}
results_df=pd.DataFrame(data=results_d)
print(results_df)
results_df.to_excel("O4_BBH_test_time1A.xlsx") 
