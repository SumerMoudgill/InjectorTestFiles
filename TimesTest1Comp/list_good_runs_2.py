import csv
import numpy as np
from astropy.time import Time
import math

def gps2jd(timeGPS):
    t=Time(timeGPS, format="gps")
    t.format="jd"
    #print(t.value)
    return t.value
#list_1 = [3,6,9,14,16,37,38,39,40,60,61,66,78,79,81,82,84,98,100,104,105,107,108,113,125,131,139,146,162,163,177,188,189,192,193,203,204,205,207,215,217,222,223,225,227,238,241,250,259,262,266,268,272,275,276,279,280,281,288,289,290,291,292,293,294,295,297,309,320,321,328,330,331,335,336,337,342,343,344,349,350,358,369,370,383,387,388,394,395,396,397,399,403,410,411,416,425,427,428,429]

#with open ('lhbns_injectioninfo.csv', mode='r', newline = '', encoding = 'utf-8') as file:
with open ('lhbbh_injectioninfo.csv', mode='r', newline = '', encoding = 'utf-8') as file:
    reader = csv.reader(file, delimiter='\t')
    newlist = list(reader)
print(newlist)
print("line")
print(newlist[0])
print("element 22")
print(newlist[0][22])

#print(len(list_1))
list_t_gps = []
list_t_jdt = []
i=0
#for n in range(len(newlist)):
#    if newlist[n][22]=="success":
#        listelement=[]
#        t_jdt = gps2jd(newlist[n][10])
#        listelement_gps = [n, i, newlist[n][1],"H1L1"]
#        list_t_gps.append(listelement_gps)
#        listelement_jdt = [n, i, t_jdt,"H1L1"]
#        list_t_jdt.append(listelement_jdt)
#        list_t_gps.append(float(newlist[n][10]))
#        i=i+1
for n in range(len(newlist)):
    if newlist[n][22]=="success":
        listelement=[]
        t_jdt = gps2jd(newlist[n][10])
        listelement_gps = [n, i, newlist[n][10],"H1L1"]
        list_t_gps.append(listelement_gps)
        listelement_jdt = [n, i, t_jdt,"H1L1"]
        list_t_jdt.append(listelement_jdt)
        list_t_gps.append(float(newlist[n][10]))
        i=i+1


#for i in range(len(list_1)):
#    listelement = []
#    n = list_1[i]
#    t_jdt = gps2jd(newlist[n][10])
#    listelement_gps = [n, i, newlist[n][10],"H1L1"]
#    list_t_gps.append(listelement_gps)
#    listelement_jdt = [n, i, t_jdt,"H1L1"]
#    list_t_jdt.append(listelement_jdt)
#    list_t_gps.append(float(newlist[n][10]))
print("GPS")
print(list_t_gps)
print("JDT")
print(list_t_jdt)
