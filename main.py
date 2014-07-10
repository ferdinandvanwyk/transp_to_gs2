import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as interp
from netCDF4 import Dataset

in_file =  str(sys.argv[1])
output_radius =  sys.argv[2]
output_time =  sys.argv[3]
ncfile = Dataset(in_file, 'r', format='NETCDF3')

print in_file, output_radius, output_time

#Read relevant data from TRANSP file
time = ncfile.variables['TIME']
x = ncfile.variables['X']
xb = ncfile.variables['XB']
ni = ncfile.variables['NI']
ne = ncfile.variables['NE']
ti = ncfile.variables['TI']
te = ncfile.variables['TE']
fbtx = ncfile.variables['FBTX']
btx = ncfile.variables['BTX']
rmaj = ncfile.variables['RMAJM']
bpol = ncfile.variables['BPOL']
omega = ncfile.variables['OMEGA']
q = ncfile.variables['Q']
shat = ncfile.variables['SHAT']
flux_centres = ncfile.variables['RMJMP']
elongation = ncfile.variables['ELONG']
triang = ncfile.variables['TRIANG']
zeffp = ncfile.variables['ZEFFP']

#Normalization quantities
amin = 0.58044
vth = 0  

