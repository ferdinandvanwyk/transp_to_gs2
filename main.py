import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as interp
from netCDF4 import Dataset

in_file =  str(sys.argv[1])
ncfile = Dataset(in_file, 'r', format='NETCDF3')

#Read relevant data from TRANSP file
time = ncfile.variables['TIME']
print time.shape
