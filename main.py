import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as interp
from scipy.io import netcdf

in_file =  str(sys.argv[1])
ncfile = netcdf.netcdf_file(in_file, 'r')

#Read relevant data from TRANSP file
time = ncfile.variables['TIME']
