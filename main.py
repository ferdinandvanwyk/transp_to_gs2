import os, sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as interp
from netCDF4 import Dataset
import operator

#function that computes the debye length given T(keV), n(m^-3)
#Follows NRL plasma formulary
def debye(temp, dens):
  return 2.35e5*np.sqrt(temp)/np.sqrt(dens)

in_file =  str(sys.argv[1])
output_radius =  float(sys.argv[2])
output_time =  float(sys.argv[3])
ncfile = Dataset(in_file, 'r', format='NETCDF3')

#Read relevant data from TRANSP file
time = ncfile.variables['TIME'][:]

#Convert input time to an index
t_idx, min_value = min(enumerate(abs(time - output_time)), key=operator.itemgetter(1))

x = ncfile.variables['X'][t_idx,:]
xb = ncfile.variables['XB'][t_idx,:]
ni_tot = ncfile.variables['NI'][t_idx,:]
nd = ncfile.variables['ND'][t_idx,:] #Deuterium (reference species) density
nh = ncfile.variables['NH'][t_idx,:] #Hydrogen density
nimp = ncfile.variables['NIMP'][t_idx,:] #Impurity ion density
nb = ncfile.variables['BDENS'][t_idx,:] #Beam ion density
ne = ncfile.variables['NE'][t_idx,:]
ti = ncfile.variables['TI'][t_idx,:] #Combined D, H
timp = ncfile.variables['TX'][t_idx,:] #Impurity temp
tb = ncfile.variables['EBEAM_D'][t_idx,:] #Energy/temp of the beam ions
te = ncfile.variables['TE'][t_idx,:]
fbtx = ncfile.variables['FBTX'][t_idx,:]
bpbt = ncfile.variables['FBPBT'][t_idx,:]
btx = ncfile.variables['BTX'][t_idx,:]
rmaj = ncfile.variables['RMAJM'][t_idx,:]
bpol = ncfile.variables['BPOL'][t_idx,:]
omega = ncfile.variables['OMEGA'][t_idx,:]
q = ncfile.variables['Q'][t_idx,:]
shat = ncfile.variables['SHAT'][t_idx,:]
flux_centres = ncfile.variables['RMJMP'][t_idx,:]
elongation = ncfile.variables['ELONG'][t_idx,:]
triang = ncfile.variables['TRIANG'][t_idx,:]
zeffp = ncfile.variables['ZEFFP'][t_idx,:]
psi_t = ncfile.variables['TRFMP'][t_idx,:]

#Normalization quantities
boltz_jk = 1.3806488e-23
boltz_evk = 8.6173324e-5
proton_mass = 1.672621777e-27

#Need to calculate factor which relates gradients in TRANSP psi_n (=sqrt(psi_t/psi_tLCFS)) and Miller a_n (= diameter/diameter LCFS)
flux_rmaj = np.interp(output_radius, np.linspace(-1,1,rmaj.shape[0]), rmaj)
flux_idx, min_value = min(enumerate(abs(rmaj - flux_rmaj)), key=operator.itemgetter(1))
mag_axis_idx, min_value = min(enumerate(abs(bpbt)), key=operator.itemgetter(1))
a_left = (rmaj[flux_idx-1] - rmaj[mag_axis_idx - (flux_idx-1-mag_axis_idx)])/(rmaj[-1] - rmaj[0]) #diameter/diameter of LCFS
rho_miller = (rmaj[flux_idx] - rmaj[mag_axis_idx - (flux_idx-mag_axis_idx)])/(rmaj[-1] - rmaj[0]) #diameter/diameter of LCFS
a_right = (rmaj[flux_idx+1] - rmaj[mag_axis_idx - (flux_idx+1-mag_axis_idx)])/(rmaj[-1] - rmaj[0]) #diameter/diameter of LCFS
psi_left = np.sqrt(psi_t[flux_idx-1]/psi_t[-1]) #sqrt(psi_t/psi_LCFS)
rho_transp = np.sqrt(psi_t[flux_idx]/psi_t[-1]) #sqrt(psi_t/psi_LCFS)
psi_right = np.sqrt(psi_t[flux_idx+1]/psi_t[-1]) #sqrt(psi_t/psi_LCFS)
dpsi_da = (psi_right-psi_left)/(a_right-a_left) # coefficient which relates psi_n and a_n grids

# Gradient calculation requires numerical differentiation
# Use a second order central method: f'(x) = [f(x+h) - f(x-h)]/2h
# Find points closest points either side of output_radius and use those
rad_idx, min_value = min(enumerate(abs(x - output_radius)), key=operator.itemgetter(1))
radb_idx, min_value = min(enumerate(abs(xb - output_radius)), key=operator.itemgetter(1))

###################################
#Calculate equilibrium parameters #
###################################
equil = {}
equil['amin'] = (rmaj[-1] - rmaj[0])/2/100
amin = equil['amin']
vth = np.sqrt((2*np.interp(output_radius, x, ti)*1000*boltz_jk/boltz_evk)/(2*proton_mass))
equil['omega'] = np.interp(output_radius, x, omega) #rad/s
equil['btref'] = fbtx[mag_axis_idx]*btx[mag_axis_idx]
beta =  403.0*nd*1e6/1e19*ti/1000/(1e5*equil['btref']**2) #n_ref(1e19 m^-3), T_ref(keV)
beta_full =  403.0*(nd*ti + nh*ti + nimp*timp + nb*tb + ne*te)*1e6/1e19/1000/(1e5*equil['btref']**2) #n_ref(1e19 m^-3), T_ref(keV)
equil['beta'] = np.interp(output_radius, x, beta)
equil['zeff'] = np.interp(output_radius, x, zeffp)
equil['beta_prime_input'] = (beta_full[rad_idx+1]-beta_full[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])*dpsi_da #See wiki definition: not taking into account B_T variation
equil['g_exb'] = (omega[rad_idx+1]-omega[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])*(rho_miller/q[radb_idx])*(amin/vth)*dpsi_da #q defined on xb grid

f = open('gs2.in', 'w')
f.write('Equilibrium Parameters: \n')
f.write('----------------------- \n')
for name, value in sorted(equil.items()):
  f.write(name + ' = ' + str(value) + '\n')
f.write('\n')

###############################
#Calculate species parameters #
###############################

############
# ELECTRON #
############
electron = {}
electron['dens'] = np.interp(output_radius, x, ne)*1e6/1e19 #1e19m^-3
electron['mass'] = 1.0/(2.0*1836.0) #Assume D-electron plasma
electron['temp'] = np.interp(output_radius, x, te)/1000 #keV
electron['fprim'] = -(ne[rad_idx+1]-ne[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/ne[rad_idx]*dpsi_da
electron['tprim'] = -(te[rad_idx+1]-te[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/te[rad_idx]*dpsi_da

#See README for details of collision frequencies
log_e = 14.9 - 0.5*np.log(electron['dens']/10) + np.log(electron['temp']) # dens_1 already in 1e19m^-3 and temp_2 already in keV
tau_e = 1.09e16 * (electron['temp'])**1.5 / (electron['dens']*1e19*log_e) 
nu_e = 1./tau_e
electron['vnewk'] = nu_e*amin/vth

#ION SPECIES 1 - Deuterium (reference)
ion_1 = {}
ion_1['dens'] = np.interp(output_radius, x, nd)*1e6/1e19 #1e19m^-3
ion_1['mass'] = 1.0
ion_1['temp'] = np.interp(output_radius, x, ti)/1000 #keV
ion_1['fprim'] = -(nd[rad_idx+1]-nd[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/nd[rad_idx]*dpsi_da
ion_1['tprim'] = -(ti[rad_idx+1]-ti[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/ti[rad_idx]*dpsi_da

log_i1 = 17.3 - 0.5*np.log(ion_1['dens']/10) + 1.5*np.log(ion_1['temp']) # dens_1 already in 1e19m^-3 and temp_1 already in keV
tau_i1 = 6.6e17 * np.sqrt(2) * (ion_1['temp'])**1.5 / (ion_1['dens']*1e19*log_i1) # for singly charged deuterium ions => mi/mp=2 and Z=1 
nu_i1 = 1./tau_i1
ion_1['vnewk'] = nu_i1*amin/vth

#ION SPECIES 2 - Hydrogen
ion_2 = {}
ion_2['dens'] = np.interp(output_radius, x, nh)*1e6/1e19 #1e19m^-3
ion_2['mass'] = 0.5
ion_2['temp'] = np.interp(output_radius, x, ti)/1000 #keV
ion_2['fprim'] = -(nh[rad_idx+1]-nh[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/nh[rad_idx]*dpsi_da
ion_2['tprim'] = -(ti[rad_idx+1]-ti[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/ti[rad_idx]*dpsi_da

log_i2 = 17.3 - 0.5*np.log(ion_2['dens']/10) + 1.5*np.log(ion_2['temp']) # dens_1 already in 1e19m^-3 and temp_1 already in keV
tau_i2 = 6.6e17 * np.sqrt(2) * (ion_2['temp'])**1.5 / (ion_2['dens']*1e19*log_i2) # for singly charged deuterium ions => mi/mp=2 and Z=1 
nu_i2 = 1./tau_i2
ion_2['vnewk'] = nu_i2*amin/vth

#ION SPECIES 3 - Impurity
ion_3 = {}
ion_3['dens'] = np.interp(output_radius, x, nimp)*1e6/1e19 #1e19m^-3
ion_3['mass'] = 0.5
ion_3['temp'] = np.interp(output_radius, x, timp)/1000 #keV
ion_3['fprim'] = -(nimp[rad_idx+1]-nimp[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/nimp[rad_idx]*dpsi_da
ion_3['tprim'] = -(timp[rad_idx+1]-timp[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/timp[rad_idx]*dpsi_da

log_i3 = 17.3 - 0.5*np.log(ion_3['dens']/10) + 1.5*np.log(ion_3['temp']) # dens_1 already in 1e19m^-3 and temp_1 already in keV
tau_i3 = 6.6e17 * np.sqrt(2) * (ion_3['temp'])**1.5 / (ion_3['dens']*1e19*log_i3) # for singly charged deuterium ions => mi/mp=2 and Z=1 
nu_i3 = 1./tau_i3
ion_3['vnewk'] = nu_i3*amin/vth

f.write('Species Parameters: \n')
f.write('------------------- \n')
f.write('Electron: \n')
for name, value in sorted(electron.items()):
  f.write(name + ' = ' + str(value) + '\n')
f.write('\n')
f.write('Deuterium (Reference): \n')
for name, value in sorted(ion_1.items()):
  f.write(name + ' = ' + str(value) + '\n')
f.write('\n')
f.write('Hydrogen: \n')
for name, value in sorted(ion_2.items()):
  f.write(name + ' = ' + str(value) + '\n')
f.write('\n')

###################################
#Calculate geoemetry parameters #
###################################
geo = {}
geo['rhoc'] = rho_miller
geo['qinp'] = np.interp(output_radius, xb, q)
geo['shat'] =  ((q[radb_idx+1]-q[radb_idx-1])/(xb[radb_idx+1]-xb[radb_idx-1]))*(rho_miller/q[radb_idx])*dpsi_da
geo['s_hat_input'] = ((q[radb_idx+1]-q[radb_idx-1])/(xb[radb_idx+1]-xb[radb_idx-1]))*(rho_miller/q[radb_idx])*dpsi_da
geo['shift'] = (flux_centres[rad_idx+1]/100/amin-flux_centres[rad_idx-1]/100/amin)/(x[rad_idx+1]-x[rad_idx-1])*dpsi_da
geo['akappa'] = np.interp(output_radius, x, elongation)
geo['akappri'] = (elongation[rad_idx+1]-elongation[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])*dpsi_da
geo['tri'] = np.interp(output_radius, x, triang)
geo['tripri'] = (triang[rad_idx+1]-triang[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])*dpsi_da
geo['rmaj'] = rmaj[mag_axis_idx]/100+amin
geo['r_geo'] = flux_centres[-1]/100+amin

f.write('Geometry Parameters: \n')
f.write('-------------------- \n')
for name, value in sorted(geo.items()):
  f.write(name + ' = ' + str(value) + '\n')
f.write('\n')

f.write('Miscellaneous Parameters (consistent with use of Miller parameters): \n')
f.write('-------------------------------------------------------------------- \n')
f.write('irho = 2 \n')
f.write('iflux = 0 \n')
f.write('bishop = 4 \n')
f.write('local_eq = .true. \n')
