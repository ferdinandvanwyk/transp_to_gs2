import os, sys
import operator
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from netCDF4 import Dataset

#function that computes the debye length given T(keV), n(m^-3)
#Follows NRL plasma formulary
def debye(temp, dens):
  return 2.35e5*np.sqrt(temp)/np.sqrt(dens)

def plot_dash(x, y, x0, filename):
    """
    Plot the profile with a dashed line at the radial location of interest.

    Parameters
    ----------
    x : array_like
        Independent variable
    y : array_like
        Dependent variable
    x0 : float
        Position along x to plot the dashed line
    filename : str
        Name of saved plot
    """

    plt.clf()
    fig, ax = plt.subplots()
    plt.plot(x, y)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    plt.axvline(x0, color='r', ls='--')
    ax.grid(b=True, which='major')
    ax.grid(b=True, which='minor')
    plt.savefig('plot_checks/' + filename)
    plt.close(fig)

def plot_gradient(x, y, x0, filename):
    """
    Plot the profile with a dashed line at the radial location of interest.

    Parameters
    ----------
    x : array_like
        Independent variable
    y : array_like
        Dependent variable
    x0 : float
        Position along x to plot the dashed line
    filename : str
        Name of saved plot
    """

    y_grad = np.gradient(y, x[1] - x[0])
    slope = np.interp(x0, x, y_grad)
    y0 = np.interp(x0, x, y)

    plt.clf()
    fig, ax = plt.subplots()
    plt.plot(x, y)
    plt.plot(x, slope*x + (y0 - slope*x0), 'r--')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(b=True, which='major')
    ax.grid(b=True, which='minor')
    plt.savefig('plot_checks/' + filename)
    plt.close(fig)

in_file =  str(sys.argv[1])
output_radius =  float(sys.argv[2])
output_time =  float(sys.argv[3])
ncfile = Dataset(in_file, 'r', format='NETCDF3')

# Create directory for plot checks and clear if already exists
os.system('mkdir plot_checks')
os.system('rm  -f plot_checks/*')

#Read relevant data from TRANSP file
time = ncfile.variables['TIME'][:]

#Convert input time to an index
t_idx, min_value = min(enumerate(abs(time - output_time)),
                       key=operator.itemgetter(1))

x = ncfile.variables['X'][t_idx,:]
xb = ncfile.variables['XB'][t_idx,:]
ni_tot = ncfile.variables['NI'][t_idx,:]
nd = ncfile.variables['ND'][t_idx,:] #Deuterium (reference species) density

try:
    nh = ncfile.variables['NH'][t_idx,:] #Hydrogen density
    h_spec_bool = True
except KeyError:
    warnings.warn('No H species detected. Parameters will adjust accordingly.')
    h_spec_bool = False

nimp = ncfile.variables['NIMP'][t_idx,:] #Impurity ion density
nb = ncfile.variables['BDENS'][t_idx,:] #Beam ion density
ne = ncfile.variables['NE'][t_idx,:]
ti = ncfile.variables['TI'][t_idx,:] #Combined D, H
timp = ncfile.variables['TX'][t_idx,:] #Impurity temp
tb = ncfile.variables['EBEAM_D'][t_idx,:] #Energy/temp of the beam ions
te = ncfile.variables['TE'][t_idx,:]
fbtx = ncfile.variables['FBTX'][t_idx,:]
fbx = ncfile.variables['FBX'][t_idx,:]
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

# Normalization quantities
boltz_jk = 1.3806488e-23
boltz_evk = 8.6173324e-5
proton_mass = 1.672621777e-27

# Need to calculate factor which relates gradients in TRANSP psi_n
# (=sqrt(psi_t/psi_tLCFS)) and Miller a_n (= diameter/diameter LCFS)
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
vth = np.sqrt((2*np.interp(output_radius, x, ti)*boltz_jk/boltz_evk)/(2*proton_mass)) # T(eV)
equil['omega'] = np.interp(output_radius, x, omega) #rad/s
btor = fbtx*btx
b = fbx*btx
equil['bref'] = np.interp(flux_centres[-1], rmaj, btor)
beta =  403.0*nd*1e6/1e19*ti/1000/(1e5*equil['bref']**2) #n_ref(1e19 m^-3), T_ref(keV)
if h_spec_bool:
    beta_full =  403.0*(nd*ti + nh*ti + nimp*timp + nb*tb + ne*te)*1e6/1e19/1000/(1e5*equil['bref']**2) #n_ref(1e19 m^-3), T_ref(keV)
else:
    beta_full =  403.0*(nd*ti + nimp*timp + nb*tb + ne*te)*1e6/1e19/1000/(1e5*equil['bref']**2) #n_ref(1e19 m^-3), T_ref(keV)
equil['beta'] = np.interp(output_radius, x, beta)
equil['zeff'] = np.interp(output_radius, x, zeffp)
equil['beta_prime_input'] = (beta_full[rad_idx+1]-beta_full[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])*dpsi_da #See wiki definition: not taking into account B_T variation
equil['g_exb'] = (omega[rad_idx+1]-omega[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])*(rho_miller/q[radb_idx])*(amin/vth)*dpsi_da #q defined on xb grid
equil['dpsi_da'] = dpsi_da
equil['bpol_flux_tube'] = np.interp(output_radius, x, bpol)
equil['btor_flux_tube'] = btor[flux_idx]

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
# Main ION #
############
#ION SPECIES 1 - Deuterium (reference)
ion_1 = {}
ion_1['dens'] = 1.0
n_ref = np.interp(output_radius, x, nd)*1e6/1e19 #1e19m^-3
ion_1['mass'] = 1.0
ion_1['temp'] = 1.0
t_ref = np.interp(output_radius, x, ti)/1000 #keV
ion_1['fprim'] = -(nd[rad_idx+1]-nd[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/nd[rad_idx]*dpsi_da
ion_1['tprim'] = -(ti[rad_idx+1]-ti[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/ti[rad_idx]*dpsi_da

############
# ELECTRON #
############
electron = {}
electron['dens'] = np.interp(output_radius, x, ne)*1e6/1e19/n_ref #1e19m^-3
electron_dens = np.interp(output_radius, x, ne)*1e6/1e19 #1e19m^-3
electron['mass'] = 1.0/(2.0*1836.0) #Assume D-electron plasma
electron['temp'] = np.interp(output_radius, x, te)/1000/t_ref #keV
electron_temp = np.interp(output_radius, x, te)/1000 #keV
electron['fprim'] = -(ne[rad_idx+1]-ne[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/ne[rad_idx]*dpsi_da
electron['tprim'] = -(te[rad_idx+1]-te[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/te[rad_idx]*dpsi_da

###############
# Collisions  #
###############

#See README for details of collision frequencies
loglam = 24.0 - np.log(1e4*np.sqrt(0.1*n_ref)/electron_temp)
mi = 2
zi = 1
electron['vnewk'] = 3.95e-3*amin*np.sqrt(0.5*mi)*loglam*n_ref/(np.sqrt(t_ref)*electron_temp**1.5)
ion_1['vnewk'] = 9.21e-5*amin*zi**4/np.sqrt(2.)*loglam*n_ref/t_ref**2

#####################
# Other ION species #
#####################
if h_spec_bool:
    #ION SPECIES 2 - Hydrogen
    ion_2 = {}
    ion_2['dens'] = np.interp(output_radius, x, nh)*1e6/1e19/n_ref #1e19m^-3
    ion_2_dens = np.interp(output_radius, x, nh)*1e6/1e19 #1e19m^-3
    ion_2['mass'] = 0.5
    ion_2['temp'] = np.interp(output_radius, x, ti)/1000/t_ref #keV
    ion_2_temp = np.interp(output_radius, x, ti)/1000 #keV
    ion_2['fprim'] = -(nh[rad_idx+1]-nh[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/nh[rad_idx]*dpsi_da
    ion_2['tprim'] = -(ti[rad_idx+1]-ti[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/ti[rad_idx]*dpsi_da
    zi = 1
    ion_2['vnewk'] = 9.21e-5*amin*zi**4/np.sqrt(2.)*loglam*ion_2_dens/ion_2_temp**2

#ION SPECIES 3 - Impurity
ion_3 = {}
ion_3['dens'] = np.interp(output_radius, x, nimp)*1e6/1e19/n_ref #1e19m^-3
ion_3_dens = np.interp(output_radius, x, nimp)*1e6/1e19 #1e19m^-3
ion_3['mass'] = 6.0
ion_3['temp'] = np.interp(output_radius, x, timp)/1000/t_ref #keV
ion_3_temp = np.interp(output_radius, x, timp)/1000 #keV
ion_3['fprim'] = -(nimp[rad_idx+1]-nimp[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/nimp[rad_idx]*dpsi_da
ion_3['tprim'] = -(timp[rad_idx+1]-timp[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/timp[rad_idx]*dpsi_da
zi = 6
ion_3['vnewk'] = 9.21e-5*amin*zi**4/np.sqrt(2.)*loglam*ion_3_dens/ion_3_temp**2

#FAST ION SPECIES
ion_4 = {}
ion_4['dens'] = np.interp(output_radius, x, nb)*1e6/1e19/n_ref #1e19m^-3
ion_4_dens = np.interp(output_radius, x, nb)*1e6/1e19 #1e19m^-3
ion_4['mass'] = 1.0
ion_4['temp'] = np.interp(output_radius, x, tb)/1000/t_ref #keV
ion_4_temp = np.interp(output_radius, x, tb)/1000 #keV
ion_4['fprim'] = -(nb[rad_idx+1]-nb[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/nb[rad_idx]*dpsi_da
ion_4['tprim'] = -(tb[rad_idx+1]-tb[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/tb[rad_idx]*dpsi_da
zi = 1
ion_4['vnewk'] = 9.21e-5*amin*zi**4/np.sqrt(2.)*loglam*ion_4_dens/ion_4_temp**2

f.write('Species Parameters: \n')
f.write('------------------- \n')
f.write('Deuterium (Reference): \n')
for name, value in sorted(ion_1.items()):
  f.write(name + ' = ' + str(value) + '\n')
f.write('\n')
f.write('Electron: \n')
for name, value in sorted(electron.items()):
  f.write(name + ' = ' + str(value) + '\n')
f.write('\n')
if h_spec_bool:
    f.write('Hydrogen: \n')
    for name, value in sorted(ion_2.items()):
      f.write(name + ' = ' + str(value) + '\n')
    f.write('\n')
f.write('Impurity: \n')
for name, value in sorted(ion_3.items()):
  f.write(name + ' = ' + str(value) + '\n')
f.write('\n')
f.write('Beam Ions: \n')
for name, value in sorted(ion_4.items()):
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
geo['tri'] = np.arcsin(np.interp(output_radius, x, triang))
geo['tripri'] = (np.arcsin(triang[rad_idx+1])-np.arcsin(triang[rad_idx-1]))/(x[rad_idx+1]-x[rad_idx-1])*dpsi_da
geo['rmaj'] = rmaj[mag_axis_idx]/100/amin
geo['r_geo'] = flux_centres[-1]/100/amin

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
f.write('\n')

#################################################
# Calculate reference values for normalizations #
#################################################
e = 1.60217657e-19
ref = {}
ref['vref'] = vth
ref['lref'] = amin
ref['nref'] = n_ref*1e19
ref['tref'] = t_ref*1000/boltz_evk
ref['bref'] = equil['bref']
larmor_freq_ref = e * ref['bref'] / (2*proton_mass)
ref['rhoref'] = vth / larmor_freq_ref

f.write('Reference Parameters: \n')
f.write('--------------------- \n')
for name, value in sorted(ref.items()):
  f.write(name + ' = ' + str(value) + '\n')
f.write('\n')

####################
# Diagnostic plots #
####################

plot_dash(x, omega, output_radius, 'omega.png')
plot_dash(rmaj, btor, rmaj[mag_axis_idx], 'bref.png')
plot_dash(x, beta, output_radius, 'beta_ref.png')
plot_dash(x, beta_full, output_radius, 'beta.png')
plot_dash(x, ti, output_radius, 'ti.png')
plot_dash(x, te, output_radius, 'te.png')
plot_dash(x, nd, output_radius, 'n_D.png')
plot_dash(x, nb, output_radius, 'n_beam.png')
plot_dash(x, nimp, output_radius, 'n_impurity.png')
plot_dash(x, ne, output_radius, 'ne.png')
plot_dash(xb, q, output_radius, 'q.png')
plot_dash(xb, elongation, output_radius, 'akappa.png')
plot_dash(xb, triang, output_radius, 'tri.png')

plot_gradient(x, beta_full, output_radius, 'beta_prime.png')
plot_gradient(x, omega, output_radius, 'g_exb.png')
plot_gradient(x, ti, output_radius, 'tprim_1.png')
plot_gradient(x, te, output_radius, 'tprim_2.png')
plot_gradient(x, nd, output_radius, 'fprim_D.png')
plot_gradient(x, nb, output_radius, 'fprim_beam.png')
plot_gradient(x, nimp, output_radius, 'fprim_impurity.png')
plot_gradient(x, ne, output_radius, 'fprim_2.png')
plot_gradient(xb, q, output_radius, 'shat.png')
plot_gradient(xb, elongation, output_radius, 'akappri.png')
plot_gradient(xb, triang, output_radius, 'tripri.png')

if h_spec_bool:
    plot_dash(x, nh, output_radius, 'n_H.png')
    plot_gradient(x, nh, output_radius, 'fprim_H.png')
