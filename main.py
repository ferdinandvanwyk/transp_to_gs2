import os
import sys
import operator
import warnings

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from netCDF4 import Dataset


def debye(temp, dens):
    """
    Function that computes the debye length given T(keV), n(m^-3). Follows
    NRL plasma formulary.
    """
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
    plt.savefig(plot_dir + '/' + filename)
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
    plt.savefig(plot_dir + '/' + filename)
    plt.close(fig)


def R(rmaj, rhoc, th, tri):
        return rmaj + rhoc * np.cos(th + tri*np.sin(th))


def Z(akappa, rhoc, th):
        return akappa * rhoc*np.sin(th)

in_file = str(sys.argv[1])
output_radius = float(sys.argv[2])
output_time = float(sys.argv[3])
ncfile = Dataset(in_file, 'r', format='NETCDF3')

# Create directory for plot checks and clear if already exists
plot_dir = 'plot_checks_rho_' + str(output_radius) + '_time_' + str(output_time) 
os.system('mkdir ' + plot_dir)
os.system('rm -f ' + plot_dir + '/*')

# Read relevant data from TRANSP file
time = ncfile.variables['TIME'][:]

# Convert input time to an index
t_idx, min_value = min(enumerate(abs(time - output_time)),
                       key=operator.itemgetter(1))

x = ncfile.variables['X'][t_idx, :]
xb = ncfile.variables['XB'][t_idx, :]
ni_tot = ncfile.variables['NI'][t_idx, :]
nd = ncfile.variables['ND'][t_idx, :]  # Deuterium (reference species) density

try:
    nh = ncfile.variables['NH'][t_idx, :]  # Hydrogen density
    h_spec_bool = True
except KeyError:
    warnings.warn('No H species detected. Parameters will adjust accordingly.')
    h_spec_bool = False

nimp = ncfile.variables['NIMP'][t_idx, :]  # Impurity ion density
nb = ncfile.variables['BDENS'][t_idx, :]  # Beam ion density
ne = ncfile.variables['NE'][t_idx, :]
ti = ncfile.variables['TI'][t_idx, :]  # Combined D, H
timp = ncfile.variables['TX'][t_idx, :]  # Impurity temp
tb = ncfile.variables['EBEAM_D'][t_idx, :]  # Energy/temp of the beam ions
te = ncfile.variables['TE'][t_idx, :]
fbtx = ncfile.variables['FBTX'][t_idx, :]
fbx = ncfile.variables['FBX'][t_idx, :]
bpbt = ncfile.variables['FBPBT'][t_idx, :]
btx = ncfile.variables['BTX'][t_idx, :]
rmaj = ncfile.variables['RMAJM'][t_idx, :]
bpol = ncfile.variables['BPOL'][t_idx, :]
omega = ncfile.variables['OMEGA'][t_idx, :]
q = ncfile.variables['Q'][t_idx, :]
shat = ncfile.variables['SHAT'][t_idx, :]
flux_centres = ncfile.variables['RMJMP'][t_idx, :]
elongation = ncfile.variables['ELONG'][t_idx, :]
triang = ncfile.variables['TRIANG'][t_idx, :]
zeffp = ncfile.variables['ZEFFP'][t_idx, :]
psi_t = ncfile.variables['TRFMP'][t_idx, :]
psi_p = ncfile.variables['PLFLX'][t_idx, :]

# Normalization quantities
boltz_jk = 1.3806488e-23
boltz_evk = 8.6173324e-5
proton_mass = 1.672621777e-27

# Need to calculate factor which relates gradients in TRANSP psi_n
# (=sqrt(psi_t/psi_tLCFS)) and Miller a_n (= diameter/diameter LCFS)
flux_rmaj = np.interp(output_radius, np.linspace(-1, 1, rmaj.shape[0]), rmaj)
flux_idx, min_value = min(enumerate(abs(rmaj - flux_rmaj)),
                          key=operator.itemgetter(1))
mag_axis_idx, min_value = min(enumerate(abs(bpbt)), key=operator.itemgetter(1))
a_left = (rmaj[flux_idx-1] - rmaj[mag_axis_idx - (flux_idx-1-mag_axis_idx)])/ \
         (rmaj[-1] - rmaj[0])  # diameter/diameter of LCFS
rho_miller = (rmaj[flux_idx] - rmaj[mag_axis_idx - (flux_idx-mag_axis_idx)])/ \
             (rmaj[-1] - rmaj[0])  # diameter/diameter of LCFS
a_right = (rmaj[flux_idx+1] - rmaj[mag_axis_idx - (flux_idx+1-mag_axis_idx)])/\
          (rmaj[-1] - rmaj[0])  # diameter/diameter of LCFS
psi_left = np.sqrt(psi_t[flux_idx-1]/psi_t[-1])  # sqrt(psi_t/psi_LCFS)
rho_transp = np.sqrt(psi_t[flux_idx]/psi_t[-1])  # sqrt(psi_t/psi_LCFS)
psi_right = np.sqrt(psi_t[flux_idx+1]/psi_t[-1])  # sqrt(psi_t/psi_LCFS)
# Coefficient which relates psi_n and a_n grids
dpsi_da = (psi_right-psi_left)/(a_right-a_left)

# Gradient calculation requires numerical differentiation
# Use a second order central method: f'(x) = [f(x+h) - f(x-h)]/2h
# Find points closest points either side of output_radius and use those
rad_idx, min_value = min(enumerate(abs(x - output_radius)),
                         key=operator.itemgetter(1))
radb_idx, min_value = min(enumerate(abs(xb - output_radius)),
                          key=operator.itemgetter(1))

####################################
# Calculate equilibrium parameters #
####################################
equil = {}
equil['amin'] = (rmaj[-1] - rmaj[0])/2/100
amin = equil['amin']
vth = np.sqrt((2*np.interp(output_radius, x, ti)*boltz_jk/boltz_evk)/ \
              (2*proton_mass))  # T(eV)
equil['omega'] = np.interp(output_radius, x, omega)  # rad/s
btor = fbtx*btx
b = fbx*btx
equil['bref'] = btor[mag_axis_idx]
# n_ref(1e19 m^-3), T_ref(keV):
beta = 403.0*nd*1e6/1e19*ti/1000/(1e5*equil['bref']**2)
if h_spec_bool:
    beta_full = 403.0*(nd*ti + nh*ti + nimp*timp + nb*tb + ne*te)* \
                    1e6/1e19/1000/(1e5*equil['bref']**2)
else:
    beta_full = 403.0*(nd*ti + nimp*timp + nb*tb + ne*te)*1e6/1e19/1000/ \
                    (1e5*equil['bref']**2)
equil['beta'] = np.interp(output_radius, x, beta)
equil['zeff'] = np.interp(output_radius, x, zeffp)
# See wiki definition: not taking into account B_T variation
equil['beta_prime_input'] = (beta_full[rad_idx+1]-beta_full[rad_idx-1])/ \
                            (x[rad_idx+1]-x[rad_idx-1])*dpsi_da
equil['g_exb'] = (omega[rad_idx+1]-omega[rad_idx-1])/ \
                    (x[rad_idx+1]-x[rad_idx-1])*(rho_miller/q[radb_idx])* \
                    (amin/vth)*dpsi_da  # q defined on xb grid
equil['dpsi_da'] = dpsi_da
equil['psi_tor'] = rho_transp
equil['psi_pol'] = np.sqrt((psi_p[flux_idx - mag_axis_idx] - psi_p[0])/(psi_p[-1] - psi_p[0]))
equil['bpol_flux_tube'] = np.interp(output_radius, x, bpol)
equil['btor_flux_tube'] = btor[flux_idx]
equil['mach'] = amin * equil['omega'] / vth

outfile_name = 'gs2_rho_' + str(output_radius) + '_time_' + str(output_time) + '.in' 
f = open(outfile_name, 'w')
f.write('Equilibrium Parameters: \n')
f.write('----------------------- \n')
for name, value in sorted(equil.items()):
    f.write(name + ' = ' + str(value) + '\n')
f.write('\n')

################################
# Calculate species parameters #
################################

############
# Main ION #
############
# ION SPECIES 1 - Deuterium (reference)
ion_1 = {}
ion_1['dens_1'] = 1.0
n_ref = np.interp(output_radius, x, nd)*1e6/1e19  # 1e19m^-3
ion_1['mass_1'] = 1.0
ion_1['temp_1'] = 1.0
t_ref = np.interp(output_radius, x, ti)/1000  # keV
ion_1['fprim_1'] = -(nd[rad_idx+1]-nd[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/ \
                     nd[rad_idx]*dpsi_da
ion_1['tprim_1'] = -(ti[rad_idx+1]-ti[rad_idx-1])/(x[rad_idx+1]-x[rad_idx-1])/ \
                     ti[rad_idx]*dpsi_da

############
# ELECTRON #
############
electron = {}
electron['dens_2'] = np.interp(output_radius, x, ne)*1e6/1e19/n_ref  # 1e19m^-3
electron_dens = np.interp(output_radius, x, ne)*1e6/1e19  # 1e19m^-3
electron['mass_2'] = 1.0/(2.0*1836.0)  # Assume D-electron plasma
electron['temp_2'] = np.interp(output_radius, x, te)/1000/t_ref  # keV
electron_temp = np.interp(output_radius, x, te)/1000  # keV
electron['fprim_2'] = -(ne[rad_idx+1]-ne[rad_idx-1])/ \
                        (x[rad_idx+1]-x[rad_idx-1])/ne[rad_idx]*dpsi_da
electron['tprim_2'] = -(te[rad_idx+1]-te[rad_idx-1])/ \
                        (x[rad_idx+1]-x[rad_idx-1])/te[rad_idx]*dpsi_da

###############
# Collisions  #
###############

# See README for details of collision frequencies
loglam = 24.0 - np.log(1e4*np.sqrt(0.1*n_ref)/electron_temp)
mi = 2
zi = 1
electron['vnewk_2'] = 3.95e-3*amin*np.sqrt(0.5*mi)*loglam*n_ref/ \
                        (np.sqrt(t_ref)*electron_temp**1.5)
ion_1['vnewk_1'] = 9.21e-5*amin*zi**4/np.sqrt(2.)*loglam*n_ref/t_ref**2

#####################
# Other ION species #
#####################
if h_spec_bool:
    # ION SPECIES 2 - Hydrogen
    ion_2 = {}
    ion_2['dens_3'] = np.interp(output_radius, x, nh)*1e6/1e19/n_ref  # 1e19m^-3
    ion_2_dens = np.interp(output_radius, x, nh)*1e6/1e19  # 1e19m^-3
    ion_2['mass_3'] = 0.5
    ion_2['temp_3'] = np.interp(output_radius, x, ti)/1000/t_ref  # keV
    ion_2_temp = np.interp(output_radius, x, ti)/1000  # keV
    ion_2['fprim_3'] = -(nh[rad_idx+1]-nh[rad_idx-1])/ \
                      (x[rad_idx+1]-x[rad_idx-1])/nh[rad_idx]*dpsi_da
    ion_2['tprim_3'] = -(ti[rad_idx+1]-ti[rad_idx-1])/ \
                      (x[rad_idx+1]-x[rad_idx-1])/ti[rad_idx]*dpsi_da
    zi = 1
    ion_2['vnewk_3'] = 9.21e-5*amin*zi**4/np.sqrt(2.)*loglam*ion_2_dens/ \
                     ion_2_temp**2

# ION SPECIES 3 - Impurity
ion_3 = {}
ion_3['dens_4'] = np.interp(output_radius, x, nimp)*1e6/1e19/n_ref  # 1e19m^-3
ion_3_dens = np.interp(output_radius, x, nimp)*1e6/1e19  # 1e19m^-3
ion_3['mass_4'] = 6.0
ion_3['temp_4'] = np.interp(output_radius, x, timp)/1000/t_ref  # keV
ion_3_temp = np.interp(output_radius, x, timp)/1000  # keV
ion_3['fprim_4'] = -(nimp[rad_idx+1]-nimp[rad_idx-1])/ \
                     (x[rad_idx+1]-x[rad_idx-1])/nimp[rad_idx]*dpsi_da
ion_3['tprim_4'] = -(timp[rad_idx+1]-timp[rad_idx-1])/ \
                     (x[rad_idx+1]-x[rad_idx-1])/timp[rad_idx]*dpsi_da
zi = 6
ion_3['vnewk_4'] = 9.21e-5*amin*zi**4/np.sqrt(2.)*loglam*ion_3_dens/ion_3_temp**2

# FAST ION SPECIES
ion_4 = {}
ion_4['dens_5'] = np.interp(output_radius, x, nb)*1e6/1e19/n_ref  # 1e19m^-3
ion_4_dens = np.interp(output_radius, x, nb)*1e6/1e19  # 1e19m^-3
ion_4['mass_5'] = 1.0
ion_4['temp_5'] = np.interp(output_radius, x, tb)/1000/t_ref  # keV
ion_4_temp = np.interp(output_radius, x, tb)/1000  # keV
ion_4['fprim_5'] = -(nb[rad_idx+1]-nb[rad_idx-1])/ \
                     (x[rad_idx+1]-x[rad_idx-1])/nb[rad_idx]*dpsi_da
ion_4['tprim_5'] = -(tb[rad_idx+1]-tb[rad_idx-1])/ \
                     (x[rad_idx+1]-x[rad_idx-1])/tb[rad_idx]*dpsi_da
zi = 1
ion_4['vnewk_5'] = 9.21e-5*amin*zi**4/np.sqrt(2.)*loglam*ion_4_dens/ion_4_temp**2

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

##################################
# Calculate geoemetry parameters #
##################################
geo = {}
geo['rhoc'] = rho_miller
geo['qinp'] = np.interp(output_radius, xb, q)
geo['shat'] =  ((q[radb_idx+1]-q[radb_idx-1])/(xb[radb_idx+1]-xb[radb_idx-1]))* \
                   (rho_miller/q[radb_idx])*dpsi_da
geo['s_hat_input'] = ((q[radb_idx+1] - q[radb_idx-1])/ \
                      (xb[radb_idx+1] - xb[radb_idx-1]))* \
                     (rho_miller/q[radb_idx])*dpsi_da
geo['shift'] = (flux_centres[rad_idx+1]/100/amin-flux_centres[rad_idx-1]/100/amin)/ \
                   (x[rad_idx+1]-x[rad_idx-1])*dpsi_da
geo['akappa'] = np.interp(output_radius, x, elongation)
geo['akappri'] = (elongation[rad_idx+1]-elongation[rad_idx-1])/ \
                     (x[rad_idx+1]-x[rad_idx-1])*dpsi_da
geo['tri'] = np.arcsin(np.interp(output_radius, x, triang))
geo['tripri'] = (np.arcsin(triang[rad_idx+1])-np.arcsin(triang[rad_idx-1]))/ \
                    (x[rad_idx+1]-x[rad_idx-1])*dpsi_da
geo['rmaj'] = (rmaj[flux_idx] + rmaj[mag_axis_idx - (flux_idx-mag_axis_idx)])/100/2/amin
geo['r_geo'] = (rmaj[-1] + rmaj[0])/2/100/amin

f.write('Geometry Parameters: \n')
f.write('-------------------- \n')
for name, value in sorted(geo.items()):
    f.write(name + ' = ' + str(value) + '\n')
f.write('\n')

f.write('Miscellaneous Parameters (consistent with use of Miller parameters):\n')
f.write('--------------------------------------------------------------------\n')
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

# Poloidal plot of flux tube
theta = np.linspace(-np.pi, np.pi, 100)
fig, ax = plt.subplots(1, 1)

plt.plot(amin*R(geo['rmaj'], rho_miller, theta, geo['tri']),
         amin*Z(geo['akappa'], rho_miller, theta))
plt.xlim(0, 6*amin)
plt.title(r'Flux Surface at $\rho_c = {:.3f}$'.format(rho_miller))
ax.set_aspect('equal')
ax.xaxis.set_minor_locator(
    mpl.ticker.MultipleLocator((plt.xticks()[0][1] - \
                                plt.xticks()[0][0]) / 2.0))
ax.yaxis.set_minor_locator(
    mpl.ticker.MultipleLocator((plt.yticks()[0][1] - \
                                plt.yticks()[0][0]) / 2.0))
ax.grid(True, 'major', color='0.92', linestyle='-', linewidth=1.4)
ax.grid(True, 'minor', color='0.92', linestyle='-', linewidth=0.7)
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.set_axisbelow(True)
plt.xlabel(r'$R$ (m)')
plt.ylabel(r'$Z$ (m)')
plt.savefig(plot_dir + '/flux_tube.png')
