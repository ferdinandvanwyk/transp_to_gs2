import os
import sys
import operator
import warnings

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from netCDF4 import Dataset
import f90nml as nml


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
if output_radius > 1 or output_radius < 0:
    raise ValueError('Please pick and output radius in the range [0,1].')
output_time = float(sys.argv[3])
try:
    template = str(sys.argv[4])
except:
    template = None
ncfile = Dataset(in_file, 'r', format='NETCDF3')

# Create directory for plot checks and clear if already exists
plot_dir = 'plot_checks_rho_' + str(output_radius) + '_time_' + str(output_time) 
os.system('mkdir ' + plot_dir)
os.system('rm -f ' + plot_dir + '/*')

####################
# Read TRANSP file #
####################
time = ncfile.variables['TIME'][:]

# Convert input time to an index
if output_time < np.min(time) or output_time >  np.max(time):
    raise ValueError('Please pick a output time in the range [{:.3f},{:.3f}]'.format(np.min(time), np.max(time)))
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
e = 1.60217657e-19

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

##################################
# Create new namelist dictionary #
##################################

if template:
    gs2 = nml.read(template)
else:
    gs2 = {}
 
gs2_keys = ['dist_fn_knobs', 'parameters', 'theta_grid_parameters',
            'theta_grid_eik_knobs', 'species_parameters_1',
            'species_parameters_2', 'species_parameters_3',
            'species_parameters_4', 'species_parameters_5', 
            'miscellaneous']

for k in gs2_keys:
    if k not in gs2.keys():
        gs2[k] = {}

####################################
# Calculate equilibrium parameters #
####################################
gs2['miscellaneous']['amin'] = (rmaj[-1] - rmaj[0])/2/100
amin = gs2['miscellaneous']['amin']
vth = np.sqrt((2*np.interp(output_radius, x, ti)*boltz_jk/boltz_evk)/ \
              (2*proton_mass))  # T(eV)
gs2['miscellaneous']['omega'] = np.interp(output_radius, x, omega)  # rad/s
btor = fbtx*btx
b = fbx*btx
gs2['miscellaneous']['bref'] = btor[mag_axis_idx]
# n_ref(1e19 m^-3), T_ref(keV):
beta = 403.0*nd*1e6/1e19*ti/1000/(1e5*gs2['miscellaneous']['bref']**2)
if h_spec_bool:
    beta_full = 403.0*(nd*ti + nh*ti + nimp*timp + nb*tb + ne*te)* \
                    1e6/1e19/1000/(1e5*gs2['miscellaneous']['bref']**2)
else:
    beta_full = 403.0*(nd*ti + nimp*timp + nb*tb + ne*te)*1e6/1e19/1000/ \
                    (1e5*gs2['miscellaneous']['bref']**2)
gs2['parameters']['beta'] = np.interp(output_radius, x, beta)
gs2['parameters']['zeff'] = np.interp(output_radius, x, zeffp)
# See wiki definition: not taking into account B_T variation
gs2['theta_grid_eik_knobs']['beta_prime_input'] = (beta_full[rad_idx+1]-beta_full[rad_idx-1])/ \
                            (x[rad_idx+1]-x[rad_idx-1])*dpsi_da
gs2['dist_fn_knobs']['g_exb'] = (omega[rad_idx+1]-omega[rad_idx-1])/ \
                    (x[rad_idx+1]-x[rad_idx-1])*(rho_miller/q[radb_idx])* \
                    (amin/vth)*dpsi_da  # q defined on xb grid
gs2['miscellaneous']['dpsi_da'] = dpsi_da
gs2['miscellaneous']['psi_tor'] = rho_transp
gs2['miscellaneous']['psi_pol'] = np.sqrt((psi_p[flux_idx - mag_axis_idx] - 
                                      psi_p[0])/(psi_p[-1] - psi_p[0]))
gs2['miscellaneous']['bpol_flux_tube'] = np.interp(output_radius, x, bpol)
gs2['miscellaneous']['btor_flux_tube'] = btor[flux_idx]
gs2['miscellaneous']['mach'] = amin * gs2['miscellaneous']['omega'] / vth

################################
# Calculate species parameters #
################################

############
# Main ION #
############
gs2['species_parameters_1']['dens'] = 1.0
n_ref = np.interp(output_radius, x, nd)*1e6/1e19  # 1e19m^-3
gs2['species_parameters_1']['mass'] = 1.0
gs2['species_parameters_1']['temp'] = 1.0
t_ref = np.interp(output_radius, x, ti)/1000  # keV
gs2['species_parameters_1']['fprim'] = -(nd[rad_idx+1]-nd[rad_idx-1])/ \
                                         (x[rad_idx+1]-x[rad_idx-1])/ \
                                         nd[rad_idx]*dpsi_da
gs2['species_parameters_1']['tprim'] = -(ti[rad_idx+1]-ti[rad_idx-1])/ \
                                         (x[rad_idx+1]-x[rad_idx-1])/ \
                                         ti[rad_idx]*dpsi_da
gs2['species_parameters_1']['type'] = 'ion'
gs2['species_parameters_1']['z'] = 1

############
# ELECTRON #
############
electron = {}
gs2['species_parameters_2']['dens'] = np.interp(output_radius, x, ne)*1e6/1e19/n_ref  # 1e19m^-3
electron_dens = np.interp(output_radius, x, ne)*1e6/1e19  # 1e19m^-3
gs2['species_parameters_2']['mass'] = 1.0/(2.0*1836.0)  # Assume D-gs2['species_parameters_2'] plasma
gs2['species_parameters_2']['temp'] = np.interp(output_radius, x, te)/1000/t_ref  # keV
electron_temp = np.interp(output_radius, x, te)/1000  # keV
gs2['species_parameters_2']['fprim'] = -(ne[rad_idx+1]-ne[rad_idx-1])/ \
                        (x[rad_idx+1]-x[rad_idx-1])/ne[rad_idx]*dpsi_da
gs2['species_parameters_2']['tprim'] = -(te[rad_idx+1]-te[rad_idx-1])/ \
                        (x[rad_idx+1]-x[rad_idx-1])/te[rad_idx]*dpsi_da
gs2['species_parameters_2']['type'] = 'electron'
gs2['species_parameters_2']['z'] = -1

###############
# Collisions  #
###############

# See README for details of collision frequencies
loglam = 24.0 - np.log(1e4*np.sqrt(0.1*n_ref)/electron_temp)
mi = 2
zi = 1
gs2['species_parameters_2']['vnewk'] = 3.95e-3*amin*np.sqrt(0.5*mi)*loglam*n_ref/ \
                        (np.sqrt(t_ref)*electron_temp**1.5)
gs2['species_parameters_1']['vnewk'] = 9.21e-5*amin*zi**4/np.sqrt(2.)*loglam*n_ref/t_ref**2

#####################
# Other ION species #
#####################
if h_spec_bool:
    # ION SPECIES 2 - Hydrogen
    gs2['species_parameters_3']['dens_3'] = np.interp(output_radius, x, nh)*1e6/1e19/n_ref  # 1e19m^-3
    ion_2_dens = np.interp(output_radius, x, nh)*1e6/1e19  # 1e19m^-3
    gs2['species_parameters_3']['mass_3'] = 0.5
    gs2['species_parameters_3']['temp_3'] = np.interp(output_radius, x, ti)/1000/t_ref  # keV
    ion_2_temp = np.interp(output_radius, x, ti)/1000  # keV
    gs2['species_parameters_3']['fprim_3'] = -(nh[rad_idx+1]-nh[rad_idx-1])/ \
                      (x[rad_idx+1]-x[rad_idx-1])/nh[rad_idx]*dpsi_da
    gs2['species_parameters_3']['tprim_3'] = -(ti[rad_idx+1]-ti[rad_idx-1])/ \
                      (x[rad_idx+1]-x[rad_idx-1])/ti[rad_idx]*dpsi_da
    zi = 1
    gs2['species_parameters_3']['z'] = 1
    gs2['species_parameters_3']['vnewk_3'] = 9.21e-5*amin*zi**4/np.sqrt(2.)*loglam*ion_2_dens/ \
                     ion_2_temp**2
    gs2['species_parameters_3']['type'] = 'ion'

# ION SPECIES 3 - Impurity
gs2['species_parameters_4']['dens_4'] = np.interp(output_radius, x, nimp)*1e6/1e19/n_ref  # 1e19m^-3
ion_3_dens = np.interp(output_radius, x, nimp)*1e6/1e19  # 1e19m^-3
gs2['species_parameters_4']['mass_4'] = 6.0
gs2['species_parameters_4']['temp_4'] = np.interp(output_radius, x, timp)/1000/t_ref  # keV
ion_3_temp = np.interp(output_radius, x, timp)/1000  # keV
gs2['species_parameters_4']['fprim_4'] = -(nimp[rad_idx+1]-nimp[rad_idx-1])/ \
                     (x[rad_idx+1]-x[rad_idx-1])/nimp[rad_idx]*dpsi_da
gs2['species_parameters_4']['tprim_4'] = -(timp[rad_idx+1]-timp[rad_idx-1])/ \
                     (x[rad_idx+1]-x[rad_idx-1])/timp[rad_idx]*dpsi_da
zi = 6
gs2['species_parameters_4']['z'] = 6
gs2['species_parameters_4']['vnewk_4'] = 9.21e-5*amin*zi**4/np.sqrt(2.)*loglam*ion_3_dens/ion_3_temp**2
gs2['species_parameters_4']['type'] = 'ion'

# FAST ION SPECIES
gs2['species_parameters_5']['dens_5'] = np.interp(output_radius, x, nb)*1e6/1e19/n_ref  # 1e19m^-3
ion_4_dens = np.interp(output_radius, x, nb)*1e6/1e19  # 1e19m^-3
gs2['species_parameters_5']['mass_5'] = 1.0
gs2['species_parameters_5']['temp_5'] = np.interp(output_radius, x, tb)/1000/t_ref  # keV
ion_4_temp = np.interp(output_radius, x, tb)/1000  # keV
gs2['species_parameters_5']['fprim_5'] = -(nb[rad_idx+1]-nb[rad_idx-1])/ \
                     (x[rad_idx+1]-x[rad_idx-1])/nb[rad_idx]*dpsi_da
gs2['species_parameters_5']['tprim_5'] = -(tb[rad_idx+1]-tb[rad_idx-1])/ \
                     (x[rad_idx+1]-x[rad_idx-1])/tb[rad_idx]*dpsi_da
zi = 1
gs2['species_parameters_5']['z'] = 1
gs2['species_parameters_5']['vnewk_5'] = 9.21e-5*amin*zi**4/np.sqrt(2.)*loglam*ion_4_dens/ion_4_temp**2
gs2['species_parameters_5']['type'] = 'beam'

##################################
# Calculate geoemetry parameters #
##################################
gs2['theta_grid_parameters']['rhoc'] = rho_miller
gs2['theta_grid_parameters']['qinp'] = np.interp(output_radius, xb, q)
gs2['theta_grid_parameters']['shat'] =  ((q[radb_idx+1]-q[radb_idx-1])/(xb[radb_idx+1]-xb[radb_idx-1]))* \
                   (rho_miller/q[radb_idx])*dpsi_da
gs2['theta_grid_parameters']['s_hat_input'] = ((q[radb_idx+1] - q[radb_idx-1])/ \
                      (xb[radb_idx+1] - xb[radb_idx-1]))* \
                     (rho_miller/q[radb_idx])*dpsi_da
gs2['theta_grid_parameters']['shift'] = (flux_centres[rad_idx+1]/100/amin-flux_centres[rad_idx-1]/100/amin)/ \
                   (x[rad_idx+1]-x[rad_idx-1])*dpsi_da
gs2['theta_grid_parameters']['akappa'] = np.interp(output_radius, x, elongation)
gs2['theta_grid_parameters']['akappri'] = (elongation[rad_idx+1]-elongation[rad_idx-1])/ \
                     (x[rad_idx+1]-x[rad_idx-1])*dpsi_da
gs2['theta_grid_parameters']['tri'] = np.arcsin(np.interp(output_radius, x, triang))
gs2['theta_grid_parameters']['tripri'] = (np.arcsin(triang[rad_idx+1])-np.arcsin(triang[rad_idx-1]))/ \
                    (x[rad_idx+1]-x[rad_idx-1])*dpsi_da
gs2['theta_grid_parameters']['rmaj'] = (rmaj[flux_idx] + rmaj[mag_axis_idx - (flux_idx-mag_axis_idx)])/100/2/amin
gs2['theta_grid_parameters']['r_geo'] = (rmaj[-1] + rmaj[0])/2/100/amin

gs2['theta_grid_eik_knobs']['irho'] = 2
gs2['theta_grid_eik_knobs']['iflux'] = 0
gs2['theta_grid_eik_knobs']['bishop'] = 4
gs2['theta_grid_eik_knobs']['local_eq'] = '.true.'

#################################################
# Calculate reference values for normalizations #
#################################################

gs2['miscellaneous']['vref'] = vth
gs2['miscellaneous']['lref'] = amin
gs2['miscellaneous']['nref'] = n_ref*1e19
gs2['miscellaneous']['tref'] = t_ref*1000/boltz_evk
larmor_freq_ref = e * gs2['miscellaneous']['bref'] / (2*proton_mass)
gs2['miscellaneous']['rhoref'] = vth / larmor_freq_ref

######################
# Write the namelist #
######################

outfile_name = 'gs2_rho_' + str(output_radius) + '_time_' + str(output_time) + '.in' 
for k1 in gs2.keys():
    for k2 in gs2[k1].keys():
        if type(gs2[k1][k2]) == np.float32 or type(gs2[k1][k2]) == np.float64:
            gs2[k1][k2] = float(gs2[k1][k2])
            
gs2 = nml.Namelist(gs2)

if outfile_name in os.listdir():
    os.system('rm -rf ' + outfile_name)

gs2.write(outfile_name)

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

plt.plot(amin*R(gs2['theta_grid_parameters']['rmaj'], rho_miller, theta, gs2['theta_grid_parameters']['tri']),
         amin*Z(gs2['theta_grid_parameters']['akappa'], rho_miller, theta))
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
