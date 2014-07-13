transp_read
===========

Reads a TRANSP file and outputs parameters relevant to a GS2 simulation

Input:
-----
* Takes the name of the TRANSP file as first command line input
* Second command line input is radial location of output (will interpolate)
* Third command line input is time for output (will not interpolate)
* Input command therefore looks like: python main.py (filename) (radius) (time)

Output:
-------
* Produces file called gs2.in
* Equilibrium parameters, e.g. temperatures, densities, gradients, flows, shears
* Geometric Miller parameters
* Miscellaneous parameters which need to be set when using Miller parameters (only option for this code so far).

Requirements:
-------------

* Python (at least 2.7)
* NetCDF4 
* HDF5
* Python NetCDF4 Interface: https://pypi.python.org/pypi/netCDF4
* Numpy + SciPy + Matplotlib

Caveats
-------

* The sign of the flow shear may not be correct! This program simply outputs the gradient of the flow without changing sign. The real sign of the flow shear will depend on the sign convention used in the experiment. <Reference something here> 
* Use python-netcdf4 interface since it seems to be the only thing that reads the TRANSP file. Could try SciPy NetCDF functions as documented here: http://docs.scipy.org/doc/scipy/reference/generated/scipy.io.netcdf.netcdf_file.html 
* Electron and ion densities and gradients are simply read and outputted. These may be different to each other due to multiple ion species being present in the experiment. In a two species GS2 simulation however, must set ni = ne and fprim_1 = fprim_2 to satisy quasineutrality.

Plotting Routine
----------------
A rudimentarty plotting routine is also included and works in the following way:

* It can be called as: python plot.py (independent variable) (dependent variable) (index of other dimension)
* Independent variables are e.g. TIME, TIME3, X, XB, RMAJM
* Dependent variable are a function of one or two of these variables.
* Most variables are a function of two independent variables so the standard behaviour requires specifying the index of the other independent variable. For example, if you wanted to plot Q(TIME, XB) as a function of time at the radial location idx = 50: 'python main.py TIME Q 50'. However if you wanted to plot Q as a function of radius for a given time index (= 98): 'python main.py XB Q 98'

Collisions
----------

*Equations for collision times/frequencies taken from Wesson Section 2.15 (T(keV), n(m^-3), mi=ion mass, mp=proton mass):
  * nu_i = 1/tau_i = 1 / [6.60e17 ((mi/mp)^1/2 Ti^3/2) / (n Z^4 ln(Lambda))] s^-1
  * nu_e = 1/tau_e = 1 / [1.09e16 (Te^3/2) / (n Z^2 ln(Lambda))] s^-1
  * log(Lambda) = Coulomb logarithm
    *log_i = 17.3 - 0.5*np.log(n/1e20) + 1.5*np.log(ti) 
    *log_e = 14.9 - 0.5*np.log(n/1e20) + np.log(te) 
  * lambda = debye length = 2.35e5 (T/n) m
*For the purposes of collision frequency calculations, ni = ne = n to satisfy quasineutrality. 





