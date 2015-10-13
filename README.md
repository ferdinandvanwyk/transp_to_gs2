transp_to_gs2
===========

Reads a TRANSP output file and outputs parameters relevant to a GS2 simulation

Input:
-----
* Takes the name of the TRANSP file as first command line input.
* Second command line input is radial location of output (will interpolate). This will be a value between 0 and 1 specifying the value of the square root of the toroidal flux.
* Third command line input is time for output (will determine nearest time index, i.e. won't interpolate).
* Input command therefore looks like: 
``` 
$ python main.py <filename> <radius> <time>
```

Output:
-------
* Produces file called gs2.in
* Equilibrium parameters, e.g. temperatures, densities, gradients (see note 
  below), flows, shears
* Geometric Miller parameters (with appropriate corrections for GS2, see 
  http://arxiv.org/abs/1403.3293)
* Miscellaneous parameters which need to be set when using Miller parameters 
  (only option for this code so far).

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
* Electron and ion densities and temperatures are normalized to the reference species (D) values. 

Plotting Routine
----------------
A rudimentarty plotting routine is also included and works in the following way:

* It can be called as: python plot.py (filename) (independent variable) (dependent variable) (index of other dimension)
* Independent variables are e.g. TIME, TIME3, X, XB, RMAJM
* Dependent variable are a function of one or two of these variables.
* Most variables are a function of two independent variables so the standard behaviour requires specifying the index of the other independent variable. For example, if you wanted to plot Q(TIME, XB) as a function of time at the radial location idx = 50: 'python plot.py (filename) TIME Q 50'. However if you wanted to plot Q as a function of radius for a given time index (= 98): 'python plot.py (filename) XB Q 98'

Radial Grids and Gradients
--------------------------
* TRANSP uses rho = psi_n= sqrt(psi_t / psi_LCFS) where psi_t is the toroidal flux and psi_LCFS is the toroidal flux of the last closed flux surface.
* Miller uses rho = a_n = diameter/diameter of LCFS where diameter refers to the diameter of the flux surface at the height of the magnetic axis.
* Can calculate a geometric coefficient which relates gradients on these two radial grids as follows:
  * Essentially want to calculate dpsi_n / da_n.
  * TRANSP file contains psi_t(rmaj) = TRFMP(RMAJM, TIME) which can be used to calculate psi_n at output_radius as well as either side
  * Can use variable RMAJM(RMAJM, TIME) to calculate the a_n at output_time.
  * Use finite difference equation to calculate required differential coefficient.
* Once coefficient has been obtained, can multiply all gradient variables by this number to find what the gradient on a_n grid.
* Simple variables such as temperature, density, etc. don't require any changes since they are the same in both grids however gradients are defined on the a_n grid when used with Miller. 

Collisions
----------

* Equations for collision times/frequencies taken from Wesson Section 2.15 (T(keV), n(m^-3), mi=ion mass, mp=proton mass):
  * nu_i = 1/tau_i = 1 / [6.60e17 ((mi/mp)^1/2 Ti^3/2) / (n Z^4 ln(Lambda))] s^-1
  * nu_e = 1/tau_e = 1 / [1.09e16 (Te^3/2) / (n Z^2 ln(Lambda))] s^-1
  * log(Lambda) = Coulomb logarithm
    *log_i = 17.3 - 0.5*np.log(n/1e20) + 1.5*np.log(Ti) 
    *log_e = 14.9 - 0.5*np.log(n/1e20) + np.log(Te) 
  * lambda = debye length = 2.35e5 (T/n) m
* For the purposes of collision frequency calculations, ni = ne = n to satisfy quasineutrality. 
* From above explanation it should be clear that collisions between ion species is not being calculated. Only the collision frequency of the ion species with itself is being calculated. To do a GS2 simulation with many species you will need to recalculate collision frequencies consistent with the number of species included in the simulation. The default for this code works for ion simulations with AE and kinetic ions and electrons.

