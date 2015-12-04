transp_to_gs2
=============

Reads a TRANSP output file and outputs parameters relevant to a GS2 simulation

Input:
------

* Function takes in three command line parameters with an optional fourth:
   * Path to TRANSP file
   * Radial location (will interpolate): value of the square root of the 
     normalized toroidal flux in range [0,1].
   * Time (will determine nearest time index, i.e. won't interpolate).
   * (optional) Template GS2 input file (relevant parameters will be replaced)
   * Input command therefore looks like:

.. code:: bash
        
   $ python main.py </path/to/transp/file> <radius> <time> (<gs2 input file>)
     
Caveats
-------

* **The sign of the flow shear may not be correct!** This program simply outputs 
  the gradient of the flow without changing sign. The real sign of the flow 
  shear will depend on the sign convention used in the experiment.

Output:
-------
* Produces a file called ``gs2_rho_<radius>_time_<time>.in``
* Equilibrium parameters, e.g. temperatures, densities, gradients (see note
  below), flows, shears
* Geometric Miller parameters (with appropriate corrections for GS2, see
  http://arxiv.org/abs/1403.3293)
* Miscellaneous parameters which need to be set when using Miller parameters
  (only option for this code so far).
* Various useful reference parameters.

Requirements:
-------------

* Python (at least 2.7)
* NetCDF4
* HDF5
* Python NetCDF4 Interface: https://pypi.python.org/pypi/netCDF4
* Python packages are listed in the ``requirements.txt`` file

Radial Grids and Gradients
--------------------------

* TRANSP uses rho = psi_tor = sqrt(psi_tor / psi_tor,LCFS) where psi_tor is the
  toroidal flux and psi_tor,LCFS is the toroidal flux of the last closed flux 
  surface.
* Miller uses rho = a_n = diameter/diameter of LCFS where diameter refers to 
  the diameter of the flux surface at the height of the magnetic axis.
* Can calculate a geometric coefficient which relates gradients on these two radial grids as follows:
   * Essentially want to calculate dpsi_n / da_n.
   * TRANSP file contains psi_t(rmaj) = TRFMP(RMAJM, TIME) which can be used to
     calculate psi_n at output_radius as well as either side
   * Can use variable RMAJM(RMAJM, TIME) to calculate the a_n at output_time.
   * Use finite difference equation to calculate required differential coefficient.
* Once coefficient has been obtained, can multiply all gradient variables by 
  this number to find what the gradient on a_n grid.
* Simple variables such as temperature, density, etc. don't require any changes
  since they are the same in both grids however gradients are defined on the 
  a_n grid when used with Miller.

Collisions
----------

* The collision frequencies are calculated using the equation from the GS2 wiki
  (http://gyrokinetics.sourceforge.net/wiki/index.php/Gs2_Input_Parameters).
* More convenient equations (calculated by M. Barnes, where densities are in 1e19, temperatures are in keV, mi is the ion mass in multiples of the proton mass, zi is the ion charge which depends on the ionization):
   * log(lambda) = 24 - log(1e4 * sqrt(n_ref/10) / t_e)
   * vnewk_i = 9.21e-5 * amin * zi^4 * loglam * n_ref / (sqrt(2) * t_ref^2)
   * vnewk_e = 3.95e-3 * amin * sqrt(0.5 * mi) * loglam * n_ref / (sqrt(t_ref) *
     t_e^1.5)
* For the purposes of collision frequency calculations, ni = ne = n to satisfy
  quasineutrality.

Plotting Routine
----------------
A rudimentarty plotting routine is also included and works in the following way:

* It can be called as: 

.. code:: bash

   python plot.py </path/to/transp/file> <independent variable> <dependent variable> <index of other dimension>

* Independent variables are e.g. TIME, TIME3, X, XB, RMAJM
* Dependent variable are a function of one or two of these variables.
* Most variables are a function of two independent variables so the standard 
  behaviour requires specifying the index of the other independent variable. 
  For example, if you wanted to plot Q(TIME, XB) as a function of time at the 
  radial location idx = 50: ``python plot.py (filename) TIME Q 50``. However if 
  you wanted to plot Q as a function of radius for a given time index (= 98): 
  ``python plot.py (filename) XB Q 98``
