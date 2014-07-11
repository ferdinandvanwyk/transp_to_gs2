transp_read
===========

Reads a TRANSP file and outputs parameters relevant to a GS2 simulation

Input:
-----
* Takes the name of the TRANSP file as first command line input
* Second command line input is radial location of output (will interpolate)
* Third command line input is time for output (will not interpolate)
* Input command therefore looks like: python main.py <filename> <radius> <time>

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
* Use NetCDF since it seems to be the only thing that reads the TRANSP file. Could try SciPy NetCDF functions as documented here: http://docs.scipy.org/doc/scipy/reference/generated/scipy.io.netcdf.netcdf_file.html 
