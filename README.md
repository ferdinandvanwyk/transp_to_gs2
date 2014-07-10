transp_read
===========

Reads a TRANSP file and outputs parameters relevant to a GS2 simulation

Input:
-----
* Takes the name of the TRANSP file as first command line input
* Second command line input is radial location of output (will interpolate)
* Third command line input is time for output (will not interpolate)

Output:
-------
* Produces a file called gs2.in as output with GS2 parameters at specified radial location
