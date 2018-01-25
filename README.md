Shameless MATLAB Port of Krishnamoorthy's "Two Poisson Means" fortran code
==================================================================

Port by Andrew Leifer
leifer@princeton.edu
24 January 2018 

In their paper: Krishnamoorthy, K and Thomson, J. (2004)  A more powerful test for comparing two Poisson means. Journal of  Statistical Planning and Inference, 119, 249-267]

Krishnamoorthy and Thomson describe an improved test to determine if two Poissonion means are significantly different, now referred to as the "E-test." Krishnamoorthy released a fortran implementation here: 

Original Code http://www.ucs.louisiana.edu/~kxk4695/
and http://www.ucs.louisiana.edu/~kxk4695/statcalc/pois2pval.for

Here is a quick and dirty port of  their fortran code to MATLAB. Note I kept all of the original fortran structure and variable names and style. 
