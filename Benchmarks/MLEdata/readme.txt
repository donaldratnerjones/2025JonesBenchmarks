The files X.csv and y.csv are intended to be used for testing algorithms 
for maximizing the kriging concentrated log-likelihood function.

The data is based on the first constraint of the automotive problem, which 
depends upon 59 variables.  The matrix in X.csv is a space-filling experimental 
design with 590 points (10 times the number of variables).  The vector in
y.csv is the value of the first constraint of the automotive problem at 
these 590 points.

The files Xbranin.csv and ybranin.csv contain input and output data for a
49-point experimental design (7x7 grid) for the Branion test function.
The input data is normalized to the unit interval.  This data can be used
to study how the smoothness and multimodality of the concentrated 
log-likelihood function depend upon the correlation function used.