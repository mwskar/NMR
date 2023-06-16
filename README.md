# NMR Project

A program that takes in NMR data and calculates the number of Hydrogen Atoms.

Made as a project for Dr. Pounds' (Mercer University) Numerical Methods class.

- Order of processes:

- Read and sort data

- Adjust for baseline

- smooth data wiht filter option [Boxcar, Savitzky-Golay]

- Fit a cublic spline

- Find areas to integrate

- Integrate [Composite Newton Cotes, Romberg, Adaptive Quadurature, Gaussian]

- Find number of Hydrogen atoms



nmr.f95 reads 'nmr.in' by default


Program Input
--------------
Name of NMR file

Baseline adjustment

Tolerance for algorithms

Filter type (0=None, 1=Boxcar, 2=SG)

Size of Filter

Filter passes

Integration (0=Newton-Cotes,1=Romberg,2=Adaptive,3=Quadrature) 

output file
