# CPDINV

Software for estimating Curie depth using centroid method with wavelet spectrum and Fourier spectrum The Curie depth is an important magnetic interface for studying the geological 
structure and thermal evolution of the crust and lithosphere, and it is common to estimate the Curie-point depth by spectral analysis of magnetic anomaly data. The software 
combines the Fourier transform and wavelet transform with the centroid method, respectively, and then complete the whole program of inversion of Curie-point depth based on FORTRAN 
language, which can run stably in Linux system.

This package consists of six applications, including projection, interpolation, power spectrum calculation (Fourier transform and wavelet transform), the least squares fitting of 
the centroid method, Curie-point depth calculation, which is installed and running on a Linux system with GFORTRAN 4.8 and above version.


We refer to the USGS potential field software (Cordell et al. 1993; Phillips, 1997), Kirby et al. (2005) and Wang et al. (2010) in writing our software.

Authors: Yihong Yin and Chunfeng Li
