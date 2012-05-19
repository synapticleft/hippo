KMBOX - Kernel Methods Toolbox
a MATLAB toolbox for nonlinear signal processing and machine learning
Version 0.6, March 26th 2012
Author: Steven Van Vaerenbergh

http://sourceforge.net/p/kmbox


ABOUT
=====

The Kernel Methods Toolbox (KMBOX) is a collection of MATLAB programs that 
implement kernel-based algorithms, with a focus on adaptive filtering 
algorithms. It can be used for nonlinear signal processing and machine 
learning.

KMBOX includes implementations of algorithms such as kernel principal 
component analysis (KPCA), kernel canonical correlation analysis (KCCA) and 
kernel recursive least-squares (KRLS).

The goal of this distribution is to provide easy-to-analyze algorithm 
implementations, which reveal the inner mechanics of each algorithm and 
allow for quick modifications. The focus of these implementations is 
therefore on readability rather than speed or memory usage.

The basis of this toolbox was a set of programs written for the Ph.D. Thesis 
"Kernel Methods for Nonlinear Identification, Equalization and Separation of 
Signals", by Steven Van Vaerenbergh in 2009.

Template files are provided to encourage external authors to include their 
own code into the toolbox.

COPYRIGHT NOTICE
================

This program is free software: you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free 
Software Foundation, version 3 (as included and available at 
http://www.gnu.org/licenses).


INSTALLATION
============

Extract the zip file. Add the toolbox' main folder and subfolders to the 
MATLAB path.


USAGE
=====

The name of each function uses the prefix "km_" to minimize interference 
with other toolboxes. Usage of each function is specified in the function 
file itself.

Most algorithms have a corresponding demonstration file in the "demo" folder 
that bears the suffix "_demo". These programs can be executed without 
setting any additional parameters.

The code uses the following conventions:
- For data matrices, data is stored and accessed in row format: each data 
  point is a row in the data matrix.


EXTENDING KMBOX
===============

The files "km_template.m" in the base directory and "km_template_demo" in 
the "/demo" directory can be used as basic templates to extend KMBOX.

If you want to have your programs included in the KMBOX distribution on 
sourceforge, please email them to steven (at) gtas.dicom.unican.es. Do 
include at least one "demo" file for each algorithm.


CITING KMBOX
============

Published reports of research using this code (or a modified version) should 
cite the following document:
  Steven Van Vaerenbergh, Kernel Methods Toolbox (KMBOX): a MATLAB toolbox 
  for nonlinear signal processing and machine learning, 2010. Software 
  available at http://sourceforge.net/p/kmbox
The bibtex format is 
@misc{vanvaerenbergh_kmbox,
  author = {Steven Van Vaerenbergh},
  title = {Kernel Methods Toolbox {KMBOX}: a {MATLAB} toolbox for nonlinear 
  signal processing and machine learning},
  year = {2011},
  howpublished = {Grupo de Tratamiento Avanzado de Se\~nal, Departamento de 
  Ingenier\'ia de Comunicaciones, Universidad de Cantabria, Spain},
  note = {Software available at \url{http://sourceforge.net/p/kmbox}}
}


CONTACT
=======

Mail: steven (at) gtas.dicom.unican.es
Web site: http://www.gtas.dicom.unican.es/members/steven
Project web site: http://sourceforge.net/p/kmbox


INCLUDED ALGORITHMS
===================

1.  Kernel Ridge Regression (KRR).
2.  Principal Component Analysis (PCA).
3.  Kernel Principal Component Analysis (KPCA), as proposed in B. Scholkopf, 
    A. Smola and K.R. Muller, "Nonlinear component analysis as a kernel 
    eigenvalue problem", Neural Computation, volume 10, no. 5, pages 1299-
    1319, 1998.
4.  Approximate Linear Dependency Kernel Recursive Least-Squares (ALD-KRLS), 
    as proposed in Y. Engel, S. Mannor, and R. Meir. "The kernel recursive 
    least- squares algorithm", IEEE Transactions on Signal Processing, 
    volume 52, no. 8, pages 2275–2285, 2004.
5.  Sliding-Window Kernel Recursive Least-Squares (SW-KRLS), as proposed in 
    S. Van Vaerenbergh, J. Via, and I. Santamaria. "A sliding-window kernel 
    RLS algorithm and its application to nonlinear channel identification", 
    2006 IEEE International Conference on Acoustics, Speech, and Signal 
    Processing (ICASSP), Toulouse, France, 2006.
6.  Naive Online Regularized Risk Minimization Algorithm (NORMA), as 
    proposed in J. Kivinen, A. Smola and C. Williamson. "Online Learning 
    with Kernels", IEEE Transactions on Signal Processing, volume 52, no. 8, 
    pages 2165-2176, 2004.
7.  Fixed-Budget Kernel Recursive Least-Squares (FB-KRLS), as proposed in S. 
    Van Vaerenbergh, I. Santamaria, W. Liu and J. C. Principe, "Fixed-Budget 
    Kernel Recursive Least-Squares", 2010 IEEE International Conference on 
    Acoustics, Speech, and Signal Processing (ICASSP 2010), Dallas, Texas, 
    U.S.A., March 2010.
8.  Incomplete Cholesky Decomposition (ICD), as proposed in Francis R. Bach 
    and Michael I. Jordan. "Kernel Independent Component Analysis", Journal 
    of Machine Learning Research, volume 3, pages 1-48, 2002.
9.  Kernel Recursive Least-Squares Tracker (KRLS-T), as proposed in M. 
    Lazaro-Gredilla, S. Van Vaerenbergh and I. Santamaria, "A Bayesian 
    Approach to Tracking with Kernel Recursive Least-Squares", 2011 IEEE 
    International  Workshop on Machine Learning for Signal Processing (MLSP 
    2011), Beijing, China, September, 2011.
10. Kernel Canonical Correlation Analysis (KCCA), as proposed in  D. R. 
    Hardoon, S. Szedmak and J. Shawe-Taylor, "Canonical Correlation 
    Analysis: An Overview with Application to Learning Methods", Neural 
    Computation, Volume 16 (12), Pages 2639--2664, 2004.


HISTORY
=======

History of changes:
v0.6 (20120326):
- inclusion of a demo for kernel canonical correlation analysis
v0.5 (20120214):
- inclusion of KRLS-T
- addition of a file identifier to each file
v0.4 (20110504):
- inclusion of NORMA, fixed-budget KRLS, kernel PCA, incomplete Cholesky 
  decomposition
- inclusion of incomplete cholesky decomposition algorithm (km_kernel_icd).
- included a listing of dependencies in function headers.
- format change: dafault format for data matrices is now  one data point per 
  row (instead of one per column).
- format change: one input argument less for online algorithms
v0.3 (20101203):
- modifications of ALD-KRLS implementation.
v0.2 (20101108):
- inclusion of kernel recursive least-squares algorithms (km_krls): ALD-KRLS 
  (Approximate Linear Dependency KRLS), SW-KRLS (Sliding-Window KRLS).
- correction of minor details
v0.1 (20100908):
- original package, includes linear PCA and kernel ridge regression 
  algorithms.