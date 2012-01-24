The following notation is used throughout these packages:

x(t+1) = F*x(t) + B*u(t) + e(t), e(t) ~ N(0,S)
z(t) = a + A*x(t) + w(t), w(t)~ N(0,W)

x(t) is the k*1 state vector at time t
z(t) is the n*1 observation vector at time t
u(t) is the m*1 vector of known inputs at time t

F is the k*k state evolution matrix
B is the k*m input response matrix
S is the k*k covariance matrix for the Gaussian state innovations

A is the n*n observation matrix
a is an optional n*1 affine term in the measurement equation
W is the n*n covariance matrix for the Gaussian observation noise

Z is the T*n matrix of observations
U is the T*m matrix of known inputs

Vfilt & Vmfilt are n*n*T arrays containing the filtered prior and posterior
error covariance matrices for each time.

Vsmooth and VVsmmoth are n*n*T arrays containing smoothed error covariance
matrices for each time. Vsmooth is the estimated covariance between
contemporaneous state estimation errors. VVsmooth is the estimated covariance
between state estimation errors at times t & t-1.

All inputs suffixed with 0 (A0, W0, F0, S0, etc.) are the starting values
for the EM algorithm. These should be chosen quite carefully, as the EM algorithm
is only capable of finding local maximum, not the global one.

For more information on the EM algorithm & its application to linear state-space
systems, I recommend:

Ghahramani and Hinton, "Parameter Estimation for LDS", U. Toronto tech. report, 1996
Digalakis, Rohlicek and Ostendorf, "ML Estimation of a stochastic linear system with the EM
	algorithm and its application to speech recognition",
	IEEE Trans. Speech and Audio Proc., 1(4):431--442, 1993.
Borman, Sean, "The Expectation Maximization Algorithm, A Short Tutorial",
	Tutorial Notes, University of Notre Dame EE Department,
	http://www.seanborman.com/publications/EM_algorithm.pdf, July 2004.