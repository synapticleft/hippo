function [newVectors, whiteningMatrix, dewhiteningMatrix, zerophaseMatrix] = gwhitenv ...
    (vectors, E, D, cumVar, s_verbose);
%WHITENV - Whitenv vectors.
%
% [newVectors, whiteningMatrix, dewhiteningMatrix] = ...
%                               whitenv(vectors, E, D, verbose);
%
% Whitens the data (row vectors) and reduces dimension. Returns
% the whitened vectors (row vectors), whitening and dewhitening matrices.
%
% ARGUMENTS
%
% vectors       Data in row vectors.
% E             Eigenvector matrix from function 'pcamat'
% D             Diagonal eigenvalue matrix from function 'pcamat'
% verbose       Optional. Default is 'on'
%
% EXAMPLE
%       [E, D] = pcamat(vectors);
%       [nv, wm, dwm] = whitenv(vectors, E, D);
%
%
% This function is needed by FASTICA and FASTICAG
%
%   See also PCAMAT

% @(#)$Id: whitenv.m,v 1.3 2003/10/12 09:04:43 jarmo Exp $

% ========================================================
% Default value for 'verbose'
if nargin < 5, s_verbose = 'on'; end

% Check the optional parameter verbose;
switch lower(s_verbose)
 case 'on'
  b_verbose = 1;
 case 'off'
  b_verbose = 0;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''verbose''\n', s_verbose));
end

% ========================================================
% In some cases, rounding errors in Matlab cause negative
% eigenvalues (elements in the diagonal of D). Since it
% is difficult to know when this happens, it is difficult
% to correct it automatically. Therefore an error is 
% signalled and the correction is left to the user.
if any (diag (D) < 0),
  error (sprintf (['[ %d ] negative eigenvalues computed from the' ...
		   ' covariance matrix.\nThese are due to rounding' ...
		   ' errors in Matlab (the correct eigenvalues are\n' ...
		   'probably very small).\nTo correct the situation,' ...
		   ' please reduce the number of dimensions in the' ...
		   ' data\nby using the ''lastEig'' argument in' ...
		   ' function FASTICA, or ''Reduce dim.'' button\nin' ...
		   ' the graphical user interface.'], ...
		  sum (diag (D) < 0)));
end

% ========================================================
% Calculate the whitening and dewhitening matrices (these handle
% dimensionality simultaneously).
factors = diag(D);
noise_factors = ones(size(D,1),1);
rolloff_ind = sum(cumsum(flipud(factors))/cumVar > .995)
%rolloff_ind = sum(factors < .002)
%rolloff_ind = floor(size(D,1)/20);%rdim - sum(d(1:rdim)>variance_cutoff*p.whitening.X_noise_fraction);
noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi-.01,0,rolloff_ind))); 
%factors = factors./noise_factors;
%plot(noise_factors);pause(1);
% rolloff_ind = ceil(size(D,1)*4/5);%rdim - sum(d(1:rdim)>variance_cutoff*p.whitening.X_noise_fraction);
% noise_factors(rolloff_ind+1:end) = .5*(1+cos(linspace(0,pi-.01,size(D,1)-rolloff_ind))); 
% factors = factors.*noise_factors;
D = diag(factors./noise_factors);
whiteningMatrix = (sqrt (D)) \ E';
dewhiteningMatrix = E * sqrt (D);
rolloff_ind = max(2,sum(cumsum(flipud(factors))/cumVar < .9))
noise_factors(:) = 1;
noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi-.01,0,rolloff_ind))); 
D = diag(flipud(1./noise_factors));
zerophaseMatrix = E*inv (sqrt (D))*E';
% Project to the eigenvectors of the covariance matrix.
% Whiten the samples and reduce dimension simultaneously.
if b_verbose, fprintf ('Whitening...\n'); end
newVectors =  whiteningMatrix * vectors;
% ========================================================
% Just some security...
if ~isreal(newVectors)
  error ('Whitened vectors have imaginary values.');
end

% Print some information to user
if b_verbose
  fprintf ('Check: covariance differs from identity by [ %g ].\n', ...
    max (max (abs (cov (newVectors', 1) - eye (size (newVectors, 1))))));
end
