% KM_KCCA_DEMO Demo file for kernel canonical correlation analysis
% algorithm.
%
% This script takes two multi-dimensional variables and uses kernel CCA to
% map them onto a single latent 1-D variable. The kernel matrices are 
% decomposed using incomplete Cholesky decomposition in order to allow 
% large data sets. This demo includes 3 flavors of the KCCA generalized
% eigenvalue problem, all yielding very similar results.
%
% This file is part of the Kernel Methods Toolbox (KMBOX) for MATLAB.
% http://sourceforge.net/p/kmbox
% Author: Steven Van Vaerenbergh (steven *at* gtas.dicom.unican.es), 2012.
% Id: km_kcca_demo v1.0 2012/03/26
%
% The algorithm in this file is based on the following publications: 
% D. R. Hardoon, S. Szedmak and J. Shawe-Taylor, "Canonical Correlation 
% Analysis: An Overview with Application to Learning Methods", Neural 
% Computation, Volume 16 (12), Pages 2639--2664, 2004.
% F. R. Bach, M. I. Jordan, "Kernel Independent Component Analysis", Journal 
% of Machine Learning Research, 3, 1-48, 2002.
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation, version 3 (http://www.gnu.org/licenses).

close all; clear all
rs = 1; % random seed
randn('state',rs); rand('state',rs); %#ok<RAND>

%% PARAMETERS
N = 1000;	% number of samples
reg = 1E-5; % regularization
kerneltype = 'gauss';   % kernel type
kernelpar = 1;  % kernel parameter

%% PROGRAM
tic

% generate data
s = randn(N,1);	% latent signal
r1 = randn(N,1); r2 = randn(N,1);	% random (helper) variables

% option 1: two multi-dimensional variables that are mappable onto s
x1 = [tanh(r1-s)+0.1*r1 r1+3*s-1/10*(sin(3*s))];
x2 = [s - 2*(1-exp(-r2))./(1+exp(-r2)) r2.*s tanh(r2+s)];

% option 2: two invertible 1D nonlinear transformations
% x1 = tanh(0.8*s)+0.1*s;	% moderate saturation
% x2 = -1/10*(sin(s*3)+1.1*s*3);   % stairway

% clean up: remove mean
x1 = x1-repmat(mean(x1),N,1);
x2 = x2-repmat(mean(x2),N,1);

% normalize variance (either do this or use an anisotropic kernel)
x1 = x1*sqrt(diag(1./diag(x1'*x1)));
x2 = x2*sqrt(diag(1./diag(x2'*x2)));

% begin KCCA

% get incompletely decomposed kernel matrices. K1 \approx G1*G1'
G1 = km_kernel_icd(x1,kerneltype,kernelpar,50);
G2 = km_kernel_icd(x2,kerneltype,kernelpar,50);

% remove mean. avoid standard calculation N0 = eye(N)-1/N*ones(N);
G1 = G1-repmat(mean(G1),N,1);
G2 = G2-repmat(mean(G2),N,1);

% ones and zeros
N1 = size(G1,2); N2 = size(G2,2);
Z11 = zeros(N1); Z22 = zeros(N2); Z12 = zeros(N1,N2);
I11 = eye(N1); I22 = eye(N2);

% 3 GEV options, all of them are fairly equivalent

% % option 1: standard Hardoon
% R = [Z11 G1'*G1*G1'*G2; G2'*G2*G2'*G1 Z22];
% D = [G1'*G1*G1'*G1+reg*I11 Z12; Z12' G2'*G2*G2'*G2+reg*I22];

% option 2: simplified Hardoon
R = [Z11 G1'*G2; G2'*G1 Z22];
D = [G1'*G1+reg*I11 Z12; Z12' G2'*G2+reg*I22];

% % option 3: Kettenring-like generalizable formulation
% R = 1/2*[G1'*G1 G1'*G2; G2'*G1 G2'*G2];
% D = [G1'*G1+reg*I11 Z12; Z12' G2'*G2+reg*I22];

% solve generalized eigenvalue problem
[alphas,betas] = eig(R,D);
[betass,ind] = sort(real(diag(betas)));
alpha = alphas(:,ind(end)); alpha = alpha/norm(alpha);
beta = betass(end);

% expansion coefficients
alpha1 = alpha(1:N1);
alpha2 = alpha(N1+1:end);

% estimates of latent variable
y1 = G1*alpha1;
y2 = G2*alpha2;

% scale the estimated signals to compare without the scalar ambiguity
scaling = sqrt(var(s))/sqrt(var(y1))*sign(s(1))*sign(y1(1));

% mean square errors
error1 = s-scaling*y1;
error2 = s-scaling*y2;
MSE1 = sumsqr(error1)/N;
MSE2 = sumsqr(error2)/N;

toc
%% OUTPUT

figure; hold all
plot(s)
plot(scaling*y1);
plot(scaling*y2);
legend('latent variable','projection 1','projection 2')

fprintf('Canonical correlation: %f\n',beta)
fprintf('MSE1: %f\n',MSE1);
fprintf('MSE2: %f\n',MSE2);
fprintf('\n')

% figure;plot(sort(diag(real(betas))))  % check eigenvalues
% figure;plot(sort(diag(betas))); % imaginary part due to numerical error