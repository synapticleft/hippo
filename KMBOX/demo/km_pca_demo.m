% KM_PCA_DEMO Principal component analysis on a simple two-dimensional
% data set. The first principal component is retrieved and plotted,
% together with both principal directions.
%
% Author: Steven Van Vaerenbergh (steven *at* gtas.dicom.unican.es), 2010.
% Id: km_pca_demo.m v1.0
% This file is part of the Kernel Methods Toolbox (KMBOX) for MATLAB.
% http://sourceforge.net/p/kmbox
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation, version 3 (as included and available at
% http://www.gnu.org/licenses).

close all
clear all

%% PARAMETERS

N = 500;			% number of data points
R = [-.9 .4; .1 .2];	% covariance matrix

%% PROGRAM
tic

X = randn(N,2)*R;	% correlated two-dimensional data

[E,v,Xp] = km_pca(X,1);		% obtain eigenvector matrix E, eigenvalues v and principal components Xp

toc
%% OUTPUT

figure; hold on
plot(X(:,1),X(:,2),'.')
plot(E(1,1)*Xp,E(2,1)*Xp,'.r')
plot([0 E(1,1)],[0 E(2,1)],'g','LineWidth',4)
plot([0 E(1,2)],[0 E(2,2)],'k','LineWidth',4)
axis equal
legend('data','first principal components','first principal direction','second principal direction')
title('linear PCA demo')

