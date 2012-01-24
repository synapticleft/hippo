% Copyright 2007 Alexander W Blocker
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU Lesser General Public License for more details.
% You should have received a copy of the GNU Lesser General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%
% Author: Alexander W Blocker
% Uses EM algorithm to estimate parameters of linear Gaussian state-space system
% with or without control inputs and an affine term in the measurement equation.
% Also allows for diagonal restrictions on system & measurement covariance matrices.
%
function [A, W, F, S, a, B, x0, V0, loglik] = kalmanMLE( Z, A0, W0, F0, S0, x0, V0, diagS, diagW, a0, B0, U)

converged = 0;
maxIter = 1000;
thresh = 1e-5;
iter = 0;

A=A0;
W=W0;
F=F0;
S=S0;

% parse input
T = size(Z,1);
n = size(Z,2);
k = size(F,1);

if nargin < 8, diagS = 0; end;
if nargin < 9, diagW = 0; end;

if nargin < 10
    a=0;
else
    a=a0;
end

if nargin < 11
    B=0;
else
    B=B0;
end

if nargin >= 12
	if size(U, 1) ~= T
		error('Dimension error - U & Z must have same number of observations');
	end
end

% Set type of model
if (a==0 & B==0)
    type=0;
elseif B==0
    type=1;
elseif a==0
    type=2;
else
	type=3;
end

fprintf(1, 'type: %d\n', type);

% Setup data object for log-likelihood
loglik = [];

% EM loop
while (iter < maxIter & converged == 0)
    
    i = iter + 1;
    
    % Smoothing estimates, selected on model type
    switch type
        case 0
            [xsmooth, Vsmooth, VVsmooth, loglik(i)] = kalmanSmooth(Z, A, W, F, S, x0, V0);
			[A, W, F, S, x0, V0] = kalmanMStep0(Z, xsmooth, Vsmooth, VVsmooth, diagS, diagW);
        case 1
            [xsmooth, Vsmooth, VVsmooth, loglik(i)] = kalmanSmooth(Z, A, W, F, S, x0, V0, a);
			[A, W, F, S, x0, V0, a] = kalmanMStep1(Z, xsmooth, Vsmooth, VVsmooth, diagS, diagW);
        case 2
            [xsmooth, Vsmooth, VVsmooth, loglik(i)] = kalmanSmooth(Z, A, W, F, S, x0, V0, a, B, U);
			[A, W, F, S, x0, V0, B] = kalmanMStep2(Z, xsmooth, Vsmooth, VVsmooth, diagS, diagW, U);
		case 3
			[xsmooth, Vsmooth, VVsmooth, loglik(i)] = kalmanSmooth(Z, A, W, F, S, x0, V0, a, B, U);
			[A, W, F, S, x0, V0, a, B] = kalmanMStep3(Z, xsmooth, Vsmooth, VVsmooth, diagS, diagW, U);
    end
	
	% Increment counter
    iter = iter + 1;
	
	% Print status
	fprintf(1, 'iteration %d, loglik = %f\n', iter, loglik(i) );
	
    if iter > 1
        % Check for convergence - relative change criterion
        convtest = abs(loglik(iter) - loglik(iter-1)) / ...
            ( (abs(loglik(iter)) + abs(loglik(iter-1)) + eps)/2 );
        if ( (convtest < thresh) & (iter > ceil(maxIter*.01) ) )
            converged = 1;
        end
        
        % Check for reversal
		if ( (loglik(iter-1) - loglik(iter)) > 1e-3 )
			error('EM algorithm reversed - terminating.\n');
		end
    end
	
end