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
% Computes Kalman smoothed state estimates for state space model
% with or without control inputs and an affine term in observation equation
%
function [xsmooth, Vsmooth, VVsmooth, loglik] = kalmanSmooth( Z, A, W, F, S, x0, V0, a, B, U)

% parse input
T = size(Z,1);
n = size(Z,2);
k = size(F,1);

if n~=size(A,1)
	error('Dimension error - A & Z must agree.');
end

if (n~=size(W,1) || n~=size(W,2))
	error('Dimension error - W is not valid covariance matrix.');
end

if (k~=size(F,2))
	error('Dimension error - F not square.');
end

if (k~=size(A,2))
	error('Dimension error - F & A do not agree.');
end

if (k~=size(S,1) || k~=size(S,2))
	error('Dimension error - S is not valid.');
end

if nargin < 8, a=0; end
if nargin < 9, B=0; end

if nargin < 10
	U = zeros( size(Z,1), 1 );
else
	if size(U, 1) ~= T
		error('Dimension error - U & Z must have same number of observations');
	end
end

% run Kalman filter
[xfilt, Vfilt, Vmfilt, loglik] = kalmanFilter( Z, A, W, F, S, x0, V0, a, B, U);

% setup data objects for smoothed observations & covariance matrices
xsmooth = zeros( T, k);
Vsmooth = zeros( k, k, T);
VVsmooth = zeros( k, k, T);

xsmooth(T,:) = xfilt(T,:);
Vsmooth(:,:,T) = Vfilt(:,:,T);
K = Vmfilt(:,:,T) * A' / (A*Vmfilt(:,:,T)*A' + W);
VVsmooth(:,:,T) = ( eye(k) - K * A) * F * Vfilt(:,:,T-1);

% run backwards recursions
J = zeros(k,k,T);
for i=1:(T-1)
	j = T - i;
	Jt = Vfilt(:,:,j) * F' / Vmfilt(:,:,j+1);
	xsmooth( j, :) = ( xfilt(j,:)' + Jt*( xsmooth(j+1,:)' - F*xfilt(j,:)' - B*U(j,:)' ) )';
	Vsmooth( :, :, j) = Vfilt(:,:,j) + Jt*( Vsmooth(:,:,j+1) - Vmfilt(:,:,j+1) )*Jt';
    J(:,:,j) = Jt;
end

for i = 1:(T-2)
    j = T - i;
    VVsmooth(:,:,j) = Vfilt(:,:,j) * J(:,:,j-1)' ...
		+ J(:,:,j) * ( VVsmooth(:,:,j+1) - F*Vfilt(:,:,j) ) * J(:,:,j-1)';
end

% Calculate log likelihood - alternate method

% loglik = 0;

% for i = 1:T
	% H = A*Vsmooth(:,:,i)*A' + W;
	% loglik = loglik - (n/2)*log(2*pi) - .5 * log( det(H) ) -.5*( Z(i,:)' - a - A*xsmooth(i,:)' )'*(H^-1)*( Z(i,:)' - a - A*xsmooth(i,:)' );
% end