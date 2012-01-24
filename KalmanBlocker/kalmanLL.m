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
% Computes log likelihood for state space model
% with or without control impulses and affine term in observation equation
%
function loglik = kalmanLL( Z, A, W, F, S, x0, V0, a, B, U)

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

% Initialize log likelihood
loglik = 0;

for i=1:T
	% Time update
	if i==1
		xhatm = F*x0 + B*U(i,:)';
		Vm = F*V0*F' + S;
	else
		xhatm = F*xhat + B*U(i,:)';
		Vm = F*V*F' + S;
	end
	
	% Correction
	H = A*Vm*A' + W;
	K = Vm * A' * H^-1;
	xhat = xhatm + K*( ( Z(i,:) )' - a - A*xhatm );
	V = ( eye(k) - K*A) * Vm;
	
	% Update log likelihood
	loglik = loglik - (n/2)*log(2*pi) - .5 * log( abs(det(H)) ) -.5*( Z(i,:)' - a - A*xhatm )'*inv(H)*( Z(i,:)' - a - A*xhatm );
end