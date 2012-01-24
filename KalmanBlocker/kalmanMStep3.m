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
% M-step for EM algorithm to estimate parameters of linear Gaussian state space
% system with control inputs and an affine term in the measurement equation.
%
function [A, W, F, S, x0, V0, a, B] = kalmanMStep3(Z, xsmooth, Vsmooth, VVsmooth, diagS, diagW, U)

T = size(Z,1);

zt = 0;
xt = 0;
ztxt = 0;
ztzt = 0;
Pt1 = 0;
Pttm1 = 0;
xtut = 0;
xtm1ut = 0;
utut = 0;

for i=1:T
	zt = zt + Z(i,:)';
	xt = xt + xsmooth(i,:)';
	ztxt = ztxt + Z(i,:)' * xsmooth(i,:);
	ztzt = ztzt + Z(i,:)' * Z(i,:);
	Pt1 = Pt1 + Vsmooth(:,:,i) + xsmooth(i,:)' * xsmooth(i,:);
end

Pt2 = Pt1 - Vsmooth(:,:,1) - xsmooth(1,:)' * xsmooth(1,:);
Ptm1 = Pt1 - Vsmooth(:,:,T) - xsmooth(T,:)' * xsmooth(T,:);

for i=2:T
	Pttm1 = Pttm1 + VVsmooth(:,:,i) + xsmooth(i,:)' * xsmooth(i-1,:);
	xtut = xtut + xsmooth(i,:)' * U(i,:);
	xtm1ut = xtm1ut + xsmooth(i-1,:)' * U(i,:);
	utut = utut + U(i,:)' * U(i,:);
end

a = (1 / T) * (zt - (ztxt/Pt1)*xt) * (1 - (1/T)*(xt'/Pt1)*xt)^-1;
A = (ztxt - a*xt')/ Pt1;
W = (1/T) * (ztzt - A*ztxt' - a*zt');
if diagW==1, W = diag(diag(W)); end;
B = ( xtut - (Pttm1/Ptm1)*xtm1ut ) / (utut + (xtm1ut'/Ptm1)*xtm1ut);
F = (Pttm1 - B*xtm1ut') / Ptm1;
S = (1/(T-1)) * (Pt2 - F*Pttm1' - B*xtut');
if diagS==1, S = diag(diag(S)); end;
x0 = xsmooth(1,:)';
V0 = Vsmooth(:,:,1);