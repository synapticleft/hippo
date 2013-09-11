function [A W Z] = cfastica(mixedsig)

% Complex FastICA
% Ella Bingham 1999
% Neural Networks Research Centre, Helsinki University of Technology

% This is simple Matlab code for computing FastICA on complex valued signals.
% The algorithm is reported in:
% Ella Bingham and Aapo Hyvärinen, "A fast fixed-point algorithm for 
% independent component analysis of complex valued signals", International 
% Journal of Neural Systems, Vol. 10, No. 1 (February, 2000) 1-8.

% When using the code, please refer to the abovementioned publication.

% Nonlinearity G(y) = log(eps+y) is used in this code; this corresponds to 
% G_2 in the above paper with eps=a_2.

% Some bugs corrected on Oct 2003, thanks to Ioannis Andrianakis for pointing them out

% Compatible with GPLv2 or later. http://www.gnu.org/licenses/gpl-2.0.txt

% Copyright (C) 1999  Ella Bingham
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps = 1; % epsilon in G
% Whitening of x:

% [Ex, Dx] = eig(cov(mixedsig'));
% d = flipud(diag(Dx));
% cumVar = sum(d);
% maxLastEig = sum(cumsum(d)/cumVar < .9999999)
% Dx = Dx(end-maxLastEig+1:end,end-maxLastEig+1:end);
% Ex = Ex(:,end-maxLastEig+1:end);
% factors = diag(Dx);
% noise_factors = ones(size(Dx,1),1);
% rolloff_ind = sum(cumsum(flipud(factors))/cumVar > .999999)
% noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi-.01,0,rolloff_ind))); 
% Dx = diag(factors./noise_factors);
% whiteningMatrix = sqrt(inv(Dx)) * Ex';%eye(size(mixedsig,1));%
% x = whiteningMatrix * mixedsig;
% dewhiteningMatrix = Ex * sqrt (Dx);%whiteningMatrix;%
% rolloff_ind = 2;
% noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi-.01,0,rolloff_ind))); 
% zerophaseMatrix = Ex*sqrt(diag(flipud(noise_factors)))*Ex';%inv (sqrt (D))*E';
% verbose = 'on';
% interactivePCA = 'off';
% [E, D, cumVar]= gpcamat(mixedsig, 1, size(mixedsig,1), interactivePCA, verbose);
% [x, whiteningMatrix, dewhiteningMatrix, zerophaseMatrix] = gwhitenv(mixedsig, E, D, cumVar, verbose);
%x = x - mean(x,2)*ones(1,m);

% Condition in Theorem 1 should be < 0 when maximising and > 0 when 
% minimising E{G(|w^Hx|^2)}. 
[x,whiteningMatrix,dewhiteningMatrix,zerophaseMatrix] = whiten(mixedsig);
%x = mixedsig;whiteningMatrix = eye(size(x,1));dewhiteningMatrix = eye(size(x,1));zerophaseMatrix = eye(size(x,1));
n = size(x,1);
%  C = cov(x');
  maxcounter = 100;
  counter = 0;
  W = randn(n,n) + 1i*randn(n,n);
  WOld = zeros(size(W));
  minAbsCos = 1;
  while counter < maxcounter && (abs(1-minAbsCos) > .0001 || counter == 0)
      WOld = W;
      Wx = W'*x;
      aWx2 = abs(Wx).^2;
      gWx = 1./(eps + aWx2);
      W = x*(Wx.*gWx)'/size(x,2) - bsxfun(@times,W,mean(gWx + aWx2.*(-gWx.^2),2).');
%       tic;
%     for j = 1:n
%       gWx(j,:) = 1./(eps + abs(W(:,j)'*x).^2);
%       dgWx(j,:) = -1./(eps + abs(W(:,j)'*x).^2).^2;
%       W(:,j) = mean(x .* (ones(n,1)*conj(W(:,j)'*x)) .* (ones(n,1)*gWx(j,:)),2) - mean(gWx(j,:) + abs(W(:,j)'*x).^2 .* dgWx(j,:)) * W(:,j);
%       %W(:,j) = mean(bsxfun(@times,x,conj(Wx(j,:)) .* gWx(j,:)),2) - mean(gWx(j,:) + aWx2(j,:) .* dgWx(j,:)) * W(:,j);
%     end
%     toc
%     figure;scatter(real(W1(:)),real(W(:)));hold all;scatter(imag(W1(:)),imag(W(:)),'r');drawnow;
    % Symmetric decorrelation:
    %W = W * sqrtm(inv(W'*W));
    [E,D] = eig(W'*W);%*C
    W = W * E * inv(sqrt(D)) * E';
    minAbsCos = min(abs(diag(W' * WOld)));
    meanAbsCos = mean(abs(diag(W' * WOld)));
    fprintf('Step no. %d, change in value of estimate: %.3g %.3g\n', counter, 1 - minAbsCos,1-meanAbsCos);
    counter = counter + 1;
    %imagesc(dewhiteningMatrix*W,[8 4]));drawnow;
  end
	A = dewhiteningMatrix * W;  
  	W = W' * whiteningMatrix;
Z = zerophaseMatrix*A;


% abs((Q*A)'*W) should be a permutation matrix; Figure 1 in the IJNS paper
% measures the error in this as 
% 	absQAHW = abs((Q*A)'*W);
%	maximum = max(absQAHW);
%	SE = sum(absQAHW.^2) - maximum.^2 + (ones(1,n)-maximum).^2;

function newMatrix = selcol(oldMatrix, maskVector)
% newMatrix = selcol(oldMatrix, maskVector);
%
% Selects the columns of the matrix that marked by one in the given vector.
% The maskVector is a column vector.
% 15.3.1998
if size(maskVector, 1) ~= size(oldMatrix, 2),
  error ('The mask vector and matrix are of uncompatible size.');
end

numTaken = 0;

for i = 1 : size (maskVector, 1),
  if maskVector(i, 1) == 1,
    takingMask(1, numTaken + 1) = i;
    numTaken = numTaken + 1;
  end
end

newMatrix = oldMatrix(:, takingMask);
