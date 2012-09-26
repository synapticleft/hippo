%function [phi,whitenMatrix,dewhitenMatrix,zerophaseMatrix] = run_lfp_real(pos,v,Xf,thresh)

%% driver script for learning a model
%%%   
%%%   m - struct containing the basis function variables
%%%     m.patch_sz - image domain patch size (side length of square patch)
%%%     m.M - whitened domain dimensions
%%%     m.N - firstlayer basis function dimensions
%%%     m.L - phase transformation basis functions dimensions
%%%     m.K - amplitude basis function dimensions
%%%     m.t - number of learning iterations
%%%     
%%%     m.A - first layer complex basis functions (m.M x m.N)
%%%     m.D - second layer transformation components (m.N x m.L)
%%%     m.B - second layer amplitude components (m.N x m.K)
%%%     
%%%   p - struct containing the learning, inference, and other parameters
%%%     p.firstlayer - first layer complex basis function parameters
%%%     p.ampmodel - second layer amplitude component parameters
%%%     p.phasetrans - second layer phase transformation parameters
%%%     
%%%   Summary of learning proceedure:
%%%       1. Initialize parameters
%%%       2. Estimate whitening transform
%%%       3. Learn first layer complex basis functions
%%%       4. Infer large batch of first layer coefficients
%%%       5. Learn second layer phase tranformation components
%%%       6. Learn second layer amplitude components
%%%       7. Display the results

pos(pos == -1) = nan;
reject = 0;
for i = 1:4
    reject = reject | min([0; diff(pos(:,i))],flipud([0; diff(flipud(pos(:,i)))])) < -20;
end
pos(reject,:) = nan;
if size(v,1) < size(pos,1)
    pos = pos(1:size(v,1),:);
end
for i = 1:4
    nanInds = find(~isnan(pos(:,i)));
    pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
end
nanInds = isnan(pos(:,1)) | isnan(pos(:,3));
pos = pos(~nanInds,:);v = v(~nanInds,:);Xf = Xf(:,~nanInds);
vel = angVel(pos);
vel = [0; vel(:,1)];
pos = bsxfun(@minus,pos,mean(pos));
%[a,~,~] = svd(pos(:,1:2),'econ');pos = a;
for i = 1:2    
    pos(:,i) = pos(:,i) - min(pos(:,i));
    pos(:,i) = pos(:,i)/(max(pos(:,i)));
    pos(:,i) = min(pos(:,i),.9999);
end
vel = filtLow(vel,1250/32,1);
vel = vel/max(vel);inds = vel > thresh;
offSet = 1;
 Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1).')));...
      [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))].'];
%   Xf = [bsxfun(@times,Xf,v(:,1).');...
%    [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))].'];
Xf = Xf(:,inds);
Xf = [real(Xf);imag(Xf)];
Xf = bsxfun(@minus,Xf,mean(Xf,2));
Xf = zscore(Xf,0,2);

%% Initialize parameters

%% Whitening transform
verbose = 'on';
[E, D,cumVar]= gpcamat(Xf, 1, size(Xf,1),'off', verbose);
% d = diag(D);
% d = sort(d,'descend');
% figure;plot(log10(1-cumsum(d)/sum(d)));hold all;plot(log10(d));
% return
[Xf, whitenMatrix, dewhitenMatrix, zerophaseMatrix] = gwhitenv(Xf, E, D,cumVar, verbose);
unittest
% % nic = 60;
% % a = zeros(1,nic)';
% % b = ones(1,nic)';
% % neig = size(Xf,1);
% % W0 = randn(nic,neig)/10;
% % phi = hminu(Xf', nic, a, b, '', W0);%
% nHidden = 60; nInput = size(Xf, 1);
% phi = randn(nHidden, nInput); 
% options.Method  = 'lbfgs';
% options.maxIter = 100;	    % Maximum number of iterations of L-BFGS to run 
% options.display = 'on';
% [phi, cost] = minFunc( @icaScoreMatching, ...
%                         phi(:), options, Xf, ...
%                         nHidden, nInput);
% phi = reshape(phi, nHidden, nInput);
