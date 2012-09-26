function [A,W,t] = icaDemo(v,Xf)
%% convert demodulated complex data to real valued w/ 2x dimensionality, 
%% run fastICA, then bin and render activations.
% Inputs:
% Xf = filtered data (for example 64x100000 complex valued matrix is 64 electrodes at 100000 sampled timepoints.
% v = time course of 1st PC (determined by [~,~,v] = svds(Xf,1); -- svds is a matlab builtin function).
% Outputs:
% A = ICA mixing matrix, contains projections of each component onto electrode space.
% W = ICA unmixing matrix, extracts component activations from original electrode space.
% t = IC activations
% CAVEAT: electrode space is 2*(x+1) dimensions, where x is the number of
% electrodes. +1 because, in addition to electrode values, I feed in a value 
% that captures the evolution of v; 2* because of real and imaginary components fed into ICA.
offSet = 1;
if 1
    % Method 1 - demodulation - multiplies Xf by phase of v
    Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).');...
        [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))].'];
else
    % Method 2 - cross product - multiplies Xf by v -- works similarly, can be
    % more stable when v is of small amplitude and its phase is poorly defined.
    % This method can be used instead if the math is easier for you.
    Xf = [bsxfun(@times,Xf,v(:,1).');[zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))].'];
end
Xf = [real(Xf);imag(Xf)];
rdim = size(Xf,1);
[A,W] = gfastica(zscore(Xf,0,2),'lastEig',rdim,'g','tanh','approach','symm','stabilization','on');
t = W*zscore(Xf,0,2);
