function [A,W,t] = icaDemo(v,Xf,pos)
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

%%Processing of position information
pos(pos == -1) = nan;
reject = 0;
for i = 1:4
    reject = reject | min([0; diff(pos(:,i))],flipud([0; diff(flipud(pos(:,i)))])) < -20;
end
pos(reject,:) = nan;
if size(Xf,2) < size(pos,1)
    pos = pos(1:size(Xf,2),:);
end
for i = 1:4
    nanInds = find(~isnan(pos(:,i)));
    pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
end
nanInds = isnan(pos(:,1)) | isnan(pos(:,3));
pos = pos(~nanInds,:);Xf = Xf(:,~nanInds); v= v(~nanInds,:);
vel = angVel(pos);
vel = [0; vel(:,1)];
pos = bsxfun(@minus,pos,mean(pos));
[a,~,~] = svd(pos(:,1:2),'econ');pos = a;
for i = 1:2    
    pos(:,i) = pos(:,i) - min(pos(:,i));
    pos(:,i) = pos(:,i)/(max(pos(:,i)));
    pos(:,i) = min(pos(:,i),.9999);
end
%%THE filtLow function requires some functions from the signal processing
%%toolbox, but is not particularly necessary.
vel = filtLow(vel,1250/32,1);
vel = vel/max(vel);
inds = vel > thresh;
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
Xf = Xf(:,inds);pos = pos(inds,:);v = v(inds);
Xf = [real(Xf);imag(Xf)];
rdim = size(Xf,1);
[A,W] = gfastica(zscore(Xf,0,2),'lastEig',rdim,'g','tanh','approach','symm','stabilization','on');
t = W*zscore(Xf,0,2);

