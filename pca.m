function [signals,PC,V] = pca(data,version)
% data - MxN matrix of input data
% (M dimensions, N trials)
% signals - MxN matrix of projected data
% PC - each column is a PC
% V - Mx1 matrix of variances

%[M,N] = size(data);
mn = mean(data,2);% subtract off the mean for each dimension
data = bsxfun(@minus,data,mn);%data - repmat(mn,1,N);

if version == 1 %covariance
%    covariance = 1 / (N-1) * data * data';
    [PC, V] = eig(cov(data'));%ariance);
    V = diag(V);
    [junk, rindices] = sort(-1*V);
    V = V(rindices);
    PC = PC(:,rindices);
    % project the original data set
%    signals = PC' * data;
else %svd
    Y = data' / sqrt(N-1);
    [u,S,PC] = svd(Y,'econ');
    % calculate the variances
    S = diag(S);
    V = S .* S;
    % project the original data
%    signals = PC' * data;
end