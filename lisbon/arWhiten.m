function [X,V] = arWhiten(X,inds)
%Whiten matrix by finding AR residuals

sz = size(X,1)/inds;
refInds = zeros(1,size(X,1));
refInds(1:sz:end) = 1;refInds = logical(refInds);
V = X(refInds,:)/X(~refInds,:);
X = X(refInds,:) - V*X(~refInds,:);