function [Z,V,E,D] = zca2(X)
epsilon = 1e-10;
% Calculate the eigenvalues and eigenvectors of the new covariance matrix.
covarianceMatrix = X*X'/size(X,2);
% [E, D] = eig(covarianceMatrix);
% 
% % Sort the eigenvalues  and recompute matrices
% [~,order] = sort(diag(-D));
% E = E(:,order);
% d = diag(D); 
% dsqrtinv = real((d + eps).^(-0.5)); %ones(size(d))*
% Dsqrtinv = diag(dsqrtinv(order));
% D = diag(d(order));
% V = E*Dsqrtinv*E';
% Z = V*X;


[U,S] = svd(covarianceMatrix);
V = U * diag(1./sqrt(diag(S) + epsilon)) * U';
Z = V * X;