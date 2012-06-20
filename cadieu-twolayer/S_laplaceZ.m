% S_laplace.m - sparse cost function - laplace distribution
%
% function sparse_cost = S_laplace(a,beta,sigma)
%

function sparse_cost = S_laplaceZ(z,beta)

%sparse_cost = beta*sqrt(z.*conj(z));

if numel(beta) == 2
    sparse_cost = [sum(beta(1)*abs(z(1,:))); sum(beta(2)*abs(z(2:end,:)),2)];
else
    sparse_cost = beta*abs(z(:));
end