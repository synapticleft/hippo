% S_laplace.m - sparse cost function - laplace distribution
%
% function sparse_cost = S_laplace(a,beta,sigma)
%

function sparse_cost = S_laplace(a,beta)
if numel(beta) == 2
    sparse_cost = [sum(beta(1)*abs(a(1,:))); sum(beta(2)*abs(a(2:end,:)),2)];
else
    sparse_cost = beta*abs(a(:));
end