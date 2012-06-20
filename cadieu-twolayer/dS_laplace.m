% dS_laplace.m - derivative sparse cost function - laplace distribution
%
% function sparse_cost = dS_laplace(a,beta)
%

function sparse_cost = dS_laplace(a,beta)
if numel(beta) == 2
    sparse_cost = [beta(1)*sign(a(1,:)); beta(2)*sign(a(2:end,:))];
else
    sparse_cost = beta*sign(a);
end