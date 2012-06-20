% dS_laplace.m - derivative sparse cost function - laplace distribution
%
% function sparse_cost = dS_laplace(a,beta)
%

function sparse_cost = dS_laplaceZ(z,beta)

sparse_cost = z./sqrt(z.*conj(z));
if beta == 0
    sparse_cost(:) = 0;
end