% dS_laplace.m - derivative sparse cost function - laplace distribution
%
% function sparse_cost = dS_laplace(a,beta)
%

function sparse_cost = dS_l1l2(a,beta)

sparse_cost = beta*bsxfun(@rdivide,a,sqrt(sum(a.*conj(a),2)));