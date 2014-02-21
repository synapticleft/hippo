% S_cauchy.m - sparse cost function - cauchy distribution
%
% function sparse_cost = S_cauchy(u,beta,sigma)
%

function sparse_cost = S_cauchyZ(a,beta,sigma)

sparse_cost = beta*log(1+(abs(a)/sigma).^2);