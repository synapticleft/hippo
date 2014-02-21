% dS_cauchy.m - derivative of sparse cost function - cauchy distribution
%
% function sparse_cost = dS_cauchy(u,beta,sigma)
%

function sparse_cost = dS_cauchyZ(z,beta,sigma)
%a = abs(z);
sparse_cost = 2*beta*(1./(1+(abs(z)/sigma).^2)).*z/(sigma^2);%(a./sigma.^2).*z;