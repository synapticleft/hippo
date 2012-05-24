function p = myGprob(x, m, C)
% GAUSSIAN_PROB Evaluate a multivariate Gaussian density.
% p = gaussian_prob(X, m, C)
% p(i) = N(X(:,i), m, C) where C = covariance matrix and each COLUMN of x is a datavector
% If X has size dxN, then p has size Nx1, where N = number of examples
x = abs(x);
if length(m)==1 % scalar
  x = x(:)';
end
[d N] = size(x);
%assert(length(m)==d); % slow
m = m(:);
M = m*ones(1,N); % replicate the mean across columns
denom = (2*pi)^(d/2)*sqrt(abs(det(C)));
%mahal = sum(((x-M)'*inv(C)).*(x-M)',2);   % Chris Bregler's trick
mahal = sum(((x-M)'/C).*(x-M)',2);   % Chris Bregler's trick
%if any(mahal<0)
%  warning('mahal < 0 => C is not psd')
%end
p = -0.5*mahal - log(denom);