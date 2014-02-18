function [f, g] = objfun_a(x0, phi, I, lambda)

prior = 'cauchy';%'laplace';%

[L M] = size(phi);

B = size(I,2);
a = reshape(x0,M,B);

E = I - phi*a;

f_recon = 0.5*sum(E(:).^2);

switch prior
    case 'laplace'
        f_sparse = lambda(1)*sum(abs(a(:)));
    case 'cauchy'
        f_sparse = lambda(1)*sum(log(1+(a(:)/lambda(2)).^2));
    case 'gaussian'
        f_sparse = lambda(1)*sum(a(:).^2);
end

%f_sparse = lambda*sum(abs((a(:))));

f = f_recon + f_sparse;

switch prior
    case 'laplace'
        d_sparse = lambda(1)*sign(a);
    case 'cauchy'
        d_sparse = lambda(1)*2*(1./(1+(a/lambda(2)).^2)).*(a./lambda(2).^2);
    case 'gaussian'
        d_sparse = lambda(1)*2*a;
end

df = -(phi'*E) + d_sparse;%lambda*sign(a);
g = df(:);
