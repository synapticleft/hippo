function [f,g] = objfun_a_conv(a0,I,phi,lambda)

prior = 'laplace';%'gaussian';%'cauchy';%
[N J R] = size(phi);
 
S = size(I,2);
a = reshape(a0, J, S+R-1);
 
EI = zeros(N,S);
%EI1 = [];
%for i = 1:N
%    EI(i,:) = sum(fastConv(a,squeeze(phi(i,:,:))));
%end
%figure;plot(EI1');
for t = 1:R
    EI = EI + phi(:,:,R-t+1)*a(:,t:t+S-1);
end
%figure;plot(EI');return
E = I - EI;

f_residual = 0.5*sum(E(:).^2);
switch prior
    case 'laplace'
        f_sparse = lambda(1)*sum(abs(a(:)));
    case 'cauchy'
        f_sparse = lambda(1)*sum(log(1+(a(:)/lambda(2)).^2));
    case 'gaussian'
        f_sparse = lambda(1)*sum(a(:).^2);
end

%fprintf('f_residual %.4f f_sparse %.4f\n', f_residual, f_sparse);

f = f_residual + f_sparse;

da = zeros(size(a));
for t = 1:R
    srt = R+1-t;
    fin = S+R-t;

    da(:,srt:fin) = da(:,srt:fin) - phi(:,:,t)'*E;
end


if 0
    figure(8);
    imagesc(a); colorbar;
    drawnow

    figure(11); clf;
        
    X_n = reshape(I, S, N);
    EI_n = reshape(EI, S, N);
    E_n = reshape(E, S, N);

    mn = min([I(:) ; EI(:) ; E(:)]);
    mx = max([I(:) ; EI(:) ; E(:)]);
    for n = 1:N
        subp(N,1,n);
        plot(1:S,X_n(:,n),'b-', ...
             1:S,EI_n(:,n),'g-.', ...
             1:S,E_n(:,n),'r-');
        axis([1 S mn mx]);
        legend('Source', 'Estimate', 'Error');
        title(sprintf('Source %d', n));
    end
end

switch prior
    case 'laplace'
        d_sparse = lambda(1)*sign(a);
    case 'cauchy'
        d_sparse = lambda(1)*2*(1./(1+(a/lambda(2)).^2)).*(a./lambda(2).^2);
    case 'gaussian'
        d_sparse = lambda(1)*2*a;
end

da = da + d_sparse;

g = da(:);

