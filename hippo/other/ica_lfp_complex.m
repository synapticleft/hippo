function [weights signal zpweights W X] = ica_lfp_complex(X,sources,f1,f2,fs)

X = X - repmat(mean(X,2),[1 size(X,2)]);
X = X./repmat(std(X,1,2),[1 size(X,2)]);

if exist('f1','var') && ~isempty(f1)
    if ~exist('fs','var') || isempty(fs)
        fs = 1250;
    end
    X = hipFilter(X,f1,f2,fs);
end

if isreal(X)%sum(imag(X(:)) == 0)
    X = hilbert(X')';
end
tic;

%%[weights sphere] = runica(signal,'maxsteps',50);%,'interupt','on'

[n,T]	= size(X);
if ~exist('sources','var')
    m=n;
else
    m = sources;
end

%%% whitening
if m<n, %assumes white noise
 	[U,D] 	= eig((X*X')/T); 
	[puiss,k]=sort(diag(D));
 	ibl = ones(m,1);%sqrt(puiss(n-m+1:n)-mean(puiss(1:n-m)));
 	bl 	= ones(m,1) ./ ibl ;
 	W	= diag(bl)*U(1:n,k(n-m+1:n))';
 	IW 	= U(1:n,k(n-m+1:n))*diag(ibl);
    ZP = U(1:n,k(n-m+1:n))*diag(bl)*U(1:n,k(n-m+1:n))';
else    %assumes no noise
 	IW 	= sqrtm((X*X')/T);
 	W	= inv(IW);
    ZP = eye(m);
end;
X	= W*X;
%fast!
[weights signal] = jade_complex(X,m);%

%also slow but converges, noncirc VERY slow
%[weights signal] = ACMNsym(X,'mle_circ');

%also slow; none converge; does converge on whitened?
%[weights signal] = nonCircComplexFastICAsym(X,'sqrt');

%%medium speed; doesnt converge;
%[weights signal] = fast_complex(X);

%%ICA_EBM - slow and doesnt converge
%[weights shat signal] = complex_ICA_EBM(X);

%%%estimation of the mixing matrix and signal separation
signal = weights'*signal;
weights = IW*weights;
zpweights = ZP*weights;
toc