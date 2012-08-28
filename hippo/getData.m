function [X whiteningMatrix dewhiteningMatrix] = getData(fname,sz,rind,dec,rdim,elecs)

data_root = '/media/work/hippocampus/';%'/media/Expansion Drive/redwood/';%
%fname = 'hippo.h5';%'96elec.h5';%
padding = 0;%20000;
info = hdf5info([data_root fname]);
nSamples = info.GroupHierarchy.Datasets(1).Dims;
%nSamples(1) = 64;
if ~exist('dec','var') || isempty(dec)
    dec = 1;
end
if isempty(sz)
    sz = nSamples(2)
    sz = dec*(floor(sz/dec))
end
%elseif exist('dec','var')
%    sz = sz*dec;
%end
if  ~exist('elecs','var') || isempty(elecs)
    elecs = 1:nSamples(1);
end
   
if ~exist('rind','var') || isempty(rind)
    rind = 0;
end
%rind = padding + ceil(rand*(nSamples(2) - padding*2 - sz));

%chunk = complex(double(h5varget([data_root fname],'/hReal',[0 rind-1],[nSamples(1) sz])),...
%      double(h5varget([data_root fname],'/hImag',[0 rind-1],[nSamples(1) sz])));
%X = double(h5varget([data_root fname],'/hReal',[0 rind],[nSamples(1) sz]));%-1

if dec > 1
    X1 = zeros(numel(elecs),round(sz/dec));
    for i = elecs
        tic;
        temp = double(h5varget([data_root fname],'/hReal',[i-1 rind],[1 sz]));%nSamples(1)
        toc
        tic;
        X1(i-min(elecs)+1,:) = decimate(temp,dec);%X(i,:),dec);
        toc
        fprintf([num2str(i-min(elecs)+1) ' ']);
    end
    X = X1;clear X1;
else
    X = double(h5varget([data_root fname],'/hReal',[0 rind],[nSamples(1) sz]));
end
X = bsxfun(@minus,X,mean(X,2));
X = bsxfun(@rdivide,X,std(X,0,2));

if nargout > 1
    % Calculate the eigenvalues and eigenvectors of covariance matrix.
    fprintf ('Calculating covariance...\n');
    covarianceMatrix = X*X'/size(X,2);
    [E, D] = eig(covarianceMatrix);
    
    % Sort the eigenvalues and select subset, and whiten
    fprintf('Reducing dimensionality and whitening...\n');
    [~,order] = sort(diag(-D));
    E = E(:,order(1:rdim));
    d = diag(D);
    d = real(d.^(-0.5));
    D = diag(d(order(1:rdim)));
    X = D*E'*X;
    
    whiteningMatrix = D*E';
    dewhiteningMatrix = E*D^(-1);
end