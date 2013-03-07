function [X whiteningMatrix dewhiteningMatrix] = getData(fname,dec,elecs,sz,rind)%,rdim)
% USAGE  >> X = getData('ec014.468.h5',1:64,[],[],32); %32 is the
% decimation factor, you may also use a smaller one to get high-freq
% information. You will need to change data_root to match the folder where
% you store the data files.

%data_root = '/media/Expansion Drive/redwood/';%'/media/work/hippocampus/';%
data_root1 = '/media/Expansion Drive/redwood/';%KenjiMizuseki/';
%fname = 'hippo.h5';%'96elec.h5';%
%padding = 0;%20000;
%info = hdf5info([data_root fname '.h5']);
%nSamples = info.GroupHierarchy.Datasets(1).Dims;
a = memmapfile([data_root1 fname '.eeg'],'Format','int16');
%nSamples(1) = 64;
%if exist('numChans','var')
%    [~,~,nSamples,numChans] = LoadBinary(input,1,numChans,[],[],[],[1 2]);
%else
%[~,~,nSamples,numChans] = LoadBinary([data_root1 fname '.eeg'],1,[],[],[],[],[1 2]);
%nSamples = [numChans nSamples];
%end
if ~exist('dec','var') || isempty(dec)
    dec = 1;
end
if  ~exist('elecs','var') || isempty(elecs)
    Par = LoadPar([data_root1 fname]);
    nSamples = [Par.nChannels numel(a.data)/Par.nChannels];
    elecs = 1:nSamples(1);
else
    nSamples(1) = max(elecs);
end
if ~exist('sz','var') || isempty(sz)
    sz = nSamples(2);
    sz = dec*(floor(sz/dec));
end
if ~exist('rind','var') || isempty(rind)
    rind = 0;
end
%rind = padding + ceil(rand*(nSamples(2) - padding*2 - sz));

%chunk = complex(double(h5varget([data_root fname],'/hReal',[0 rind-1],[nSamples(1) sz])),...
%      double(h5varget([data_root fname],'/hImag',[0 rind-1],[nSamples(1) sz])));
%X = double(h5varget([data_root fname],'/hReal',[0 rind],[nSamples(1) sz]));%-1

%if dec > 1
    X = zeros(numel(elecs),ceil(sz/dec));
    for i = 1:numel(elecs)
        %tic;
        %temp = double(h5varget([data_root fname '.h5'],'/hReal',[elecs(i)-1 rind],[1 sz]));%nSamples(1)
        %toc
        tic;
        temp = double(a.data((rind:sz-1)*nSamples(1)+elecs(i)))';
        toc
        if dec > 1
            tic;
            X(i,:) = decimate(temp,dec);%X(i,:),dec);
            toc
        else
            X(i,:) = temp;
        end
        fprintf([num2str(i) ' ']);
    end
    %X = X1;clear X1;
%else
%    X = zeros(numel(elecs),sz);
%    for i = 1:numel(elecs)
%        X(i,:) = double(h5varget([data_root fname],'/hReal',[elecs(i)-1 rind],[1 sz]));
%    end
%end
X = bsxfun(@minus,X,mean(X,2));
X = bsxfun(@rdivide,X,std(X,0,2));

% if nargout > 1
%     % Calculate the eigenvalues and eigenvectors of covariance matrix.
%     fprintf ('Calculating covariance...\n');
%     covarianceMatrix = X*X'/size(X,2);
%     [E, D] = eig(covarianceMatrix);
%     
%     % Sort the eigenvalues and select subset, and whiten
%     fprintf('Reducing dimensionality and whitening...\n');
%     [~,order] = sort(diag(-D));
%     E = E(:,order(1:rdim));
%     d = diag(D);
%     d = real(d.^(-0.5));
%     D = diag(d(order(1:rdim)));
%     X = D*E'*X;
%     
%     whiteningMatrix = D*E';
%     dewhiteningMatrix = E*D^(-1);
% end