function [S, vopts] = GSIR(X, Y, d, s)
% Input: 
% X: p X n input data matrix
% Y: response variable
% d: number of LSIR directions
% s: regularization parameter
% opts: structrue containing parameters
%       [pType]: Type of problem: 'c' for classification
%                               'r' for regression
%       [H]: number of slices, defalt: 10 for 'r' and class number for 'c'
%       [numNN]: number of nearest neighbors
%
% Output:
% S:   a structure containing
%      [edrs] : edr directions
%      [Xmean]: mean of the samples
%      [Xv]: centered LSIR variates for input data X
%      [Cov]: Covariance matrix of inverse regression
%      [Sigma]: Covariance matrix of X
% vopts: structure of parameters used
%
% last modified: 9/16/2008 (documentation)
H = 10;

[dim n] = size(X);
[~, YI] = sort(Y);
Hn = round(n/H);
cX = bsxfun(@minus,X,mean(X,2));
Sigma = zeros(size(X,1));
st = std(cX(:))*20;
c = colormap;
for i = 1:H
    if i<H
        ind = YI((i-1)*Hn+1:i*Hn);
    else
        ind =YI((H-1)*Hn+1:n);
    end
    Xi = cX(:,ind);
    [u,s,v] = svds(bsxfun(@minus,Xi,mean(Xi,2)),1);
    subplot(211);plot(ind,Xi,'.');subplot(212);%plot(ind,v,'.');pause(.1);
    %plot(abs(u),'color',c(i*6,:));hold all;pause(.1);
    [h,xout] = hist(v',linspace(-3*std(v),3*std(v),100));%linspace(-st,st,100));%
    [~,mins] = peakdet(filtfilt(gausswin(20),1,h),1);
    hs(i,:) = filtfilt(gausswin(30),1,h);
   if isempty(mins)
        mins = [1 numel(h)];
    else
        mins = [1 mins(:,1).' numel(h)];
    end
    xout(1) = min(v);xout(end) = max(v);%v = v(randperm(numel(v)));
    for j = 2:numel(mins)
        inda = v > xout(mins(j-1)) & v < xout(mins(j));
        subplot(212);plot(ind(inda),v(inda),'.');hold all;
        if sum(ind)
            Xm = mean(Xi(:,inda),2);
            %plot(Xm,'color',c(i*6,:));hold all;pause(.1);
            Sigma = Sigma + Xm*Xm'*sum(inda);
        end
    end
    pause(1);hold off;
    %for j=1:ni
    %    dist2j = sum((Xi-repmat(Xi(:,j),1,ni)).^2);
    %    [~, dI] = sort(dist2j);
    %    J(ind(dI(1:numNNi)),ind(j)) = 1/numNNi;
    %end
end
figure;imagesc(hs);return
%figure;subplot(211);imagesc(Sigma);
% % figure;imagesc(h1);
% % for i = 1:size(h1,1)
% %     [maxes,mins] = peakdet(h1(i,:),1);
% %     hold on;scatter(maxes(:,1),i*ones(size(maxes,1),1),'g','filled');
% %     hold on;scatter(mins(:,1),i*ones(size(mins,1),1),'r','filled');
% % end
% subplot(212);imagesc(hs);

%J = J'*J;
%cX = X - repmat(Xmean, 1, n);
%Sigma = cX*J*cX';

Xmean = mean(X, 2);
Sigma = .5*(Sigma + Sigma');
Cov = cX*cX';
Cov = .5*(Cov + Cov');
[B L] = eig(Sigma, Cov + s*eye(dim));
[~ ,LI] = sort(diag(L),'descend');
B = B(:,LI(1:d));
for i = 1:d
    B(:,i) = B(:,i)/norm(B(:,i));
    [~, maxi] = max(abs(B(:,i)));
    B(:,i) = B(:,i)*sign(B(maxi,i));
end

S.edrs = B;
S.Xmean = Xmean;
S.Xv = B'*cX;
S.Cov = Cov;
S.Sigma = Sigma;