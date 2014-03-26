function [yX Xcov Xstd samps] = multiRegress(X,sp,pos,file)%,dec)
% linear regression from wavelet-transformed LFP to spike trains
% only uses timepoints when rat is running at the moment...

dec = round(size(X,2)/size(pos,1));
Fs = floor(log2(1250/32*dec/2)*2)/2;
levels = 1:.5:9;
levels(levels > Fs) = [];
2.^[min(levels) max(levels)]

%[pos,~,fast] = fixPos(pos);
if size(pos,1) > size(X,2)
    pos = pos(1:size(X,2),:);
end
[~,~,~,~,~,b] = fixPos(pos);
b = bwmorph(b~=0,'dilate',round(1250/32));
b = resample(double(b),dec,1) > .5;
trial = bwlabel(b);clear b;
%for i =1:4
%    pos1(:,i) = resample(pos(:,i),dec,1);
%end
%if size(pos1,1) > size(X,2)
%    pos1 = pos1(1:size(X,2),:);
%end
%[~,~,~,trial] = fixPos(pos1);
%fast = resample(double(fast),dec,1) > .5;
Xstd = 0;Xcov = single(0); yX = single(0);samps = 0;
for i = 1:max(trial)
   Xtemp = zeros(numel(levels)*size(X,1),sum(trial == i), 'single');
   for j = 1:numel(levels)
       Xtemp((j-1)*size(X,1)+(1:size(X,1)),:) = morFilter(single(X(:,trial == i)),2^(levels(j)),1250/32*dec);
   end
   inds = 1:sum(trial == i);%trial == i;%fast(trial == i)'; %% alternate inds = trial == i;
   Xtemp = [real(Xtemp);imag(Xtemp)];
   Xstd = Xstd + sum(Xtemp(:,inds).*conj(Xtemp(:,inds)),2);
   Xcov = Xcov + Xtemp(:,inds)*Xtemp(:,inds)';
   yX = yX + sp(:,trial == i)*Xtemp(:,inds)';
   samps = samps + sum(inds);
   fprintf('.');
end
Xstd = sqrt(Xstd);
if exist('file','var')
    save([file 'multiRegress.mat'],'yX','Xcov','Xstd','samps');
end

%sz = size(Xcov,1)/2;
%XcovC = complex(Xcov(1:sz,1:sz)+Xcov(sz+1:end,sz+1:end),Xcov(sz+1:end,1:sz)-Xcov(1:sz,sz+1:end));