function [cc mse kern filtSig y kern1] = splitToep(y,x,lags,bin,numCross,sniffTimes)

params.Fs = 1000/bin;params.fpass = [0 20];
win = [.5 .2];params.tapers = [12 win(1) 1];
% [S,t,f] = mtspecgramc(y',win,params);
% idx = kmeans(S,3);h = hist(idx,1:3);[~, h] = max(h);
% idx = idx ~= h;
for i = 1:size(y,1)
    y(i,:) = y(i,:)/max(.01,std(y(i,:)));
end
%figure;plot(linspace(0,numel(y)*bin/1000,numel(y)),y);
%hold all;plot(t,idx,'linewidth',2);
newIdx = zeros(numel(y),1);
newIdx = sum(abs(y),1) > 0;
% fast = find(idx);
% dIdx = [0; diff(idx)];
% i = 1;
% while i < numel(idx)
%     if dIdx(i) ==1 && idx(i+1) == 1
%         start = ceil(t(i)*1000/bin);
%         while dIdx(i+1) ~= -1
%             i = i+1;
%         end
%         stop = ceil(t(i)*1000/bin);
%         newIdx(start:stop) = 1;
%     end
%     i = i+1;
% end
% %plot(linspace(0,numel(y)*bin/1000,numel(y)),newIdx);
% figure;plot(y);hold all;plot(newIdx);
    
% num = 3;fasTime = 250/bin;
% f = find(sniffTimes);
% for i = 1:(numel(f)-num)
% if f(i+num)-f(i) < fasTime*(num)
% newIdx(f(i+1):f(i+num)) = 1;
% end
% end
%

numX = size(x,1);
xx = zeros(size(x,2),lags*size(x,1));
warning off all;
for i = 1:size(x,1)
    x(i,:) = x(i,:)/std(x(i,:));
    xx(:,(i-1)*lags+(1:lags)) = fliplr(toeplitz(x(i,:),zeros(lags,1)));
end
ridge = 1000;
fs = 1000/bin;
warning off all;
newIdx= newIdx(lags:end);newIdx = logical(newIdx);
xx = xx(lags:end,:);
y = y(:,lags:end);ts = 1:numel(y);
%figure;plot(y);hold all;plot(ts(newIdx == 1),y(newIdx == 1),'r')

%y = double(newIdx)';
%figure;plot(y);
%size(y)
params.tapers = [5 9];
%[S f] = mtspectrumc(y,params);(:,1:sum(newIdx));(1:sum(newIdx),:)
xx = xx(newIdx,:);
y = y(:,newIdx);%y = y(2,:);
[cc mse kern] = pipeLine(y,xx,numCross,ridge,fs,numX);
kern = squeeze(mean(kern));
filtSig = kern*xx';
return
kern = squeeze(mean(kern,1));if size(kern,1) == 1 kern = kern';end
[S0 f0] = mtspectrumc(y(~newIdx),params);
xx = xx(newIdx,:);
y = y(:,newIdx);
[S1 f1]= mtspectrumc(y,params);
[S1 f1]= mtspectrumc(y,params);
figure;hold all;plot(f1,S1/max(S1));plot(f0,S0/max(S0));%plot(f,S/max(S));
%kTemp = sin(linspace(0,8*pi,size(xx,2)));
%y = xx * kTemp';y = (y+randn(size(y))/1)'; 
[~,kern1,snr1] = pipeLine(y,xx,numCross,ridge,fs,numX);
kern1 = squeeze(mean(kern1,1));if size(kern1,1) == 1 kern1 = kern1';end
[~,~,snr2] = pipeLine(y,xx,numCross,ridge,fs,numX,kern');
figure;plot(snr1);hold all;plot(snr2);
filtSig = kern.'*xx';