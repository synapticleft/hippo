function [kern filtSig y] = toeplags(y,x,lags,bin,numCross,subSet,offSet,ridge)

if ~exist('offSet','var')
    offSet = 0;
end
if ~exist('subSet','var')
    subSet = max(size(y));
end
if ~exist('ridge','var')
    ridge = 1;
end
y = y(:,(1:subSet));
x = x(:,offSet+(1:subSet));
for i = 1:size(y,1)
    y(i,:) = y(i,:)/max(.01,std(y(i,:)));
end
numX = size(x,1);
xx = zeros(size(x,2),lags*size(x,1));
warning off all;
for i = 1:size(x,1)
    x(i,:) = x(i,:)/std(x(i,:));
    xx(:,(i-1)*lags+(1:lags)) = fliplr(toeplitz(x(i,:),zeros(lags,1)));
end
fs = 1000/bin;
%x = decimate(x,down);
%y = decimate(y,down);
warning off all;
xx = xx(lags:end,:);
y = y(:,lags:end);
%[r,s,t,u,v] = canoncorr(y',xx);
[h r] = makeHist(y');
figure;subplot(121);imagesc(r{2},r{1},h);axis image%scatter(u(:,1),u(:,2));
%[h r] = makeHist(v);
%subplot(222);imagesc(r{2},r{1},h);%scatter(v(:,1),v(:,2));
%sPlot([u v]');
[cc,~, kern] = pipeLine(y,xx,numCross,ridge,fs,numX);
kern = squeeze(mean(kern));%if size(kern,1) == 1 kern = kern';end
%kern = reshape(kern,[numel(kern)/size(x,1) size(x,1)]);
%filtSig = kern*xx';
%kern = hilbert(kern);
%sPlot(kern');
filtSig = kern'*xx';
[h r] = makeHist(real(filtSig'));
figure(1);subplot(122);imagesc(r{2},r{1},h);axis image;%scatter(real(filtSig(1,:)),real(filtSig(2,:)));
%[h r] = makeHist([real(filtSig(1,:))' imag(filtSig(1,:)')]);
%subplot(224);imagesc(r{1},r{2},h);
%[h r] = makeHist([real(filtSig(1,:))'-y(1,:)' real(filtSig(2,:))'-y(2,:)']);%[real(filtSig(1,:))' y(1,:)']);
%subplot(222);imagesc(r{2},r{1},h);%scatter(real(filtSig(1,:))',y(1,:)');%
%[h r] = makeHist([real(filtSig(2,:))'-y(2,:)' y(2,:)']);
%subplot(224);imagesc(r{2},r{1},h);
% filtSig = real(filtSig);
% y = complex(y(1,:),y(2,:));bins{2} = -3.2:.1:3.2;bins{1} = 0:.1:3;
% [h r] = makeHist([abs(y)' angle(y)'],bins);
% %hma = mean(h);
% subplot(222);imagesc(r{2},r{1},h);
% filtSig = complex(filtSig(1,:),filtSig(2,:));
% [h r] = makeHist([abs(filtSig)' angle(filtSig)'],bins);
% %hm =mean(h);
% 
% subplot(224);imagesc(r{2},r{1},h);
% figure;sPlot([abs(filtSig(2:end)); abs(y(2:end)); -log(abs(filtSig(2:end))./abs(y(2:end))); getDiffs(filtSig,1); getDiffs(y,1)]);%plot(abs(filtSig));hold all;plot(abs(y));
%figure;subplot(211);plot(sort(abs(filtSig)),sort(abs(y)));%scatter(abs(filtSig),abs(y));
%subplot(212);plot(sort(angle(filtSig)),sort(angle(y)));%scatter(angle(filtSig),angle(y));
%figure;plot(r{2},hm./hma);
%[corr(angle(filtSig)',angle(y)') corr(abs(filtSig)',abs(y)')]
% filtSig1 = [ones(1,size(filtSig,2)); filtSig; filtSig.^2; prod(filtSig); sqrt(sum(filtSig.^2)); atan(filtSig(1,:)./filtSig(2,:))];
% numX = size(filtSig,1);
% [cc kern] = pipeLine(y,filtSig1',numCross,ridge,fs,numX);
% figure;
% % for i= 1:size(filtSig,2)
% %     plot([filtSig(1,i) y(1,i)],[filtSig(2,i) y(2,i)]);hold all;
% % end
% subplot(221);plot(sort(real(filtSig(1,:))),sort(y(1,:)));
% subplot(222);plot(sort(real(filtSig(2,:))),sort(y(2,:)));
% subplot(223);scatter(real(filtSig(1,:)),y(1,:));
% subplot(224);scatter(real(filtSig(2,:)),y(2,:));

function [h x] = makeHist(dat,ctrs)
try
    if ~exist('ctrs','var')
        ctrs{1} = -1:.02:1;ctrs{2} = ctrs{1};
    end
[h x] = hist3(dat(:,1:2),ctrs);h = log(max(1,h));%,[100 100]
catch
    h = 1;x{1} = 1;x{2} = 1;
end