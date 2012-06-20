function [A, noise, R] = fitDists3(fileNum,modelSize,numInds)
sub = 1;fs = 1000/sub;
load_data;
spikes = spktrain;
v = signals;
if ~exist('numInds','var');
    numInds = size(v,2);
end
spikes = [spikes; zeros(size(v,1)-size(spikes,1),size(spikes,2))];
spikes = spikes(:,1:numInds);
cMor = getMor(fm,sf,fs,5);z = v;
for i= 1:numInds
    v(:,i) = conv(v(:,i),cMor,'same');
    spikes(:,i) = conv(spikes(:,i),cMor,'same');
end
s = spikes(2:end,:).*exp(1i*-angle(spikes(1:end-1,:)));
spikes = spikes.*exp(1i*-angle(v));
v = v(2:end,:).*exp(1i*-angle(v(1:end-1,:)));
vm = mean(v,2);spm = mean(spikes,2);sm = mean(s,2);
figure;%plot(real(vm)./mean(abs(v),2));hold all;plot(imag(vm)./mean(abs(v),2));plot(real(sm)./mean(abs(spikes),2));plot(imag(sm)./mean(abs(spikes),2));
vm = angle(vm);vm = (vm-mean(vm))/std(vm);sm = angle(sm); sm = (sm-mean(sm))/std(sm);spm = angle(spm); spm = (spm-mean(spm))/std(spm);
cs = cumsum(sm); cs = cs-mean(cs);cs = cs/std(cs);
plot(vm);hold all;plot(sm);plot(spm);plot(cs);
%spikes = filtLow(spikes',fs,3)';
%v = filtLow(v',fs,3)';
figure;imagesc(complexIm(spikes',0,.5,5));
figure;imagesc(complexIm(v',0,1,50));