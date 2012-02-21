function runTrigger(pos,dat,thresh,lags)

Fs = 1250/32;
v = angVel(pos)';
v(isnan(v)) = 0;
vf = filtLow(v(1,:),Fs,.2);
[~, x] = findpeaks(vf,'MINPEAKHEIGHT',thresh);
%figure;plot(vf);hold all;scatter(x,vf(x));return
%x = x(1).loc;
x(x-lags<1) = [];x((x+lags+1) > size(dat,1)) = [];
%figure;plot(vf);hold all;scatter(x,vf(x),'r','filled');
allDat = zeros(numel(x),2*lags+1,size(dat,2));
allVel = zeros(numel(x),2*lags+1,size(v,1));
allAng = zeros(numel(x),2*lags+1);
allFreq = zeros(numel(x),2*lags+1,size(dat,2));
for i = -lags:lags
    allDat(:,i+lags+1,:) = dat(x+i,:);
    allVel(:,i+lags+1,:) = v(:,x+i)';
    allAng(:,i+lags+1) = circ_dist(angle(dat(x+i,1)),angle(dat(x+i,2)));
    allFreq(:,i+lags+1,:) = circ_dist(angle(dat(x+i,:)),angle(dat(x+i+1,:)));
end
temp = repmat(1:size(allAng,2),[size(allAng,1) 1]);
t1 = allAng;%log(squeeze(allDat(:,:,2)./allDat(:,:,1)));
[h b] = hist3([temp(:) t1(:)],[60 50]);
figure;imagesc(b{1},b{2},(h)');
figure;plot(squeeze(mean(abs(allDat(:,:,1))))*1000);hold all;plot(circ_std(allAng))
allAng = circ_mean(allAng);
%temp = repmat(1:size(allFreq,2),[size(allFreq,1) 1]);
%t1 = squeeze((allFreq(:,:,1)));
%[h b] = hist3([temp(:) t1(:)],[30 100]);
%figure;imagesc(b{1},b{2},(h)');
allDat = squeeze(mean(abs(allDat)));
allDat = filtLow(allDat',Fs,5,1)';
allFreq = squeeze(circ_mean(allFreq));
allFreq = filtLow(allFreq',Fs,.5,1)';
allVel = squeeze(mean(allVel));
figure;subplot(311);plot(allVel(:,1),'linewidth',2);axis tight;
set(gca,'fontsize',16);title 'velocity';
subplot(3,1,2:3);imagesc(b{1},b{2},h');
set(gca,'fontsize',16);title 'angular distance PC1 & 2'
return
subplot(412);plot(allDat);axis tight;
subplot(413);plot(allAng);axis tight;
subplot(414);plot(allFreq);axis tight;
%figure;scatter(allDat(:,1),allFreq(:,1));hold all;
%scatter(allFreq(:,1),allFreq(:,2));
%figure;scatter(allDat(2:end,2),diff(allDat(:,1)));