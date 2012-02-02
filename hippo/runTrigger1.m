function runTrigger1(pos,dat,thresh,lags)

Fs = 1250/32;
v = angVel(pos)';
v(isnan(v)) = 0;
vf = filtLow(v(1,:),Fs,.2,1);
x = findpeaks(vf,thresh);
x = x(1).loc;
x(x-lags<1) = [];x((x+lags+1) > size(dat,1)) = [];
%figure;plot(vf);hold all;scatter(x,vf(x),'r','filled');
allDat = zeros(numel(x),2*lags+1,size(dat,2));
allVel = zeros(numel(x),2*lags+1,size(v,1));
for i = -lags:lags
    allDat(:,i+lags+1,:) = dat(x+i,:);
    allVel(:,i+lags+1,:) = v(:,x+i)';
end
%temp = repmat(1:size(allAng,2),[size(allAng,1) 1]);
%t1 = log(squeeze(allDat(:,:,2)./allDat(:,:,1)));
%[h b] = hist3([temp(:) t1(:)],[60 200]);
%figure;imagesc(b{1},b{2},log(h)');
%figure;plot(squeeze(mean(abs(allDat(:,:,1))))*1000);
allDat = squeeze(mean((allDat)));
allVel = squeeze(mean(allVel));
figure;subplot(411);plot(allVel);axis tight;
subplot(412);plot(abs(allDat)/max(abs(allDat)));hold all;plot(real(allDat)/max(real(allDat)));axis tight;
subplot(413);plot(imag(allDat));axis tight;
subplot(414);plot(circ_dist(circ_mean(angle(allDat).'),angle(allDat)));axis tight;
figure;scatter(real(allDat),imag(allDat))
%figure;scatter(allDat(:,1),allFreq(:,1));hold all;
%scatter(allFreq(:,1),allFreq(:,2));
%figure;scatter(allDat(2:end,2),diff(allDat(:,1)));