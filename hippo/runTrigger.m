function runTrigger(pos,dat,thresh,lags)

Fs = 1250/32;
v = angVel(pos)';
v(isnan(v)) = 0;
vf = filtLow(v(1,:),Fs,.2,1);
x = findPeaks(vf,thresh);
x = x(1).loc;
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
allDat = squeeze(mean(abs(allDat)));
allDat = filtLow(allDat',Fs,.5,1)';
allAng = circ_mean(allAng);
allFreq = squeeze(circ_mean(allFreq));
allFreq = filtLow(allFreq',Fs,.5,1)';
allVel = squeeze(mean(allVel));
figure;subplot(411);plot(allVel);axis tight;
subplot(412);plot(allDat);%renorm(allDat')');
hold all;
%plot(diff(allDat(:,1)));
axis tight;
subplot(413);plot(allAng);axis tight;
subplot(414);plot(allFreq);axis tight;
figure;scatter(allDat(:,1),allFreq(:,1));hold all;
%scatter(allFreq(:,1),allFreq(:,2));
%figure;scatter(allDat(2:end,2),diff(allDat(:,1)));