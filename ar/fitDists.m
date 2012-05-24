function [A, noise, x, V] = fitDists(fileNum,modelSize,figsOn,numInds)

load_data;
spikes = spktrain;
v = signals;
if ~exist('numInds','var');
    numInds = size(v,2);
end
filt = 1;
sum(spikes(:))
spikes = [spikes; zeros(size(v,1)-size(spikes,1),size(spikes,2))];
spikes = spikes(:,1:numInds);
if filt 
    v = hipFilter(v',.01,200,10000)';
end
cMor = getMor(fm,sf,1000,5);
for i= 1:numInds
    z(i,:) = conv(v(:,i),cMor,'same');
    [coeffs(i,:) noise(i)] = arburg(z(i,:),modelSize);
%    [~, p(i,:)] = residue([1 zeros(1,modelSize-1)],coeffs(i,:));
end
%v = glueAll(v);
%coeffs(end+1,:) = coeffs(1,:);noise(end+1) = noise(1);
coeffs(end+1,:) = mean(coeffs);noise(end+1) = mean(noise);
%[~,p(end+1,:)] = residue([1 zeros(1,modelSize-1)],coeffs(end,:));
for i = 1:size(z,1)
    thisNoise = makeComplex(numel(v(:,1)),sqrt(noise(end)),[1 1]);
    simSig(i,:) = arData(thisNoise,coeffs(end,:));
end
if figsOn
    figure;plot(real(z(1,:)));hold all;plot(real(simSig(1,:)));
    fftz = mean(abs(fft(z,[],2)));fftAR = mean(abs(fft(simSig,[],2)));
    figure;plot(log(fftz));hold all;plot(log(fftAR));
    zAng = getDiffs(z,1);arAng = getDiffs(simSig,1);
    xs = 0:.01:pi/4;
    figure;plot(xs,log(hist(zAng(:),xs)));hold all;plot(xs,log(hist(arAng(:),xs)));
    xs = 0:.01:1;
    figure;plot(xs,log(hist(abs(z(:)),xs)));hold all;plot(xs,log(hist(abs(simSig(:)),xs)));
end
spikes = logical(spikes)';
z_fit = gmdistribution.fit([real(z(:)) imag(z(:))],1);
z_spike = z(spikes);
range = 2*std(z(:));
bounds{1} = linspace(-range,range,40);bounds{2} = bounds{1};
pZ_Data = hist3([real(z(:)) imag(z(:))],bounds);
pZ_Data = pZ_Data/sum(pZ_Data(:));
pZ_spike_Data = hist3([real(z_spike) imag(z_spike)],bounds);
pZ_spike_Data = pZ_spike_Data/sum(pZ_spike_Data(:));
pSpike_z_Data = pZ_spike_Data./pZ_Data;

z_spike_fit = gmdistribution.fit([real(z_spike) imag(z_spike)],1);
[x y] = meshgrid(bounds{1});
pZ_spike_Fit = reshape(pdf(z_spike_fit,[x(:) y(:)]),[numel(bounds{1}) numel(bounds{1})]);
pZ_Fit = reshape(pdf(z_fit,[x(:) y(:)]),[numel(bounds{1}) numel(bounds{1})]);
pZ_Fit = pZ_Fit/sum(pZ_Fit(:));
pZ_spike_Fit = pZ_spike_Fit/sum(pZ_spike_Fit(:));
pSpike_z_Fit = pZ_spike_Fit./pZ_Fit;

myFit.mu = mean(z_spike);
myFit.angle = angle(myFit.mu);
z_spike_rotate = z_spike*exp(1i*-myFit.angle);
myFit.varReal = var(real(z_spike_rotate));
myFit.varImag = var(imag(z_spike_rotate));
myFit.rotate = [cos(myFit.angle) -sin(myFit.angle); sin(myFit.angle) cos(myFit.angle)];
myFit.Sigma = myFit.rotate*[myFit.varReal 0;0 myFit.varImag]*myFit.rotate';
myFit = gmdistribution([real(myFit.mu) imag(myFit.mu)],myFit.Sigma);
pZ_spike_myFit = reshape(pdf(myFit,[x(:) y(:)]),[numel(bounds{1}) numel(bounds{1})]);
pZ_spike_myFit = pZ_spike_myFit/sum(pZ_spike_myFit(:));
pSpike_z_myFit = pZ_spike_myFit./pZ_Fit;
z_fit = gmdistribution.fit(z(:),1);
%myFit = gmdistribution.fit(z_spike,1);
%z = z(:);spikes = spikes(:);
z = z.'; spikes = spikes.';
temp = zeros(modelSize,sum(spikes(:)));
for i = 1:modelSize
    temp(i,:) = z(find(spikes)-i+1);
end
spikeFit.mu = mean(temp,2);
spikeFit.Sigma = cov(temp.');
spikeTimes = getTimes1(spikes)'+1;
[spikeFit.mu angle(spikeFit.mu)-angle(spikeFit.mu(1)) abs(spikeFit.mu)]
%spikeFit.mu = stableZ(spikeTimes,z,20,modelSize);
%spikeFit.mu = weightedW(spikes,z,modelSize);
%[spikeFit.mu angle(spikeFit.mu)-angle(spikeFit.mu(1)) abs(spikeFit.mu)]
spikeTimes1 = circshift(spikeTimes,[0 -modelSize+1]);
z = z.'; z = z(1:size(spikeTimes,1),:);
%spikeFit.mu = accumarray(spikeTimes1(:),z(:),[],@mean);
%spikeFit.mu = flipud(spikeFit.mu(1:modelSize));
angDiff = zeros(size(spikeTimes));absDiff = angDiff;zDiff = angDiff;Vs = angDiff;
for i = 1:size(spikeTimes,1)
    [x, V, A] = complexKalman(spikes.',spikeFit,z_fit,coeffs(end,:),noise(end),i,z);
    angDiff(i,:) = circ_dist(z(i,:),x(1,:));
    absDiff(i,:) = log(abs(z(i,:))./abs(x(1,:)));
    zDiff(i,:) = z(i,:) - x(1,:);
    xAll(i,:) = x(1,:);
    Vs(i,:) = squeeze(V(1,2,:));
end
figure;imagesc(imHist(spikeTimes,angDiff,(-pi:.01:pi)));
figure;scatter(real(z(i,:)),real(x(1,:)));
figure;plot(accumarray(spikeTimes(:),(angDiff(:)),[],@mean));
%figure;plot(accumarray(spikeTimes(:),absDiff(:),[],@mean));
figure;plot(accumarray(spikeTimes(:),log(abs(zDiff(:))),[],@mean));

figure;hold all;for i = 100:-1:1
    f = find(spikeTimes == i);
    ff(i) = numel(f);
    scatter(real(zDiff(f)),imag(zDiff(f)));
    zDm(i) = mean(zDiff(f));zs(i) = var(zDiff(f));
    zm(i) = mean(z(f));xm(i) = mean(xAll(f));
end
plot(zDm,'k','LineWidth',2);
figure;plot(zDm);hold all;plot(zm);plot(xm);plot(spikeFit.mu,'LineWidth',2);
figure;plot(accumarray(spikeTimes(:),log(real(Vs(:))),[],@mean));
hold on;plot(accumarray(spikeTimes(:),log(imag(Vs(:))),[],@mean),'r');plot(log(zs),'k');
%figure;plot(accumarray(spikeTimes(:),angle(Vs(:)),[],@mean));%abs
figure;plot(real(zDm));hold all;plot(imag(zDm));
figure;plot(real(z(size(spikeTimes,1),:)));hold on;plot(real(x(1,:)),'r');hold on;
scatter(find(spikes(:,size(spikeTimes,1))),zeros(1,sum(spikes(:,size(spikeTimes,1)))),'k');
figure;plot(log(abs(squeeze(V(1,1,:)))));


if figsOn
    figure;subplot(3,3,1);imagesc(log(pZ_Data));axis image;
    range1 = [min(log(pZ_Data(pZ_Data(:)>0))) max(log(pZ_Data(:)))];
    subplot(3,3,2);imagesc(log(pZ_spike_Data'));axis image;
    range2 = [min(log(pZ_spike_Data(pZ_spike_Data(:)>0))) max(log(pZ_spike_Data(:)))];
    subplot(3,3,3);imagesc(log(pSpike_z_Data'));axis image
    range3 = [min(log(pSpike_z_Data(pZ_Data(:)>0 & pZ_spike_Data(:)>0))) max(log(pSpike_z_Data(:)))];
    subplot(3,3,4);imagesc(log(pZ_Fit),range1);axis image;
    subplot(3,3,5);imagesc(log(pZ_spike_Fit),range2);axis image;
    subplot(3,3,6);imagesc(log(pSpike_z_Fit),range3);axis image;
    subplot(3,3,7);imagesc(log(pZ_Fit),range1);axis image;
    subplot(3,3,8);imagesc(log(pZ_spike_myFit),range2);axis image;
    subplot(3,3,9);imagesc(log(pSpike_z_myFit),range3);axis image;
    xs = -pi:.1:pi;
    figure;plot(xs,log(hist(angle(z_spike),xs)/numel(z_spike)));hold all;
    c = complex(x,y);
    c = angle(c);c = c-min(c(:))+.001; c = c/max(c(:)); c = ceil(c*numel(xs));
    temp = accumarray(c(:),pZ_spike_myFit(:),[],@mean);
    temp = temp/sum(temp);
    plot(xs,log(temp));
end

function z = stableZ(spikes,zs,back,p)
[r c] = find(spikes == back);
c = c - back;
z = zeros(p,1);
for i = 1:numel(r)
    if (c(i)-back-p > 0)
        z = z + zs(c(i)-back-(1:p)+1,r(i));
    end
end
z = z/numel(r);

function z = weightedW(spikes,zs,p)
zs = zs.';spikes = spikes.';
z = zeros(p,1);sum = 0;
for i = 1:size(zs,1)
    f = find(spikes(i,:)==1);
    f(end+1) = size(zs,2);
    f(f<p) = [];
    %sum = sum + f(end) - f(1);
    for j = 1:(numel(f)-1)
        factor = min(f(j+1)-f(j),inf);
        z = z + (zs(i,f(j)-(1:p)+1)*factor).';
        sum = sum + factor;
    end
end
z = z/sum;
        

function im = imHist(spikeTimes,data,bins)
spikeTimes = spikeTimes(:);data = data(:);
spikeTimes = min(50,ceil(spikeTimes/2));
range = max(spikeTimes);
im = zeros(range,numel(bins));
for i = 1:range
    im(i,:) = hist(data(spikeTimes == i),bins);
    im(i,:) = im(i,:)/sum(im(i,:));
end
im = log(max(.00001,im));

function t = getTimes1(spikes)
spikes = spikes';
s = size(spikes);
t = zeros(size(spikes));
for i = 1:s(1)
    counter = 500;
    for j= 1:s(2)
        if spikes(i,j)
            counter = 0;
        else
            counter = counter + 1;
        end
        t(i,j) = counter;
    end
end
t = t';

function t = getTimes(spikes)
s = size(spikes); spikes = spikes(:);
counter = 500;
t = zeros(size(spikes));
for i = 1:numel(spikes)
    if spikes(i)
        counter = 0;
    else
        counter = counter + 1;
    end
    t(i) = counter;
end
t = reshape(t,s);

function c = makeComplex(len,dev,ratio)
tot = sum(ratio.^2);
c = dev/sqrt(tot)*complex(randn(len,1)*ratio(1),randn(len,1)*ratio(2));

function y = arData(x,b)
y = x;
for i = 1:numel(x)
    for j = max(1,i-numel(b)+1):(i-1)
        y(i) = y(i) - y(j)*b(i-j+1);
    end
end

function dPhase = getDiffs(sig,gap)
angleA = angle(sig);
if min(size(sig)) == 1
    dPhase = circ_dist(angleA((gap+1):end),angleA(1:(end-gap)));
else
    dPhase = circ_dist(angleA(:,(gap+1):end),angleA(:,1:(end-gap)));
end

% function sig = glueAll(sig)
% temp = sig(:,1);
% %figure;
% for i = 2:size(sig,2)
%     %plot([temp(end-200:end); sig(1:200,i)]);pause(1);%[temp(end-2:end) sig(1:3,i)]
%     if temp(end) == sig(2,i)
%         temp = [temp; sig(3:end,i)];
%     else
%         temp = [temp; sig(2:end,i)];
%     end
% end
% figure;plot(sig(:));hold all;plot(temp);
% sig = temp;