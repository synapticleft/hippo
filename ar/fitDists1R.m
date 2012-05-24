function [A, noise, R, x, V] = fitDists1R(fileNum,modelSize,figsOn,numInds)

load_data;
spikes = spktrain;
v = signals;
if ~exist('numInds','var');
    numInds = size(v,2);
end
sum(spikes(:))
spikes = [spikes; zeros(size(v,1)-size(spikes,1),size(spikes,2))];
spikes = logical(spikes(:,1:numInds));
v = hipFilter(v',.01,200,10000)';
cMor = getMor(fm,sf,1000,5);
for i= 1:numInds
    z(:,i) = real(conv(v(:,i),cMor,'same'));
    [coeffs(i,:) noise(i)] = arburg(z(:,i),modelSize);
    zPhase = angle(z(:,i));
    [cPhase nPhase] = arburg(zPhase,modelSize);
end
x1 = xcorr(z(:,i),500,'coeff');
% figure;plot(x1);hold all;
% coeffs(end+1,:) = mean(coeffs);noise(end+1) = mean(noise);
% for i = 1:8
%     [c n] = arburg(z(:,end),i);
%     simSig = arData(sqrt(n)*randn(1,numel(z)),c);
%     x2 = xcorr(simSig,500,'coeff');
%     plot(x2);
% end
log(var(z(:)))
temp = zeros(modelSize,sum(spikes(:)));
for i = 1:modelSize
    temp(i,:) = z(find(spikes)-i+1);
end
spikeFit.Sigma = cov(temp');
spikeFit.mu = mean(temp,2);
spikeTimes = getTimes1(spikes)'+1;
z = z.'; z = z(1:size(spikeTimes,1),:);
angDiff = zeros(size(spikeTimes));absDiff = angDiff;zDiff = angDiff;Vs = angDiff;
for i = 1:size(spikeTimes,1)
    [x, V, A,R] = complexKalman(spikes.',spikeFit,coeffs(end,:),noise(end),i,z);
    %angDiff(i,:) = circ_dist(angle(z(i,:)),angle(x(1,:)));
    zDiff(i,:) = z(i,:) - x(1,:);
    xAll(i,:) = x(1,:);
    Vs(i,:) = squeeze(V(1,1,:));
end
R = spikeFit.Sigma;
figure;plot(z(size(spikeTimes,1),:));hold on;plot(x(1,:),'r');hold on;
scatter(find(spikes(:,size(spikeTimes,1))),zeros(1,sum(spikes(:,size(spikeTimes,1)))),'k');
%figure;scatter(real(z(i,:)),real(x(1,:)));
%figure;plot(accumarray(spikeTimes(:),(angDiff(:)),[],@mean));
%figure;plot(accumarray(spikeTimes(:),log(abs(zDiff(:))),[],@mean));

%figure;hold all;
for i = 100:-1:1
    f = find(spikeTimes == i);
    ff(i) = numel(f);
%    scatter(real(zDiff(f)),imag(zDiff(f)));
    zDm(i) = mean(zDiff(f));zs(i) = var(z(f));%.*conj(zDiff(f))
    zm(i) = mean(z(f));xm(i) = mean(xAll(f));
end
%plot(zDm,'k','LineWidth',2);

%%%%FIG 1
[vT vSS] = simVar(A,noise(end),modelSize,R);
[h xs] = histISI(spikes);
h = cumsum(h);h = h/max(h);
f = min(find(h >.95));
hIm(:,:,1) = ones(2,numel(h));
hIm(:,:,2) = ones(2,1)*h;
hIm(:,:,3) = hIm(:,:,2);
ys = [log(min(vT)) log(max(vT))];
figure;hold all;imagesc(xs,ys,hIm);
plot(log(covAllZ(spikeTimes,z,400)),'Linewidth',3);
plot(real(log(vT)),'Linewidth',3);%plot([1 500],ones(2,1)*log(vSS(1,1)),'k--');
plot([f f],ys+[-.1 .1],'k--','Linewidth',2);
axis tight;
legend({'Actual','AR3'},'Location','Southeast');legend boxoff;
set(gca,'fontsize',16);
ylabel('Log Variance');
xlabel('Time after spike');
set(gca,'xlim',[1 400]);
% figure;plot(real(z(size(spikeTimes,1),:)));hold on;plot(real(x(1,:)),'r');hold on;
% scatter(find(spikes(:,size(spikeTimes,1))),zeros(1,sum(spikes(:,size(spikeTimes,1)))),'k');

function [h xs] = histISI(spikes)
xs = 0:1:700;
h = zeros(size(xs));
for i = 1:size(spikes,2)
    f = find(spikes(:,i));
    f = diff(f);
    h = hist(f,xs) + h;
end

function [v P] = simVar(A,noise,p,R)
E = zeros(p); E(1,1) = noise;
P = R;
numIt = 2000;
v = zeros(numIt,1);
for i = 1:numIt
    P = A*P*A' + E;
    v(i) = P(1,1);
end
v = real(v);
vecP = (eye(p^2,p^2) - kron(A,A))\E(:);
P = reshape(vecP,p,p);
[log(max(v)) log(P(1,1))]

function z = covAllZ(spikes,zs,forw)
[r c] = find(spikes == 1);
r((c+forw-1) > size(spikes,2)) = [];
c((c+forw-1) > size(spikes,2)) = [];
z = zeros(numel(r),forw);
for i = 1:numel(r)
    z(i,:) = zs(r(i),c(i)+(1:forw)-1);
end
z = var(z);

function z = covSomeZ(spikes,zs,back)
[r c] = find(spikes == back);
c = c - back;
r(c < 0) = []; c(c < 0) = [];
z = zeros(numel(r),back);
for i = 1:numel(r)
    z(i,:) = zs(r(i),c(i)+(1:back));
end
z = var(z);

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
    counter = 1001;
    for j= 1:s(2)
        if spikes(i,j)
            counter = 0;
        else
            counter = counter + 1;
        end
        t(i,j) = counter;
    end
end
t(t > 1000) = 500;
t = t';

function dPhase = getDiffs(sig,gap)
angleA = angle(sig);
if min(size(sig)) == 1
    dPhase = circ_dist(angleA((gap+1):end),angleA(1:(end-gap)));
else
    dPhase = circ_dist(angleA(:,(gap+1):end),angleA(:,1:(end-gap)));
end