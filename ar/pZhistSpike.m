function pZhistSpike(spikes,z)

p = 4;
spikes = logical(spikes);
[coeff noise] = arburg(z(:,1),p);
%thisNoise = makeComplex(numel(A),sqrt(noise(i)),[1 1]);
histZ = zeros(p+1,sum(spikes(:)));
for i =0:p
    histZ(i+1,:) = z(find(spikes)-i);
end
[a b c] = svd(histZ','econ');
angle(c(:,1))
figure;hold all;
for i = 1:p
    %plot(hist(log(abs(histZ(i+1,:))./abs(histZ(1,:))),-.05:.01:.05));
    plot(hist(angle(histZ(i+1,:)) - angle(histZ(1,:)),-pi:.01:pi));
end
%figure;plot(histZ);
A = [-coeff(2:end);eye(numel(coeff)-1)];A(end,:) = [];
P = zeros(numel(coeff)-1);
Q = P; Q(1,1) = noise(end);
for i = 1:p
    P = A*P*A' + Q;
    estNoise(i) = real(P(1,1));
end
[~,poles] = residue([1 zeros(1,p-1)],coeff);
poles
poles = mean(poles);
bounds{1} = linspace(-.1,.1,40);bounds{2} = linspace(-.1,.1,40);
seed = histZ(1,:);
figure;
for i = 1:p
    subplot(2,2,i);
    seed = seed*poles';
    temp = seed - histZ(i+1,:);
    [h c] = hist3([real(temp)' imag(temp)'],bounds);
    imagesc(c{2},c{1},h);%%scatter(real(temp),imag(temp));
end
figure;
for i = 1:p
    subplot(2,2,i);
    temp = makeComplex(sum(spikes(:)),sqrt(estNoise(i)),[1 1]);
    [h c] = hist3([real(temp) imag(temp)],bounds);
    imagesc(c{2},c{1},h);%%scatter(real(temp),imag(temp));
end
    

function c = makeComplex(len,dev,ratio)
tot = sum(ratio.^2);
c = dev/sqrt(tot)*complex(randn(len,1)*ratio(1),randn(len,1)*ratio(2));
