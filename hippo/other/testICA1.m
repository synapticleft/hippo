function [A,Xf] = testICA1(shanks,SNRs,freq,nSources,A) %[W,t,tacts]
%% older version of simulation
div = 1000;
    pos = [linspace(0,1,100)];% linspace(0,1,2000) linspace(0,1,4000)];
    if ~freq
        v = ones(1,numel(pos));%real(exp(1i*linspace(0,numel(pos),numel(pos))));%
    else
        v = exp(1i*linspace(0,numel(pos),numel(pos)));    
    end
    nChannels = shanks;
if ~exist('Xf','var')
    sigma = [.5 1]*3;
    %nSources = 2;%100;
    if nSources > div
        scaleSources = floor(nSources/div);
        nSources = div;
    else
        scaleSources = 1;
    end
    %nFields = poissrnd(1.5,nSources);
    [xs ys] = meshgrid(1:nChannels(2),1:nChannels(1));
    mix = zeros(nSources,prod(nChannels));
    Xf = zeros(prod(nChannels),numel(v));
    for ii = 1:scaleSources
        tic
        acts = zeros(nSources,numel(v));
        for i = 1:nSources
            temp = exp(-((xs-rand*(nChannels(2)+1)).^2/sigma(1).^2 + (ys-rand*(nChannels(1)+1)).^2/sigma(2).^2));
            mix(i,:) = temp(:);
%            acts(i,:) = complex(randn(1,numel(v)),randn(1,numel(v)));
            %for k = 1:nFields(i)
                r = rand;%/5+.4;%randn/5+.5;
                acts(i,:) = acts(i,:) + exprnd(1)*giveCmor(pos'-r,freq,max(rand^2/100,.0001)).';%exprnd(1);
            %end
        end
%        acts = filtLow(acts,numel(pos),numel(pos)/50);
        acts = acts + complex(randn(size(acts)),randn(size(acts)))*std(v)/SNRs(1);%sqrt(mean(v*v'))
        if ~freq
            acts = real(acts);
        end
        toc
        %sPlot(acts,[],0);
        Xf = Xf + mix'*acts;
    end
    Xf = bsxfun(@times,Xf,v);
end
%Xf = bsxfun(@minus,Xf,mean(Xf,2));
if ~freq
    [a,b,c] = runica(Xf);
    figure;subplot(411);plot(abs(acts)');
    subplot(412);plot(abs(whiten(acts)'));
    subplot(413);plot(Xf');
    subplot(414);plot(c');
else
Xf = Xf + complex(randn(size(Xf)),randn(size(Xf)))*std(Xf(:))/SNRs(2);
if ~exist('A','var')
    [A,~,~] = cfastica(Xf);%ACMNSym(Xf,'mle_circ');%
end
Am = mean(A);
A = bsxfun(@times,A,exp(1i*-angle(Am)));
W = pinv(A);
A1 = bsxfun(@rdivide,A,sqrt(sum(A.*conj(A))));
A1 = bsxfun(@minus,A1,mean(A1.').');
%figure;showGrid(A1,nChannels);
temp = W*bsxfun(@times,Xf,conj(v));
thresh = 3;
tInds = find(sum(abs(temp) > thresh,2) > numel(v)/100);
figure;subplot(311);plot(abs(acts)');
subplot(312);plot(abs(whiten(Xf)'));axis tight;
subplot(313);plot(abs(temp(tInds,:)'));axis tight;%hold all;plot(abs(xcorr
%subplot(414);hold all;
%for i = 1:numel(tInds)
%   absInds = abs(temp(tInds(i),:)) > thresh / 2;
%   scatter(find(absInds),angle(temp(tInds(i),absInds)),abs(temp(tInds(i),absInds))*20,'filled');
%end
%params.Fs = 1000;params.tapers = [3 5];
%[S,f] = mtspectrumc(abs(temp(tInds,:))',params);
%plot(f,S);
%figure;xc = zeros(size(Xf,1),201);
%Xf = whiten(Xf);
%for i = 1:size(Xf,1)
%    xc(i,:) = xcorr(Xf(i,:),100,'coeff');
%end
%plot(abs(xc)');
end

function mor = giveCmor(pos,Fc,Fb)
mor = exp(2*1i*pi*Fc*pos).*exp(-(pos.*pos)/Fb);%((pi*Fb)^(-0.5))*