function Xf = testICA(pos,v,shanks,sourceSNR) %[W,t,tacts]
%%old version with multiple fields per neuron

sz = min(size(pos,1),numel(v));
accumbins = 100;
pos = pos(1:sz,:);v = v(1:sz);
thresh = .05;
nChannels = [8 shanks];
sigma = [.5 1]*5;
nSources = 100;
scaleSources = 1;
nFields = poissrnd(2,nSources);
%sourceSNR = inf;

[xs ys] = meshgrid(1:nChannels(2),1:nChannels(1));
A = zeros(nSources,prod(nChannels));
[pos,posd] = fixPos(pos);
vel = angVel(pos);
vel = [zeros(1,size(vel,2)); vel];
vel = filtLow(vel(:,1),1250/32,1);
vel = vel/max(vel);inds = vel > thresh;
Xf = zeros(prod(nChannels),numel(v));
v = v(inds);posd = posd(inds);
posd(posd > 1) = posd(posd > 1) + 1; 
for ii = 1:scaleSources
    tic
    acts = zeros(nSources,numel(v));
    for i = 1:nSources
        temp = exp(-((xs-rand*nChannels(1)).^2/sigma(1).^2 + (ys-rand*nChannels(2)).^2/sigma(2).^2));
        %imagesc(temp,[0 max(temp(:))]);drawnow;
        A(i,:) = temp(:);
        for k = 1:nFields(i)
            r = rand*2;
            if r > 1
                r = r+1;
                sn = 1;
            else
                sn = -1;
            end
            acts(i,:) = acts(i,:) + exprnd(1)*giveCmor(posd-r,sn*5,rand/200).';%exprnd(.005));
        end
    end
    acts = acts + complex(randn(size(acts)),randn(size(acts)))*std(v)/sourceSNR;
    toc
    sPlot(acts(:,10000:15000),[],0);
    Xf(:,inds) = Xf(:,inds) + A'*acts;
end
Xf(:,inds) = bsxfun(@times,Xf(:,inds),v');
% posd(posd > 1) = posd(posd > 1) -1;
% tacts = zeros(size(acts,1),accumbins*2);
% for i = 1:size(acts,1)
%     tacts(i,:) = accumarray([ones(numel(posd),1) 1+floor(posd*accumbins)],acts(i,:),[1 accumbins*2],@mean);
% end
% 
% sPlot(Xf,[],0);
% [A W Z alphas] = ACMNsym(Xf,'mle_circ');
% acts = W*bsxfun(@times,Xf,exp(1i*angle(v.')));
% t = zeros(size(acts,1),accumbins*3);
% for i = 1:size(acts,1)
%     t(i,:) = accumarray([ones(numel(posd),1) 1+floor(posd*accumbins)],acts(i,:),[1 accumbins*2],@mean);
% end

function mor = giveCmor(pos,Fc,Fb)
mor = exp(2*1i*pi*Fc*pos).*exp(-(pos.*pos)/Fb);%((pi*Fb)^(-0.5))*