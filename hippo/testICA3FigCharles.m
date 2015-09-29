function testICA3FigCharles()
%% This function runs ICA on simulated LFPs, trying a variety of settings
%% corresponding to the SNR of single neurons, and the smoothness of their place-fields
%% across space. All of the results are then saved in a file.
%% requires fastica (real-valued case) or ACMNsym (complex-valued/oscillatory case)

%% simulation for paper using real-valued, non-oscillatory activity
allSNR = linspace(0,.5,11);             %the range of SNR's to test
allSmooth = max(1,linspace(0,50,11));   %the range of receptive field smoothness to test
nSources = 1000;                        %number of neurons in simulation
trials = 50;                            %number of runs across virtual track
div = 100;                              % neurons are simulated in batches of this size due to memory constraints
pos = linspace(0,1,100);                %discretized virtual track
nChannels = 64;                         %number of channels in virtual electrode array
sigma = 2;                              %responses on the array are blurred by this factor to mimic LFP spreading
scaleSources = floor(nSources/div);     %number of batches of neurons for simulation
nSources = div;                         %number of neurons per batch
xcLen = round(numel(pos)/4);            %length of spatial cross-correlation (measure of RF smoothness)
xc = zeros(numel(allSNR),numel(allSmooth),2,xcLen*2+1); %spatial cross-correlation of neural RFs
useReal = 1;                            %simulate firing rates (rather than complex-valued oscillatory signals)
%% run simulation
for kk = numel(allSNR):-1:1
    for ll = numel(allSmooth):-1:1
        SNRs = allSNR(kk);smoothness = allSmooth(ll);
        gWinPos = gausswin(smoothness);
        Xf = zeros(prod(nChannels),trials*numel(pos));  %time series of simulated LFP
        mix = zeros(nSources,prod(nChannels));          %how neural activity is combined to produce LFP (mixing matrix)
        acts = zeros(nSources,trials*numel(pos));       %time series of spike trains
        for ii = 1:scaleSources
            tic;
            for i = 1:nSources
                temp = exp(-(linspace(1,prod(nChannels),prod(nChannels))-...
                    (rand*1.1-.05)*prod(nChannels)).^2/sigma(1).^2); %blurred electrical field of a single neuron
                mix(i,:) = temp(:);                                     %assigned to a row of the mixing matrix
                act = max(0,filter(gWinPos,1,randn(numel(pos)*2,1)))';  %place field of the neuron
                act = act(:,floor(numel(pos))+(1:numel(pos)));          %using the last half to remove initial filter transient
                act = max(0,act);                                       %rectify signal
                temp1 = (max(0,1+randn(1,trials)*SNRs(1))'*act)';       %multiple place field by noise in each run
                acts(i,:) = temp1(:);
            end
            actsM = squeeze(mean(reshape(acts,[size(acts,1) numel(pos) trials]),3)); %reshape the activations
            for i = 1:size(acts,1)                                                   %and calculate autocorrelation for each neuron
                xc(kk,ll,1,:) = squeeze(xc(kk,ll,1,:))' + xcov(actsM(i,:),xcLen);
            end
            Xf = Xf + mix'*acts;                                                    %add all activity to LFP
            toc
        end
        if useReal %real-valued analysis
            Xf = bsxfun(@minus,Xf,mean(Xf,2));
            [A W] = fastica(Xf,'approach','symm','g','tanh','maxNumIterations',500);
            sk = sign(skewness((W*(Xf))')');
            c = bsxfun(@times,W*(Xf),sk);
        else        %complex-valued analysis
            [A W] = ACMNsym(bsxfun(@times,Xf,exp(1i*(1:size(Xf,2)))),'mle_circ');
            c = W*Xf;
        end
        c = reshape(c,[size(c,1) numel(pos) trials]);
        c = permute(c,[1 3 2]);
        temp = c;temp(:,:) = zscore(temp(:,:),0,2);
        c = squeeze(mean(temp,2));  %average ICA feature activations across trials to derive position-dependent activity
        [mx] = max(abs(c),[],2);
        f = find(mx > 2);           %find ICA features that are place-sensitive
        fn(kk,ll) = numel(f);       %number of place sensitive features for each setting
        Xf1 = A(:,f)*(W(f,:)*Xf);
        varExp(kk,ll) = abs(corr(Xf1(:),Xf(:)))^2; %fraction of total variance explained by place-selective features
        for i = 1:numel(f) %spatial autocorrelation of each feature
            xc(kk,ll,2,:) = squeeze(xc(kk,ll,2,:))' + xcov(c(f(i),:),xcLen);
        end
        allCs{kk,ll} = c; %save all position-selective activations (for plotting individually)
    end
end
save('allSimC.mat','allCs','xc','varExp','fn');

%figure;
%    col = colormap(hsv(numel(pos)));
%    for i = 1:size(c,1)
%        [~,m] = max(c(i,:));
%        cc(i,:) = col(m,:);
%    end
%set(gca,'nextPlot','add','ColorOrder',cc(f,:));
%xs = repmat(pos,[numel(f) 1]);
%temps = std(temp(f,:,:),0,2);
%plot(xs',c(f,:)','linewidth',2);%axis tight;
%set(gca,'fontsize',16,'color',[0 0 0]);xlabel 'position'; ylabel 'activation';drawnow;
%xc = bsxfun(@rdivide,xc,max(abs(xc),[],2));
%if figsOn
%figure;
%plot(-xcLen:xcLen,xc(1,:),'b','linewidth',2);hold all;
%plot(-xcLen:xcLen,xc(2,:),'r','linewidth',2);
%legend({'Sources','Unsupervised','Supervised'});
%set(gca,'fontsize',16,'xtick',[],'ytick',[]);xlabel('distance');ylabel('autocorrelation');title('Tuning');
%end