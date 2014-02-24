function [A W t] = testICA2d(pos,SNRs,v)
%% simulation for paper using real-valued, non-oscillatory activity
%allSNR = linspace(0,.5,11);
%allSmooth = max(1,linspace(0,50,11));
pos = pos(:,1:2);
nanInds = any(isnan(pos'));
pos(nanInds,:) = [];
if exist('v','var')
    v(nanInds) = [];
end
nSources = 1000;
%trials = 50;
div = 100;
%pos = linspace(0,1,100);
nChannels = 64;
sigma = 2;
scaleSources = floor(nSources/div);
nSources = div;
pos = bsxfun(@minus,pos,min(pos)) + eps;
pos = ceil(pos);
pos1 = (pos(:,1)-1)*max(pos(:,2)) + pos(:,2);
%xcLen = round(numel(pos)/4);
%xc = zeros(numel(allSNR),numel(allSmooth),2,xcLen*2+1);
%for kk = 1%numel(allSNR):-1:1
%    for ll = 4%numel(allSmooth):-1:1
%        SNRs = allSNR(kk);smoothness = allSmooth(ll);
        %SNRs = 0;
        smoothness = 20;
        %gWinPos = gausswin(smoothness);
        Xf = zeros(prod(nChannels),size(pos,1));
        mix = zeros(nSources,prod(nChannels));
        acts = zeros(nSources,size(pos,1));
        for ii = 1:scaleSources
            tic;
            act = max(0,imfilter(randn([max(pos) nSources]),fspecial('gaussian',smoothness,smoothness/2),'symmetric'));
            act = permute(act,[3 1 2]);act = act(:,:);acts = act(:,pos1);
            noise = filtLow(complex(randn(size(acts)),randn(size(acts))),1250/4,2);
            acts = acts.*max(0,1+noise/std(noise(:))*SNRs);
            mix = randn(nSources,prod(nChannels));
%             for i = 1:nSources
%                temp = exp(-(linspace(1,prod(nChannels),prod(nChannels))-(rand*1.1-.05)*prod(nChannels)).^2/sigma(1).^2);
%                mix(i,:) = temp(:);
%                %acts(i,:) = act(pos(:,1),pos(:,2),i);
%             end
                %act = max(0,imfilter(randn(max(pos)),fspecial('gaussian',smoothness,smoothness/2),'symmetric');%gWinPos,1,randn(numel(pos)*2,1)))';
                %act = act(:,floor(numel(pos))+(1:numel(pos)));
            %    act = max(0,act);
            %    temp1 = (max(0,1+randn(1,trials)*SNRs(1))'*act)';
            %    acts(i,:) = temp1(:);
            %end
            %actsM = squeeze(mean(reshape(acts,[size(acts,1) numel(pos) trials]),3));
            %for i = 1:size(acts,1)
            %    xc(kk,ll,1,:) = squeeze(xc(kk,ll,1,:))' + xcov(actsM(i,:),xcLen);
            %end
            Xf = Xf + mix'*acts;
            toc
        end
        if exist('v','var')
            Xf = bsxfun(@times,Xf,v');
            [A,W] = ACMNsym(Xf,'mle_circ');
        else
            Xf = bsxfun(@minus,Xf,mean(Xf,2));
            [A W] = fastica(Xf,'approach','symm','g','tanh','maxNumIterations',200);
            v = ones(size(Xf,2),1);
        end
        t = runTriggerView2d(pos,Xf,2,v,W);
        if ~isreal(W)
             t(:,:) = bsxfun(@times,t(:,:),exp(1i*angle(mean(A))).');
        else
             sk = sign(skewness((W*(Xf))')');
             t(:,:) = bsxfun(@times,t(:,:),sk);
        end
%         c = bsxfun(@times,W*(Xf),sk);
%         c = reshape(c,[size(c,1) numel(pos) trials]);%b*Xf;
%         c = permute(c,[1 3 2]);
%         temp = c;temp(:,:) = zscore(temp(:,:),0,2);%temp(:) = zscore(temp(:));
%         c = squeeze(mean(temp,2));
%         [mx] = max(c,[],2);
%         f = find(mx > 2);
%         fn(kk,ll) = numel(f);
%         Xf1 = A(:,f)*(W(f,:)*Xf);
%         varExp(kk,ll) = corr(Xf1(:),Xf(:))^2;
%         for i = 1:numel(f)
%             xc(kk,ll,2,:) = squeeze(xc(kk,ll,2,:))' + xcov(c(f(i),:),xcLen);
%         end
%        allCs{kk,ll} = c;
%    end
%end
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