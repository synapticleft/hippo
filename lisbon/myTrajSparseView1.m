function [XFull,AFull,reg,stimOr] = myTrajSparseView1(fn,phi,lambda)%dewhiteningMatrix,phi) whiteningMatrix
%% run convolutional sparse coding, then bin and render activations.
%% uses cwt to orthogonalize temporal dependencies.
load(fn,'data','center_ILD');
inds = [39 31];%31 [15 16];29 30 31 
indStim = 40;
choice = [data{2:end,7}];
f = choice ~= 3;
for i = 2:size(data,1)
    if ~numel(data{i,inds(1)})
        f(i-1) = 0;
    end
end
ord = 4;
X = [];reg = [];XFull = [];seeds = [];stimOr = [];
scales = 2.^(1:3);%(1:.5:2.5);%3.5
wname = 'haar';%'morl';%rbio3.1 and gaus1 are alternatives
%for i = 2:size(data,1)
    %if f(i-1)
    f = find(f);
    for i = 1:numel(f)
        trialData = reshape([data{f(i)+1,inds}],[numel(data{f(i)+1,indStim}) numel(inds)])';
        %fs = [find(~isnan(data{i,indStim}),1) find(~isnan(data{i,indStim}),1,'last')];
        trajcwt = [];
        for j =  1:size(trialData,1)
            %trajcwt = toeplitz(zeros(1,ord),trialData(j,:))];%[trajcwt; cwt(trialData(j,:),scales,wname)];
            trajcwt = [trajcwt; diff(trialData(j,:))];
        end
        %         seedTemp = 0;
        %         for j = 1:size(phi,2)
        %             for k = 1:size(phi,1)
        %                 seedTemp(j,:) = seedTemp(j,:) + conv(trajcwt(k,:),squeeze(phi(k,j,:)));
        %             end
        %         end
        %         seeds = [seeds seedTemp(:,~isna];
        trajcwt = trajcwt(:,~isnan(data{f(i)+1,indStim}));
        reg = [reg f(i)*ones(1,size(trajcwt,2))];
        X = [X trajcwt];
        stimDat = data{f(i)+1,indStim} - center_ILD;
        %stimcwt = cwt(stimDat(~isnan(data{i,6})),scales,wname);
        XFull = [XFull trialData(:,~isnan(data{f(i)+1,indStim}))];%(f)];
        %stim = [stim stimcwt];
        stimOr = [stimOr stimDat(~isnan(data{f(i)+1,indStim}))];
        %trialNum(i) = size(X,2);
        %    end
    end
%f = find(f);
%XFull = X;
X = zscore(X,0,2);
%X = whiten(X);
%X = whiten(X);
J = size(phi,2);		% number of basis functions for source generation
R = size(phi,3);%20;		% number of time points in basis functions generating sources
N = size(X,1);		% number of sources
randn('seed',1);
rand('seed',1);

Jrows = 4;
paramstr = sprintf('%s_J=%03d_R=%03d_N=%03d_%s', 'hippo', J, R, N, datestr(now,30));
update = 1;
if ~exist('phi','var') || isempty(phi)
phi = randn(N,J,R);
for j = 1:J
    phi(:,j,:) = phi(:,j,:) / sqrt(sum(sum(phi(:,j,:).^2)));
end
end
randn('seed',2);
rand('seed',2);
     %lambda = [.5 .15];%%[.1 .5]
     cols = colormap('jet');%[1 0 0; 0 1 0];
AFull = [];     
useInds = 1:J;
%f1 = figure;f2 = figure;
for t = 1:numel(f)
    %% select data
    Xsamp = X(:,reg == f(t));
    %% compute the map estimate
    S = size(Xsamp,2);
    P = S+R-1;	% number of selection locations for source generation
    a0 = zeros(J,P);%randn(J, P);
    opts.x0 = a0(:);
    callF   = @(x) objfun_a_conv(x,Xsamp,phi,lambda);
    a1 = lbfgsb(callF,Inf*ones(1,J*P)',Inf*ones(1,J*P)',opts);
    a1 = reshape(a1, J, P);a1 = a1(:,R/2+(1:S));
    AFull = [AFull a1];
%    figure(f1);imagesc(a1);
%      figure(f2);
%     for j = 1:J
%         subplot(4,ceil(J/4),j);
%         inds = reg == f(t);
%         temp = [XFull(end,inds); 1:sum(inds)];
%         ff = abs(a1(useInds(j),:)) > .1;
%         if sum(ff)
%             actI = max(1,min(64,32+round(20*a1(useInds(j),ff))));
%         %plot(temp(1,min(numel(ff),max(1,find(ff,1)+(-10:10)))),temp(2,min(numel(ff),max(1,find(ff,1)+(-10:10)))),'k--','linewidth',.5);hold all;
%         scatter(temp(1,ff),temp(2,ff),abs(a1(useInds(j),ff))*40,cols(actI,:),'filled');hold all;axis tight;%sign(a1(useInds(j),ff))/2+1.5
%         end
%     end
     drawnow;
end