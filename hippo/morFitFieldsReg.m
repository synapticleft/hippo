function [allCorr ccs stds] = morFitFieldsReg(lfp,sp,pos,file)% allWs 
%% from filtered fits of spike trains, make trials x position activity maps

accumbins = 100;
dec = 32;
decs = 1;
levels = 1:.5:9;
lfp = single(lfp);
sp = single(sp);
[~,pos,isFast,runNum] = fixPos(pos);
isFast1 = logical(round(interp([double(isFast); 0],dec/decs)));
pos1 = interp(pos,dec/decs);
isFast1 = isFast1(1:size(lfp,2));
runNum1 = max(1,min(max(runNum),round(interp(runNum,dec/decs))));
pos1 = max(1,min(accumbins*2,round(pos1*accumbins)));
pos1 = pos1(isFast1);
runNum1 = runNum1(isFast1)';
lambdas = [0 4.^(-8:4)];
ind = [5 17 1:4 6:16];
    
for i = 1:numel(levels)
    y = morFilter(sp,2^levels(ind(i)),1250/decs);
    X = morFilter(lfp,2^levels(ind(i)),1250/decs);
    X = double(X(:,isFast1));y = double(y(:,isFast1));
    stds(ind(i),:) = std(X,0,2);
    X = bsxfun(@rdivide,X,stds(ind(i),:)');
    allCorr(ind(i),:,:) = y*X'/size(X,2);
    [~,ccs(1,ind(i),:,:)] = ridgeScan(y,X,[],lambdas);%allWs(1,i,:,:,:) 
    [~,ccs(2,ind(i),:,:)] = ridgeScan(y,X,pos1<=accumbins,lambdas);
    [~,ccs(3,ind(i),:,:)] = ridgeScan(y,X,mod(pos1,2) == 1,lambdas);
    [~,ccs(4,ind(i),:,:)] = ridgeScan(y,X,mod(runNum1,2) == 1,lambdas);
%     allRuns(1,i,:,:,:) = makeBins(y,pos1,runNum1);
%     allRuns(2,i,:,:,:) = makeBins(squeeze(allWs(2,i,:,:))*X,pos1,runNum1);
%     allRuns(3,i,:,:,:) = makeBins(squeeze(allWs(3,i,:,:))*X,pos1,runNum1);
%     allRuns(4,i,:,:,:) = makeBins(squeeze(allWs(4,i,:,:))*X,pos1,runNum1);
%     allRuns(5,i,:,:,:) = makeBins(squeeze(allWs(5,i,:,:))*X,pos1,runNum1);
    ind(i)
    save([file 'SpkFields.mat'],'ccs');
end

%save([file num2str(i) 'SpkFields.mat'],'allCorr','ccs');%,'allWs','allRuns','stds'

function [Ws,ccs] = ridgeScan(y,X,train,lambdas)%,pos,runNum),runs
if isempty(train)
    test = true(1,size(X,2));%1:size(X,2);
    train = test;
else
    test = xor(true(size(X,2),1),train);%setdiff(1:size(X,2),train);
end
yX = y(:,train)*X(:,train)';
XX = X(:,train)*X(:,train)';
for i = numel(lambdas):-1:1
    Ws(i,:,:) = yX/(XX + eye(size(XX,1))*sum(train)*lambdas(i));
    ccs(i,:) = myMSE(y(:,test),squeeze(Ws(i,:,:))*X(:,test));
end
yX = y(:,test)*X(:,test)';
XX = X(:,test)*X(:,test)';
for i = numel(lambdas):-1:1
    Ws1(i,:,:) = yX/(XX + eye(size(XX,1))*sum(test)*lambdas(i));
    ccs(i,:) = ccs(i,:) + myMSE(y(:,train),squeeze(Ws1(i,:,:))*X(:,train)).';
end
Ws = (Ws+Ws1)/2;
ccs = ccs/2;
%runs = makeBins(y,pos,runNum);

function mse = myMSE(y,yhat)
mse = sum(abs(y-yhat).^2,2)./sum(y.*conj(y),2);

function cc = myCorr(x,y)
cc = sum(x.*conj(y),2)./(sqrt(sum(x.*conj(x),2)).*sqrt(sum(y.*conj(y),2)));

function binData = makeBins(y,pos,runNum)
y = double(abs(y));
for i = size(y,1):-1:1
    binData(i,:,:) = accumarray([pos runNum'],y(i,:),[max(pos) max(runNum)],@mean);
end