function [allCorr allWs ccs stds] = morFitFieldsReg(lfp,sp,pos,file)
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
    
for i = numel(levels):-1:1
    y = morFilter(sp,2^levels(i),1250/decs);
    X = morFilter(lfp,2^levels(i),1250/decs);
    X = double(X(:,isFast1));y = double(y(:,isFast1));
    stds(i,:) = std(X,0,2);
    X = bsxfun(@rdivide,X,stds(i,:)');
    allCorr(i,:,:) = y*X'/size(X,2);
    [allWs(1,i,:,:,:) ccs(1,i,:,:)] = ridgeScan(y,X,[],lambdas);
    [allWs(2,i,:,:,:) ccs(2,i,:,:)] = ridgeScan(y,X,pos1<=accumbins,lambdas);
    [allWs(3,i,:,:,:) ccs(3,i,:,:)] = ridgeScan(y,X,mod(pos1,2) == 1,lambdas);
    [allWs(4,i,:,:,:) ccs(4,i,:,:)] = ridgeScan(y,X,mod(runNum1,2) == 1,lambdas);
    %allWs(2,i,:,:) = y*X'/(X*X' + 0*eye(size(X,1)));%y/X;
    %inds = pos1 <=accumbins;
    %allWs(3,i,:,:) = y(:,inds)*X(:,inds)'/(X(:,inds)*X(:,inds)' + 0*eye(size(X,1)));
    %inds = mod(pos1,2) == 1;
    %allWs(4,i,:,:) = y(:,inds)*X(:,inds)'/(X(:,inds)*X(:,inds)' + 0*eye(size(X,1)));
    %inds = mod(runNum1,2) == 1;
    %allWs(5,i,:,:) = y(:,inds)*X(:,inds)'/(X(:,inds)*X(:,inds)' + 0*eye(size(X,1)));
    %ccs(1,i,:) = myCorr(y,squeeze(allWs(2,i,:,:))*X);
    %ccs(2,i,:) = myCorr(y(:,pos1 > accumbins),squeeze(allWs(3,i,:,:))*X(:,pos1 > accumbins));
    %ccs(3,i,:) = myCorr(y(:,mod(pos1,2) == 0),squeeze(allWs(4,i,:,:))*X(:,mod(pos1,2) == 0));
    %ccs(4,i,:) = myCorr(y(:,mod(runNum1,2) == 0),squeeze(allWs(5,i,:,:))*X(:,mod(runNum1,2) == 0));
%     allRuns(1,i,:,:,:) = makeBins(y,pos1,runNum1);
%     allRuns(2,i,:,:,:) = makeBins(squeeze(allWs(2,i,:,:))*X,pos1,runNum1);
%     allRuns(3,i,:,:,:) = makeBins(squeeze(allWs(3,i,:,:))*X,pos1,runNum1);
%     allRuns(4,i,:,:,:) = makeBins(squeeze(allWs(4,i,:,:))*X,pos1,runNum1);
%     allRuns(5,i,:,:,:) = makeBins(squeeze(allWs(5,i,:,:))*X,pos1,runNum1);
    i
end

%save([file 'SpkFields.mat'],'allCorr','allRuns','allWs','ccs','stds');

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
    ccs(i,:) = myCorr(y(:,test),squeeze(Ws(i,:,:))*X(:,test));
end
yX = y(:,test)*X(:,test)';
XX = X(:,test)*X(:,test)';
for i = numel(lambdas):-1:1
    Ws1(i,:,:) = yX/(XX + eye(size(XX,1))*sum(test)*lambdas(i));
    ccs(i,:) = ccs(i,:) + myCorr(y(:,train),squeeze(Ws1(i,:,:))*X(:,train)).';
end
Ws = (Ws+Ws1)/2;
ccs = ccs/2;
%runs = makeBins(y,pos,runNum);

function cc = myCorr(x,y)
cc = sum(x.*conj(y),2)./(sqrt(sum(x.*conj(x),2)).*sqrt(sum(y.*conj(y),2)));

function binData = makeBins(y,pos,runNum)
y = double(abs(y));
for i = size(y,1):-1:1
    binData(i,:,:) = accumarray([pos runNum'],y(i,:),[max(pos) max(runNum)],@mean);
end