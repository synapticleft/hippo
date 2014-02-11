function morFitFields(lfp,sp,pos,file)
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
runNum1 = max(1,min(max(runNum),round(interp(runNum,dec/decs))));
pos1 = max(1,min(accumbins*2,round(pos1*accumbins)));
pos1 = pos1(isFast1);
runNum1 = runNum1(isFast1);
    
for i = numel(levels):-1:1
    y = morFilter(sp,2^levels(i),1250/decs);
    X = morFilter(lfp,2^levels(i),1250/decs);
    X = double(X(:,isFast1));y = double(y(:,isFast1));
    allWs(1,i,:,:) = y*X';
    allWs(2,i,:,:) = y*X'/(X*X' + 0*eye(size(X,1)));%y/X;
    inds = pos1 <=accumbins;
    allWs(3,i,:,:) = y(:,inds)*X(:,inds)'/(X(:,inds)*X(:,inds)' + 0*eye(size(X,1)));
    inds = mod(pos1,2) == 1;
    allWs(4,i,:,:) = y(:,inds)*X(:,inds)'/(X(:,inds)*X(:,inds)' + 0*eye(size(X,1)));
    inds = mod(runNum1,2) == 1;
    allWs(5,i,:,:) = y(:,inds)*X(:,inds)'/(X(:,inds)*X(:,inds)' + 0*eye(size(X,1)));
    ccs(1,i,:) = myCorr(y,squeeze(allWs(2,i,:,:))*X);
    ccs(2,i,:) = myCorr(y(:,pos1 > accumbins),squeeze(allWs(3,i,:,:))*X(:,pos1 > accumbins));
    ccs(3,i,:) = myCorr(y(:,mod(pos1,2) == 0),squeeze(allWs(4,i,:,:))*X(:,mod(pos1,2) == 0));
    ccs(4,i,:) = myCorr(y(:,mod(runNum1,2) == 0),squeeze(allWs(5,i,:,:))*X(:,mod(runNum1,2) == 0));
    allRuns(1,i,:,:,:) = makeBins(y,pos1,runNum1);
    allRuns(2,i,:,:,:) = makeBins(squeeze(allWs(2,i,:,:))*X,pos1,runNum1);
    allRuns(3,i,:,:,:) = makeBins(squeeze(allWs(3,i,:,:))*X,pos1,runNum1);
    allRuns(4,i,:,:,:) = makeBins(squeeze(allWs(4,i,:,:))*X,pos1,runNum1);
    allRuns(5,i,:,:,:) = makeBins(squeeze(allWs(5,i,:,:))*X,pos1,runNum1);
%     for j = size(y,1):-1:1
%         absDec = double(abs(y(j,:)));
%         allRuns(1,i,j,:,:) = accumarray([pos1 runNum1'],absDec,[max(pos1) max(runNum)],@mean);
%         absDec = double(abs(squeeze(allWs(2,i,j,:)).'*X));
%         allRuns(2,i,j,:,:) = accumarray([pos1 runNum1'],absDec,[max(pos1) max(runNum)],@mean);
%         absDec = double(abs(squeeze(allWs(3,i,j,:)).'*X));
%         allRuns(3,i,j,:,:) = accumarray([pos1 runNum1'],absDec,[max(pos1) max(runNum)],@mean);
%         absDec = double(abs(squeeze(allWs(4,i,j,:)).'*X));
%         allRuns(4,i,j,:,:) = accumarray([pos1 runNum1'],absDec,[max(pos1) max(runNum)],@mean);
%         %showGrid(squeeze(allRuns(:,i,j,:,:)));drawnow;
%     end
    i
end

save([file 'SpkFields.mat'],'allRuns','allWs','ccs');


function cc = myCorr(x,y)
cc = sum(x.*conj(y),2)./(sqrt(sum(x.*conj(x),2)).*sqrt(sum(y.*conj(y),2)));

function binData = makeBins(y,pos,runNum)
y = double(abs(y));
for i = size(y,1):-1:1
    binData(i,:,:) = accumarray([pos runNum'],y(i,:),[max(pos) max(runNum)],@mean);
end