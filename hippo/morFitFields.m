function morFitFields(lfp,sp,pos,file)
%% from filtered fits of spike trains, make trials x position activity maps

accumbins = 100;
dec = 32;
decs = 4;
levels = 2:.5:7;

[~,pos,isFast,runNum] = fixPos(pos);
isFast1 = logical(round(interp([double(isFast); 0],dec/decs)));
pos1 = interp(pos,dec/decs);
runNum1 = max(1,min(max(runNum),round(interp(runNum,dec/decs))));
pos1 = max(1,min(accumbins*2,round(pos1*accumbins)));
pos1 = pos1(isFast1);
    
for i = numel(levels):-1:1
    y = morFilter(sp,2^levels(i),1250/decs);
    X = morFilter(lfp,2^levels(i),1250/decs);
    X = X(:,isFast1);y = y(:,isFast1);
    allWs(4,i,:,:) = y*X';
    allWs(3,i,:,:) = y*X'/(X*X' + 0*eye(size(X,1)));%y/X;
    inds = pos1 <=accumbins;
    allWs(2,i,:,:) = y(:,inds)*X(:,inds)'/(X(:,inds)*X(:,inds)' + 0*eye(size(X,1)));
    inds = mod(pos1,2) == 1;
    allWs(1,i,:,:) = y(:,inds)*X(:,inds)'/(X(:,inds)*X(:,inds)' + 0*eye(size(X,1)));
    for j = size(y,1):-1:1
        ccs(3,i,j) = corr(y(j,:)',(squeeze(allWs(3,i,j,:)).'*X)');
        ccs(2,i,j) = corr(y(j,pos1 > accumbins)',(squeeze(allWs(2,i,j,:)).'*X(:,pos1 > accumbins))');
        ccs(1,i,j) = corr(y(j,mod(pos1,2) == 0)',(squeeze(allWs(1,i,j,:)).'*X(:,mod(pos1,2) == 0))');
        absDec = double(abs(squeeze(allWs(3,i,j,:)).'*X));
        allRuns(4,i,j,:,:) = accumarray([pos1 runNum1(isFast1)'],absDec,[max(pos1) max(runNum)],@mean);
        absDec = double(abs(squeeze(allWs(2,i,j,:)).'*X));
        allRuns(3,i,j,:,:) = accumarray([pos1 runNum1(isFast1)'],absDec,[max(pos1) max(runNum)],@mean);
        absDec = double(abs(squeeze(allWs(1,i,j,:)).'*X));
        allRuns(2,i,j,:,:) = accumarray([pos1 runNum1(isFast1)'],absDec,[max(pos1) max(runNum)],@mean);
        absDec = double(abs(y(j,:)));
        allRuns(1,i,j,:,:) = accumarray([pos1 runNum1(isFast1)'],absDec,[max(pos1) max(runNum)],@mean);
        %showGrid(squeeze(allRuns(:,i,j,:,:)));drawnow;
    end
    i
end

save([file 'SpkFields.mat'],'allRuns','allWs','ccs');

%%%%%%%%%%%%%%
%d = '/media/work/hippocampus/ec013.670/';
%m = memmapfile([[d file] 'CwtSpk.dat']);
%dims = [30 1560096];
%dims(3) = numel(m.Data)/prod(dims)/4;

%m = memmapfile([[d file] 'CwtSpk.dat'],'Format',{'single' dims 'X'});


%allRuns = zeros(numel(levels),dims(3),max(pos),max(runNum));
%runNum = round(interp(runNum,dec));
%isFast = logical(round(interp([double(isFast); 0],dec)));
%pos = interp(pos,dec);

%for i = 1:numel(levels)
%    X = permute(m.Data(1).X([i i+numel(levels)],:,:),[2 1 3]);
%    for j = 1:dims(3)
%        absDec = decimate(squeeze(abs(double(complex(X(:,1,j),X(:,2,j))))),dec);
%        allRuns(i,j,:,:) = accumarray([pos(isFast) runNum(isFast)'],absDec(isFast),[max(pos) max(runNum)],@mean);
%    end
%    i
%end

%save([file 'SpkFields.mat'],'allRuns');