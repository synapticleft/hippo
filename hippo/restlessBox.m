function [allAct posMean] = restlessBox(act,pos,bounds,frameCol)
%% scatter activations in waiting box, with colors corresponding to the binned activation renderings

thresh = 3;
dec = round(size(pos,1)/size(act,2));
for i =1:2
    posd(:,i) = decimate(pos(:,i),dec);
end
posd = posd(1:size(act,2),:);
bounded = posd(:,1) > bounds(1) & posd(:,1) < bounds(2) & ...
    posd(:,2) > bounds(3) & posd(:,2) < bounds(4);
numExp = 13*8/dec;
bounded = bwmorph(bounded,'dilate',numExp);
for i =1:numExp
    bounded = bounded - bwmorph(bounded,'endpoints');
end
bounded = logical(bounded);
%figure;scatter(posd(:,1),posd(:,2),'filled');
%hold all;scatter(posd(bounded,1),posd(bounded,2),'filled','r');
bounded = bwlabel(bounded);
%xdim = ceil(sqrt(max(bounded)));
%ydim = ceil(max(bounded)/xdim);
%col = [0 0 1; 0 1 0; 1 0 0; 1 1 0; 0 1 1;1 0 1];
%f1 = figure;
%f2 = figure;
%allAct = zeros(size(act,1),max(bounded));
posMean = sideSeq(bounded,posd,mean(bounds(3:4)));
%posMean = posMean > 0;
for i = 1:2
    figure;hold all;%subplot(2,1,i);
    inds = bounded & posMean == i-1;
    %tempMag = abs(act(:,inds));
    %allAct(:,i) = sum(max(tempMag - thresh,0),2);
    for j = 1:size(act,1)
    [ys xs] = meshgrid(find(inds),j);%1:size(act,1));
    tempMag = abs(act(j,inds));
    scatter(posd(ys(tempMag > thresh),2),posd(ys(tempMag > thresh),1),...
        tempMag(tempMag > thresh)*3,getCol(xs,tempMag>thresh,frameCol),'filled');
    drawnow;
    end
    set(gca,'xlim',max(0,bounds(3:4)),'ylim',bounds(1:2),'color',[0 0 0]);
    axis square
end
return
for i = 1:max(bounded)
    tempMag = abs(act(:,bounded == i));
    allAct(:,i) = sum(max(tempMag - thresh,0),2);
    [ys xs] = meshgrid(1:sum(bounded == i),1:size(act,1));
    offSet = find(bounded == i,1);
    %scatter(posd(offSet+ys(tempMag > thresh),1),posd(offSet+ys(tempMag > thresh),2),...
    %    tempMag(tempMag > thresh)*5,col(mod(xs(tempMag > thresh),size(col,1))+1,:),'filled');
    figure(f2);hold all;
    %subplot(1,2,1+(posMean(i)<0));hold all;
    if exist('down','var')
        inds = zeros(size(act,1),1);
        inds(down) = 1;
        tempMag1 = tempMag;tempMag1(~inds,:) = 0;tempMag(logical(inds),:) = 0;    
        scatter(posd(offSet+ys(tempMag1 > thresh),1),posd(offSet+ys(tempMag1 > thresh),2),...
        tempMag1(tempMag1 > thresh)*10,getCol(xs,tempMag1,thresh,frameCol),'v','filled');
    end
    if exist('up','var')
        inds = zeros(size(act,1),1);
        inds(up) = 1;
        tempMag1 = tempMag;tempMag1(~inds,:) = 0;tempMag(logical(inds),:) = 0;    
        scatter(posd(offSet+ys(tempMag1 > thresh),1),posd(offSet+ys(tempMag1 > thresh),2),...
        tempMag1(tempMag1 > thresh)*10,getCol(xs,tempMag1,thresh,frameCol),'^','filled');
    end
    scatter(posd(offSet+ys(tempMag > thresh),1),posd(offSet+ys(tempMag > thresh),2),...
        tempMag(tempMag > thresh)*10,getCol(xs,tempMag,thresh,frameCol),'s','filled');
    %figure(f1);
    %subplot(1,2,1+(posMean(i)<0));hold all;
    %plot(posd(offSet+(0:5),1),posd(offSet+(0:5),2),'k--');
    %scatter(posd(offSet+ys(tempMag > thresh),1),posd(offSet+ys(tempMag > thresh),2),...
    %    tempMag(tempMag > thresh)*5,col(mod(xs(tempMag > thresh),size(col,1))+1,:),'filled');
    %%scatter(posd(bounded == i,1),posd(bounded == i,2),20,cc(ceil(linspace(.01,64,size(tempMag,2))),:),'filled');
    drawnow;
end
%for i = 1:2
%    subplot(1,2,i);set(gca,'xlim',max(0,bounds(1:2)),'ylim',bounds(3:4));
%end
figure(f2);set(gca,'xlim',max(0,bounds(1:2)),'ylim',bounds(3:4),'color',[0 0 0]);

function colVec = getCol(xs,inds,frameCol)    
%if ~exist('frameCol','var') || isempty(frameCol)
%        colVec = col(mod(xs(tempMag > thresh),size(col,1))+1,:);
%    else
        colVec = frameCol(xs(inds),:);
%    end

function [runDir posSeq] = sideSeq(runInds,pos,ref)
posSeq = zeros(1,max(runInds)+1);
runDir = zeros(size(runInds));
for i = 1:max(runInds)+1
    if i == 1
        inds = 1:find(runInds == 1,1)-1;
    elseif i == max(runInds)+1
        inds = find(runInds == i-1, 1, 'last' )+1:size(pos,1);
    else
        inds = (find(runInds == i-1)+1):(find(runInds == i)-1);
    end
    temp = [max(pos(inds,2))-ref min(pos(inds,2))-ref];
    if abs(temp(1)) > abs(temp(2))
        posSeq(i) = temp(1);
    else
        posSeq(i) = temp(2);
    end
    runDir(runInds == i) = posSeq(i) > 0;
end