function [coeffs y yHat] = evolvePredictTraj2(fn,lambda,inds,timePast,offSet,trialHist,diff,startAlign,finalChoice) %sessNorm,
% fit LDA and regularized linear regression of start- and end- aligned
% trials, look at performance, and weights.

% choices: 1) re-normalize measures in each session to correct for drifts
% (sessNorm) NOW I DO THIS FOR EYE DATA BUT NOT OTHERS
% 2) regularization of fit (lambda)
% 3) fit choice, previous choice, or right answer
% 4) number of time steps in the past
% 5) which measurements to include
% 6) separate each difficulty level or combine them all

diff = [-diff diff];
if numel(inds) == 1 || exist('finalChoice','var')
    [allData allOut allOutShift] = preProcessRoberto(fn,inds,timePast,offSet,[],diff(:),startAlign,1);
else
    [allData allOut allOutShift] = preProcessRoberto(fn,inds,timePast,offSet,[],diff(:),startAlign);
end

xcAccum = 0;
xyAccum = 0;
%allData(:,:,end) = allData(randperm(size(allData,1)),:,end);
%figure;hold on;
for i = 1:size(allData,1)
    y = allData(i,1:end,end);
    inds = find(y(:) ~= y(1) & (1:numel(y(:)))' > timePast);%s
    %scatter([i i],[inds(1) inds(end)],'b');
    if ~exist('finalChoice','var')
        x = squeeze(allData(i,:,1:end-1));
    else
        x = squeeze(allData(i,:,3:end-1));
    end
    xx = zeros(timePast*size(x,2),numel(inds)+timePast);
    warning off all;
    for j = 1:size(x,2)
        %try
            xx((j-1)*timePast+(1:timePast),:) = toeplitz(x(inds(1)-timePast:inds(end),j),zeros(timePast,1))';
        %catch
            %size(xx((j-1)*timePast+(1:timePast),:))
            %size(toeplitz(x(inds(1)-timePast:inds(end),j),zeros(timePast,1))')
            %scatter([i i],[inds(1) inds(end)],'r');
        %[inds(1) inds(end)]
        %end
    end
    xx = xx(:,timePast+1:end);
    if exist('finalChoice','var')
        xx = [squeeze(allData(i,inds,1:2))'; xx];
    end
    xyAccum = xyAccum + y(inds)*xx';
    xcAccum = xcAccum + xx*xx';
end
coeffs = xyAccum/(xcAccum+lambda*diag(diag(xcAccum)));

for i = 1:size(allData,1)
    if ~exist('finalChoice','var')
        x = squeeze(allData(i,:,1:end-1));
    else
        x = squeeze(allData(i,:,3:end-1));
    end
    y = allData(i,1:end,end) ;
    inds = find(y(:) ~= y(1) & (1:numel(y(:)))' > timePast);
    xx = zeros(timePast*size(x,2),numel(inds)+timePast);
    for j = 1:size(x,2)
        xx((j-1)*timePast+(1:timePast),:) = toeplitz(x(inds(1)-timePast:inds(end),j),zeros(timePast,1))';
    end
    xx = xx(:,timePast+1:end);
    if exist('finalChoice','var')
        xx = [squeeze(allData(i,inds,1:2))'; xx];
    end
    yHat(i,1:size(xx,2)) = coeffs*xx;
    %xyAccum = xyAccum + y(inds)*xx';
    %xcAccum = xcAccum + xx*xx';
end
y = allData(:,1:end,end);
return

for i = timePast+1:size(allData,2)
    if numel(inds) == 1 || exist('finalChoice','var')
        X = [allData(:,i-timePast+1:i,3:end-1)];
        X = [squeeze(allData(:,i,1:2)) X(:,:)];
    else
        X = allData(:,i-timePast+1:i,1:end-1);
        X = X(:,:);
    end
    if trialHist
        X = [X allOutShift(:,1)];% circshift([answer; binAnswer; choice;correct]',[-1 0])];
    end
    y = (squeeze(allData(:,i,end)))';%zscore(squeeze(allData(:,i,end)))';
    X = (X);
    for k = 1:numel(lambda)
        [cc(k,i-timePast),mse(k,i-timePast),coeff(k,i-timePast,:) stats(i)] = ridgeCross(y,X,1,lambda(k),100);
    end
    yHat(i,:) = X*squeeze(coeff(k,i-timePast,:));
%     coeff(i-timePast,:) = (X'*y')'/(X'*X + lambda*size(X,1)*eye(size(X,2)));%y/X';
%     yHat = coeff(i-timePast,:)*X';
%     mse(i-timePast) = mean(sqrt((y-yHat).^2));%sum(yHat'+1 ~=info(inds,class))/sum(inds);
%     temp = corrcoef(y,yHat);
%     cc(i-timePast) = temp(1,2);
    %if 0 % figsOn
    %if (i-timePast == 90) || (i-timePast == 120)
    %    figure;scatter(y,yHat,'filled');
    %end
    %end
end
y = squeeze(allData(:,:,end));
yHat = yHat';