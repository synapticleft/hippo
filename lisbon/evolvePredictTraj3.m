function [stats] = evolvePredictTraj3(fn,lambda,sigma,inds,timePast,offSet,trialHist,diff,startAlign,finalChoice) %sessNorm,
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


scramble = randperm(size(allData,1));
allData = allData(scramble,:,:);
trainInds = 1:floor(size(allData,1)*.95);
testInds = floor(size(allData,1)*.95)+1:size(allData,1);
bootStrap = 100;
for i = bootStrap:-1:1
        scrambles(i,:) = randperm(numel(trainInds));
        %ys(i,:) = y(scramble);
end
scrambles = [1:numel(trainInds); scrambles];
for i = size(allData,2):-1:timePast+1
    if numel(inds) == 1 || exist('finalChoice','var')
        X = allData(trainInds,i-timePast+1:i,3:end-1);
        X = [squeeze(allData(trainInds,i,1:2)) X(:,:)];
    else
        X = allData(trainInds,i-timePast+1:i,1:end-1);
        X = X(:,:);
    end
    if trialHist
        X = [X allOutShift(:,1)];% circshift([answer; binAnswer; choice;correct]',[-1 0])];
    end
    y = (squeeze(allData(trainInds,i,end)))';%zscore(squeeze(allData(:,i,end)))';
    XX(i,:,:) = X'*X;
    Xy(i,:,:) = X'*y(scrambles)';
end

for i = timePast+1:size(allData,2)
    if numel(inds) == 1 || exist('finalChoice','var')
        X = allData(trainInds,i-timePast+1:i,3:end-1);
        X = [squeeze(allData(trainInds,i,1:2)) X(:,:)];
    else
        X = allData(trainInds,i-timePast+1:i,1:end-1);
        X = X(:,:);
    end
    if trialHist
        X = [X allOutShift(:,1)];% circshift([answer; binAnswer; choice;correct]',[-1 0])];
    end
    y = (squeeze(allData(trainInds,i,end)))';
    weights = exp(-((1:size(XX,1))-i).^2/sigma.^2);
    XXtemp = squeeze(sum(bsxfun(@times,XX(weights>.05,:,:),weights(weights >.05)'),1));
    Xytemp = squeeze(sum(bsxfun(@times,Xy(weights>.05,:,:),weights(weights>.05)'),1));
    %XXtemp = reshape(weights*XX(:,:),size(XX,2),[]);
    %Xytemp = reshape(weights*Xy(:,:),size(XX,2),[]);
    kernB = Xytemp'/(XXtemp + lambda*diag(diag(XXtemp)));
    yEst = squeeze(kernB)*X';
    y = zscore(y(scrambles),[],2);
    yEst = zscore(yEst,[],2);
    mseB = mean((y-yEst).^2,2);
    %mseB = mseB./mean(yEst.^2,2);
    ccB = mean(y.*yEst,2);%zscore(y(scrambles)).*zscore(yEst),2);
    stats.kern(i-timePast,:) = kernB(1,:);
    stats.cc(i-timePast) = ccB(1);
    stats.mse(i-timePast) = mseB(1);
    stats.kernL(i-timePast,:) = squeeze(prctile(kernB(2:end,:),2.5))';
    stats.kernH(i-timePast,:) = squeeze(prctile(kernB(2:end,:),97.5))';
    stats.ccL(i-timePast) = prctile(ccB(2:end),2.5);
    stats.ccH(i-timePast) = prctile(ccB(2:end),97.5);
    stats.mseL(i-timePast) = prctile(mseB(2:end),2.5);
    stats.mseH(i-timePast) = prctile(mseB(2:end),97.5);
end


% %    X = (X);
% %    for k = 1:numel(lambda)
% %        [cc(k,i-timePast),mse(k,i-timePast),coeff(k,i-timePast,:) stats(i)] = ridgeCross(y,X,1,lambda(k),100);
% %    end
%     if numel(inds) == 1 || exist('finalChoice','var')
%         X = allData(testInds,i-timePast+1:i,3:end-1);
%         X = [squeeze(allData(testInds,i,1:2)) X(:,:)];
%     else
%         X = allData(testInds,i-timePast+1:i,1:end-1);
%         X = X(:,:);
%     end
%     yHat(i,:) = X(:,:)*squeeze(coeff(k,i-timePast,:));
%     %if 0 % figsOn
%     %if (i-timePast == 90) || (i-timePast == 120)
%     %    figure;scatter(y,yHat,'filled');
%     %end
%     %end
% end
% y = squeeze(allData(testInds,:,end));
% yHat = yHat';
% coeff = squeeze(coeff);

if 1%figsOn
%coeff = reshape(coeff,size(coeff,1),timePast,[]);
figure;imagesc((stats.kern./(stats.kernH-stats.kernL)*2)',[-2 2]);%coeff(:,:)'./([stats.kernH]-[stats.kernL])*2,[-2 2]);
%figure;plot(coeff(:,:));
figure;plot(stats.mse);hold all;plot(stats.cc);plot([stats.ccL; stats.ccH]');%max(mse)-mse);hold all;plot(cc);ylim([0 1]);
%[~,m] = min(mse);
%figure;plot(squeeze(coeff(90,:,:)));
end