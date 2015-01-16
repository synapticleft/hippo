function [coeff,cc,mse,stats,y,yHat] = evolvePredictTraj1(fn,lambda,inds,timePast,offSet,trialHist,diff,startAlign,finalChoice) %sessNorm,
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
    y = (squeeze(allData(trainInds,i,end)))';%zscore(squeeze(allData(:,i,end)))';
    X = (X);
    for k = 1:numel(lambda)
        [cc(k,i-timePast),mse(k,i-timePast),coeff(k,i-timePast,:) stats(i)] = ridgeCross(y,X,1,lambda(k),100);
    end
    if numel(inds) == 1 || exist('finalChoice','var')
        X = allData(testInds,i-timePast+1:i,3:end-1);
        X = [squeeze(allData(testInds,i,1:2)) X(:,:)];
    else
        X = allData(testInds,i-timePast+1:i,1:end-1);
        X = X(:,:);
    end
    yHat(i,:) = X(:,:)*squeeze(coeff(k,i-timePast,:));
    %if 0 % figsOn
    %if (i-timePast == 90) || (i-timePast == 120)
    %    figure;scatter(y,yHat,'filled');
    %end
    %end
end
y = squeeze(allData(testInds,:,end));
yHat = yHat';
coeff = squeeze(coeff);

if 1%figsOn
%coeff = reshape(coeff,size(coeff,1),timePast,[]);
figure;imagesc(coeff(:,:)'./([stats.kernH]-[stats.kernL])*2,[-2 2]);
%figure;plot(coeff(:,:));
figure;plot(max(mse)-mse);hold all;plot(cc);ylim([0 1]);plot([stats.ccL]);plot([stats.ccH]);
[~,m] = min(mse);
%figure;plot(squeeze(coeff(90,:,:)));
end