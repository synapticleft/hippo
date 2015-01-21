function [yOrig1 yFits1 allData kern] = evolvePredictTraj4(fn,lambda,inds,timePast,offSet,trialHist,diff,startAlign,finalChoice) %sessNorm,
% variant of evolvePredictTraj1 that incorporates smoothness of kernels
% over time.

% choices: 1) re-normalize measures in each session to correct for drifts
% (sessNorm) NOW I DO THIS FOR EYE DATA BUT NOT OTHERS
% 2) regularization of fit (lambda)
% 3) fit choice, previous choice, or right answer
% 4) number of time steps in the past
% 5) which measurements to include
% 6) separate each difficulty level or combine them all
%frac = 1;%.95
diff = [-diff diff];
if numel(inds) == 1 || (exist('finalChoice','var') && finalChoice == 1)
    [allData allOut allOutShift cILD] = preProcessRoberto(fn,inds,timePast,offSet,[],diff(:),startAlign,1);
else
    [allData allOut allOutShift cILD] = preProcessRoberto(fn,inds,timePast,offSet,[],diff(:),startAlign);
end
inds = groupILDs(allData,cILD,1);
yOrig1 = [];yFits1 = [];
% if frac == 1
%    scramble = 1:size(allData,1);
% else
%     scramble = randperm(size(allData,1));
% end
% allData = allData(scramble,:,:);
% for i = 1:size(allData,2)
%     f = find(squeeze(allData(i,:,end)) ~= 0);
%     allData(i,f(1):f(end),end) = filtfilt(gausswin(5),sum(gausswin(5)),squeeze(allData(i,f(1):f(end),end)));
% end
figure;
for j = 1:max(inds)
    trainInds = j ~= inds;%logical(ones(1,numel(inds)));%
    testInds = j == inds;
    yOrig = zeros(sum(testInds),size(allData,2)-timePast);
    yFits = yOrig;kern = [];
    for i = size(allData,2):-1:timePast+1
        if numel(inds) == 1 || exist('finalChoice','var')
            X = allData(:,i-timePast+1:i,3:end-1);
            X = [squeeze(allData(:,i,1:2)) X(:,:)];
        else
            X = allData(:,i-timePast+1:i,1:end-1);
            X = X(:,:);
        end
        if trialHist
            X = [X allOutShift(:,1)];% circshift([answer; binAnswer; choice;correct]',[-1 0])];
        end
        y = (squeeze(allData(:,i,end)))';%zscore(squeeze(allData(:,i,end)))';
        XXtemp = X(trainInds,:)'*X(trainInds,:);
        Xytemp = X(trainInds,:)'*y(trainInds)';
        kernB = Xytemp'/(XXtemp + lambda*diag(diag(XXtemp)));
        yEst = squeeze(kernB)*X(testInds,:)';
        %y = zscore(y,[],2);%y(scrambles);%
        %yEst = zscore(yEst,[],2);
        yOrig(:,i-timePast) = y(testInds);
        yFits(:,i-timePast) = yEst;
        kern(i-timePast,:) = kernB;
        %mseB = mean((y-yEst).^2,2);
        %ccB = mean(y.*yEst,2);%zscore(y(scrambles)).*zscore(yEst),2);
        %stats.kern(i-timePast,:) = kernB;
        %stats.cc(i-timePast) = ccB;
        %stats.mse(i-timePast) = mseB;
    end
    yFits1 = [yFits1; yFits];yOrig1 = [yOrig1; yOrig];
    subplot(2,3,1);plot(yOrig','b');hold on;plot(mean(yFits),'r','linewidth',2);hold off;
    subplot(2,3,4);plot(yFits');
    subplot(2,3,2);plot(squeeze(allData(testInds,:,1))');
    subplot(2,3,5);plot(squeeze(allData(testInds,:,2))');
    subplot(2,3,3);imagesc(kern',[-1 1]*.2);drawnow;
    %input('');
end