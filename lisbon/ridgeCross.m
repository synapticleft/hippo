function [cc, mse ,kern] = ridgeCross(y,x,numCross,ridges,bootStrap) % stats
%% cross-validated ridge regression. If you bootstrap, no cross-validation performed, use numCross = 1

samples = size(x,1);
scramble = randperm(size(x,1));
x = x(scramble,:);y = y(scramble,:);
cc = zeros(numCross,numel(ridges),size(y,2));mse = cc;
kern = zeros(numCross,numel(ridges),size(x,2),size(y,2));
for i = 1:numCross
    testInds = max(1,ceil([(i-1) i]/numCross*samples));
    testInds = testInds(1):testInds(2);
    trainInds = setdiff(1:samples,testInds);
    if numCross == 1 trainInds = testInds; end
    xTest = x(testInds,:);
    xTrain = x(trainInds,:);
    yTest = y(testInds,:);
    yTrain = y(trainInds,:);
    xCov = xTrain'*xTrain;
    for j = 1:numel(ridges)
        kern(i,j,:,:) = conj(xCov + ridges(j)*diag(diag(xCov)))\(yTrain'*xTrain).';
        %if exist('useKern','var')
        %   yEst = xTest*useKern';
        if size(y,2) == 1
           yEst = xTest*squeeze(kern(i,j,:,:))';
        else
            yEst = xTest*conj(squeeze(kern(i,j,:,:)));
        end
        cc(i,j,:) = diag(corr(yTest,yEst));
        mse(i,j,:) = mean(abs(yTest-yEst).^2)/max(.001,var(yTest(:)));
    end
end
kern = squeeze(mean(kern,1));
cc = squeeze(mean(cc,1));
mse = squeeze(mean(mse,1));
% if exist('bootStrap','var')
%     xCov = x'*x;
%     x = (xCov + ridge*diag(diag(xCov)))\x';
%     %vy = max(.001,var(y));
%     ys = zeros(bootStrap,numel(y));
%     for i = 1:bootStrap
%         scramble = randperm(numel(y));
%         ys(i,:) = y(scramble);
%     end
%     kernB = ys*x';
%     yEst = kernB*x;
%     ys = zscore(ys);
%     ccB = mean(ys.*zscore(yEst),2);
%     mseB = mean((ys-yEst).^2,2);
%     stats.kernL = squeeze(prctile(kernB,2.5))';
%     stats.kernH = squeeze(prctile(kernB,97.5))';
%     stats.ccL = prctile(ccB,2.5);
%     stats.ccH = prctile(ccB,97.5);
%     stats.mseL = prctile(mseB,2.5);
%     stats.mseH = prctile(mseB,97.5);
% end