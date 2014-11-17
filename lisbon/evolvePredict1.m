function cc = evolvePredict1(fn,lambda,inds,timePast,trialHist,figsOn) %sessNorm,
% fit LDA and regularized linear regression of start- and end- aligned
% trials, look at performance, and weights.

% choices: 1) re-normalize measures in each session to correct for drifts
% (sessNorm) NOW I DO THIS FOR EYE DATA BUT NOT OTHERS
% 2) regularization of fit (lambda)
% 3) fit choice, previous choice, or right answer
% 4) number of time steps in the past
% 5) which measurements to include
% 6) separate each difficulty level or combine them all

[allData allOut allOutShift] = preProcessRoberto(fn,inds,timePast,0,[]);%,[-.5 .5]);

y = allOut(:,1)';

for i = timePast+1:size(allData,2)
    X = allData(:,i-timePast:i,:);
    X = X(:,:);
    if trialHist
        X = [X allOutShift];%circshift([answer; binAnswer; choice]',[1 0]) circshift([answer; binAnswer; choice]',[-1 0])];
    end
    X = zscore(X);
    coeff(i-timePast,:) = (X'*y')'/(X'*X + lambda*size(X,1)*eye(size(X,2)));%y/X';
    yHat = coeff(i-timePast,:)*X';
    mse(i-timePast) = mean(sqrt((y-yHat).^2));%sum(yHat'+1 ~=info(inds,class))/sum(inds);
    temp = corrcoef(y,yHat);
    cc(i-timePast) = temp(1,2);
    if 0 % figsOn
    if (i-timePast == 70) || (i-timePast == 120)
        figure;scatter(allOut(:,1),yHat,'filled');hold all;
        %title(i);
    end
    ylim([-3 3]);drawnow;
    %m(i-timePast) = getframe(gcf);
    end
end
%movie2avi(m,'timeFit.avi');
if figsOn
%coeff = reshape(coeff,size(coeff,1),timePast,[]);
figure;imagesc(coeff(:,:)');
figure;plot(max(mse)-mse);hold all;plot(cc);ylim([0 1]);
[~,m] = min(mse);
%figure;plot(squeeze(coeff(90,:,:)));
end