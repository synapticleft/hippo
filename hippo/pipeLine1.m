function [cc mse kern snr f] = pipeLine1(y,x,numCross,ridge,fs,numX,useKern)
figsOn = 0;
samples = size(x,1);
%x = bsxfun(@minus,x,mean(x));
%x = bsxfun(@rdivide,x,std(x));
%for testing
%tk = sin(linspace(0,pi,size(x,2)));
%y = (x*tk')'; y = y + randn(size(y))*std(y(:))/2;
%y = bsxfun(@minus,y,mean(y,2));
cc = zeros(numCross,size(y,1));
%params.Fs = fs;params.fpass = [0 70];params.tapers = [5 9];
kern = zeros(numCross,size(x,2),size(y,1));
for i = 1:numCross
    testInds = max(1,ceil([(i-1) i]/numCross*samples));
    testInds = testInds(1):testInds(2);
    trainInds = setdiff(1:samples,testInds);
    if numCross == 1 trainInds = testInds; end
%    [min(testInds) max(testInds) samples]
    xTest = x(testInds,:);
    xTrain = x(trainInds,:);
    yTest = y(:,testInds);
    yTrain = y(:,trainInds);
    xCov = xTrain'*xTrain;
    xCov = xCov + ridge*eye(size(xCov,1));
        kern(i,:,:) = conj(xCov)\(conj(yTrain)*xTrain).';
        if exist('useKern','var')
            yEst = xTest*useKern';
        elseif size(y,1) == 1
            yEst = xTest*squeeze(kern(i,:,:))';
        else
            yEst = xTest*conj(squeeze(kern(i,:,:)));
        end
%    [gamma,~,~,~,~,f] = coherencyc(yTest',yEst,params);
%    gammasq = gamma.*conj(gamma);
%    snr(i,:,:) = gammasq./(1-gammasq);
    cc(i,:) = diag(corr(yTest.',yEst));
    mse(i,:) = mean(abs(yTest.'-yEst).^2)/mean(abs(yTest.').^2);
%     try
%     sigSpec(i,1,:) = mtspectrumc(yTest,params);
%     sigSpec(i,2,:) = mtspectrumc(yEst',params);
%     [sigSpec(i,3,:) f] = mtspectrumc(yTest-yEst',params);
%     catch
%         sigSpec(i,1:3,:) = repmat(f,[3 1]);
%     end
    %snr(i,:,:) = mtspectrumc(yTest-yEst',params);
%      subplot(311);plot(squeeze(kern(i,:)));
%      subplot(313);plot(f,squeeze(snr(i,:)));
%      drawnow;pause(.5);
end
if figsOn
figure;
temp = squeeze(mean(kern,1));%temp = temp(:,1);
if exist('useKern','var') temp = useKern; end
temp = reshape(temp,[numel(temp)/numX numX]);
subplot(311);plot(real(temp));axis tight;
%temp1 = complex(yTest,-yEst');
subplot(312);plot(real(yTest)');hold all;plot(real(yEst));%plot(yTest'-yEst);%sPlot(temp1,0);%hold off;plot(yTest');hold all;plot(yEst);axis tight;
axis tight;snr = squeeze(mean(snr,1));
subplot(313);%plot(f,snr);
%[S f] = %x(:,1)
%%hold all;plot(f,S);
%[S f] = mtspectrumc(yTest-yEst',params);S = S;S = S/max(S)*max(snr(:));
sigSpec = squeeze(mean(sigSpec,1));
%params.trialave = 1;
plot(f,sqrt(sigSpec));hold all;plot(f,sqrt(snr/20));
[S f] = mtspectrumc(xTest(:,1)/std(xTest(:))*std(yTest(:)),params);plot(f,sqrt(S));axis tight;
% subplot(414);%params.tapers = [3 5];
% [S f] = mtspectrumc(temp,params);
% %imagesc(f,1:size(S,2),log(S)');
% plot(mean(S,2));
%figure;scatter(yTest,yEst');
end
cc = squeeze(mean(cc,1));
mse = squeeze(mean(mse,1));