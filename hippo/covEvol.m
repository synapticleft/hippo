function covEvol(data,u,step,window,dataF)

decay = .99;
[~,inds] = sort(angle(u),'descend');
data = data(inds,:);
tcov = zeros(size(data,1),2*window-1);
if exist('dataF','var')
    dataF = real(dataF(inds,:));
    tcov1 = tcov;
    lims1 = [0 0];
end
c = getCol(size(data,1));
h = figure;plot([3 3]);
lims = [0 0];
for i = 1:step:inf
    figure(h);
    tm = mean(data(:,i+(1:window)));
    for j = 1:size(data,1)
        tcov(j,:) = xcov(data(j,i+(1:window)),tm);
    end
    lims = [min(lims(1)*decay,min(tcov(:))) max(lims(2)*decay,max(tcov(:)))];
    if exist('dataF','var')
        tm = mean(dataF(:,i+(1:window)));
        for j = 1:size(data,1)
            tcov1(j,:) = xcov(dataF(j,i+(1:window)),tm);
        end
        lims1 = [min(lims1(1)*decay,min(tcov1(:))) max(lims1(2)*decay,max(tcov1(:)))];
        subplot(212);
        hold off;plot(window,0);
        %set(gca,'ColorOrder',c,'ylim',lims1);hold on;
        %plot(tcov1');
        imagesc(tcov1);
        subplot(211);
    end
        hold off;plot(window,0);
    %set(gca,'ColorOrder',c,'ylim',lims);hold on;
    %plot(tcov');
    %plot(mean(tcov + fliplr(tcov))/2,'k','linewidth',3);
    imagesc(tcov);
    drawnow;
end

function c = getCol(n)
c = repmat(linspace(0,1,n)',[1 3]);
for i = 1:3
c(:,i) = max(0,min((1-abs(c(:,i)-i/4)*2.5)*1.5,1));
end
c = fliplr(c);