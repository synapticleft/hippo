function plotMultiRun1(acts,runs)
thresh = 1;
angCol = colormap(hsv(size(acts,ndims(acts))));
h = zeros(1,size(acts,ndims(acts)));
figure;subtightplot(5,1,1:4);hold on;
if exist('runs','var')
    numRuns = max(runs);
    acts(:,[1 end]) = 0;
else
    numRuns = size(acts,1);
    acts(:,:,[1 end]) =0;
end
for i = numRuns:-1:1
    if exist('runs','var')
        t1 = acts(runs == i,:);
    else
        t1 = squeeze(acts(i,:,:));
    end
    t1 = abs(filtfilt(gausswin(5),sum(gausswin(5)),t1')');
    t1 = t1(max(t1,[],2)>thresh,:);
    t1(t1 < thresh) = nan;
    [~,peakLoc] = max(t1,[],2);
    h(peakLoc) = h(peakLoc) + 1;
    c = peakLoc;%rand(1,numel(peakLoc));%
    cc = angCol(c,:);%max(1,ceil(c*64)),:);
    scatter(peakLoc,i*ones(1,numel(peakLoc)),100,cc,'filled');drawnow;
    %set(gca,'nextPlot','add','ColorOrder',cc);%,'Color',[0 0 0],'xticklabel',[],'fontsize',16);
    %hold on;plot(i+2*t1','linewidth',2);drawnow;
end
set(gca,'color',[0 0 0],'xlim',[0 size(acts,ndims(acts))],'ylim',[0 numRuns],'xticklabel',[],'fontsize',16);
plot([1 1]*size(acts,ndims(acts))/2,[0 numRuns+2],'w--','linewidth',2);
ylabel('runs');
angC = colormap(hsv(200));
subtightplot(5,1,5);hold all;
for i = 1:size(acts,ndims(acts))
    bar(i,h(i),'facecolor',angC(i,:),'edgecolor',angC(i,:));
end
plot([1 1]*size(acts,ndims(acts))/2,[0 max(h)],'w--','linewidth',2);
set(gca,'fontsize',16);
ylabel('counts');
xlabel('position');
axis tight;set(gca,'color','k');