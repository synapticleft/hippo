function handEyeVis(fn,ILD)
figure;
if ~exist('ILD','var')
    ILD = [.5 2 4];
end
cols = colormap('jet');
[d,~,~,center_ILD] = preProcessRoberto(fn,[10 12 15 16 6],0,0,[],[-ILD ILD]);
inds = groupILDs(d,center_ILD,1);
for i = 11:max(inds)/2
fmax = 0;
%subplot(132);
subplot(1,6,[4 5]);
jj = find(inds == 2*i -1);
for j = 1:numel(jj)
f = find(squeeze(d(jj(j),:,end)) ~= 0);
fmax = max(numel(f),fmax);
scatter(d(jj(j),f,3),1:numel(f),[],cols(round(max(1,min(64,d(jj(j),f,1)*20+32))),:),'filled');%max(5,100*d(jj(j),f,2))
hold on;axis tight;
end
hold off;
%subplot(131);
subplot(1,6,[2 3]);
jj = find(inds == 2*i);
for j = 1:numel(jj)
f = find(squeeze(d(jj(j),:,end)) ~= 0);
fmax = max(numel(f),fmax);
scatter(d(jj(j),f,3),1:numel(f),[],cols(round(max(1,min(64,d(jj(j),f,1)*20+32))),:),'filled');%max(5,100*d(jj(j),f,2))
hold on;
end
hold off;axis tight;
%subplot(131);set(gca,'ylim',[0 fmax]);
%subplot(132);set(gca,'ylim',[0 fmax]);
subplot(1,6,[2 3]);set(gca,'ylim',[0 fmax]);
subplot(1,6,[4 5]);set(gca,'ylim',[0 fmax]);
%subplot(133);
for j = 1:numel(jj)    
f = find(squeeze(d(jj(j),:,end)) ~= 0);
subplot(1,6,1);
plot(squeeze(d(jj(j),f,end)),1:numel(f),'b','linewidth',2);hold all;
axis tight;set(gca,'ylim',[0 fmax]);
subplot(1,6,6);
hold all;
plot(-squeeze(d(jj(j),f,end)),1:numel(f),'r','linewidth',2);%hold off;
axis tight;set(gca,'ylim',[0 fmax]);
end
subplot(1,6,1);hold on;plot([0 0],[0 fmax],'k--');set(gca,'xtick',[],'ytick',[0 60],'yticklabel',{'0','.75'});ylabel('Time (s)');
subplot(1,6,[2 3]);hold on;plot([0 0],[0 fmax],'k--');set(gca,'xtick',[],'ytick',[]);
subplot(1,6,[4 5]);hold on;plot([0 0],[0 fmax],'k--');set(gca,'xtick',[],'ytick',[]);
subplot(1,6,6);hold on;plot([0 0],[0 fmax],'k--');set(gca,'xtick',[],'ytick',[]);

%title(i);
input('');
end