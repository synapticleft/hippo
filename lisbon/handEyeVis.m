function handEyeVis(fn,ILD)
showEye = 1;
figure;
if ~exist('ILD','var')
    ILD = [.5 2 4];
end
cols = colormap('jet');
[d,center_ILD] = preProcessRoberto1(fn,[10 11 15 16 6],0,0,[],[-ILD ILD]);
inds = groupILDs(d,center_ILD,1);
for i = 1:max(inds)/2 %11
fmax = 0;
for k = 0:1
subplot(1+showEye,6,2*k+[2 3]);
jj = find(inds == 2*i - k);
for j = 1:numel(jj)
f = find(squeeze(d(jj(j),:,end)) ~= 0);
fmax = max(numel(f),fmax);
scatter(d(jj(j),f,3),d(jj(j),f,4),[],cols(round(max(1,min(64,d(jj(j),f,1)*20+32))),:),'filled');%max(5,100*d(jj(j),f,2))
hold on;axis tight;
end
hold off;
%subplot(1+showEye,6,2*k+[2 3]);set(gca,'ylim',[0 fmax]);
end
for j = 1:numel(jj)    
f = find(squeeze(d(jj(j),:,end)) ~= 0);
subplot(1+showEye,6,1);
plot(squeeze(d(jj(j),f,end)),1:numel(f),'b','linewidth',2);hold all;
axis tight;set(gca,'ylim',[0 fmax]);
subplot(1+showEye,6,6);
plot(-squeeze(d(jj(j),f,end)),1:numel(f),'r','linewidth',2);hold all;%hold off;
axis tight;set(gca,'ylim',[0 fmax]);
end
hold off;
subplot(1+showEye,6,1);hold on;plot([0 0],[0 fmax],'k--');set(gca,'xtick',[],'ytick',[0 60],'yticklabel',{'0','.75'});ylabel('Time (s)');hold off;
%subplot(1+showEye,6,[2 3]);hold on;plot([0 0],[0 fmax],'k--');set(gca,'xtick',[],'ytick',[]);hold off;
%subplot(1+showEye,6,[4 5]);hold on;plot([0 0],[0 fmax],'k--');set(gca,'xtick',[],'ytick',[]);hold off;
subplot(1+showEye,6,6);hold on;plot([0 0],[0 fmax],'k--');set(gca,'xtick',[],'ytick',[]);hold off;
if showEye
    for k = 0:1
        for l = 1:2
    jj = find(inds == 2*i -k);
    subplot(2,4,4+k*2+l);
    for j = 1:numel(jj)
        f = find(squeeze(d(jj(j),:,end)) ~= 0);
        fmax = max(numel(f),fmax);
        if l == 1
            scatter(d(jj(j),f,l),1:numel(f),[],cols(round(max(1,min(64,d(jj(j),f,2)*20+32))),:),'filled');%max(5,100*d(jj(j),f,2))
        else
            %scatter(d(jj(j),f,l-1),d(jj(j),f,l),[],cols(min(64,1:numel(f)),:),'filled');hold on;
            %plot(d(jj(j),f,l-1),d(jj(j),f,l),'k');
            scatter(d(jj(j),f,3),1:numel(f),[],cols(round(max(1,min(64,d(jj(j),f,1)*20+32))),:),'filled');%
        end
        hold on;axis tight;
    end
    hold off;
        end
    end
end

%title(i);
input('');
end