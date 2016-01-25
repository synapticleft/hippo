function [t1 t2 vels vels1 poses] = runTriggerView1d(pos,v,Xf,accumbins,thresh,r)
%% 1d binning of ICA activations on linear track
Xf = bsxfun(@times,Xf,exp(1i*angle(v(:,1))).');
if exist('r','var')
    r1 = pinv(r);
    t = r*Xf;
else
    t = Xf;
    r1 = ones(size(Xf,1));
end
clear Xf;
[~,pos,inds,runs,vel,b] = fixPos(pos,thresh);
inds(size(t,2)+1:end) = 0;
offset = [-50 299];
t2 = zeros(size(t,1),max(runs),(range(offset)+1)*2);
vels1 = zeros(max(runs),(range(offset)+1)*2);
poses = zeros(size(vels1));
for i = 1:max(runs)
    ind = find(runs == i & b > 0,1);
    ind1 = find(runs == i & b < 0,1);
    t2(:,i,:) = [t(:,ind+(offset(1):offset(2))) t(:,ind1+(offset(1):offset(2)))];
    vels1(i,:) = [vel(ind+(offset(1):offset(2)))' vel(ind1+(offset(1):offset(2)))'];
    poses(i,:) = [pos(ind+(offset(1):offset(2)))' pos(ind1+(offset(1):offset(2)))'];
    ind1 = find(runs == i & b < 0,1,'last')-ind1-offset(1);
    ind = find(runs == i & b > 0,1,'last')-ind-offset(1);
    vels1(i,1+(ind+50:range(offset))) = 0;
    poses(i,1+(ind+50:range(offset))) = nan;
    vels1(i,range(offset)+1+ind1+50:end) = 0;
    poses(i,range(offset)+1+ind1+50:end) = nan;
end
figure;%subplot(212);
%set(gca,'color','k','nextPlot','add','ColorOrder',jet(max(runs)),'ylim',[0 6],'fontsize',16);
%plot(vels1','linewidth',2);hold all;plot([1 1]*range(offset),[-1 6],'w','linewidth',4);
subplot(221);
%subplot(211);
set(gca,'color','k','nextPlot','add','ColorOrder',jet(double(max(runs))),'ylim',[0 1],...
    'fontsize',16,'xtick',[0 350],'xticklabel',{'0',num2str(round(350*32/1250*10)/10)},...
    'ytick',[0 1],'yticklabel',[0 250]);
poses(:,range(offset)+1:end) = poses(:,range(offset)+1:end)+1;
poses = poses - 1;
plot(poses','linewidth',2);hold all;plot([1 1]*range(offset),[0 1],'w','linewidth',4);
ylabel('Position (cm)');
xlabel('Time (s)');%title('Velocity');
subplot(222);
imagesc(vels1*1250/32/2.7,[-10 90]);hold all;plot([1 1]*range(offset),[0 1+max(runs)],'w','linewidth',4);
set(gca,'fontsize',16,'xtick',[]);
ylabel('Trial #');title('Velocity');
%%%%%%%%%
posd = max(1,ceil(pos*accumbins));
t1 = zeros(size(t,1),max(runs),accumbins(1)*2);
vels = accumarray([runs(inds); posd(inds)']',vel(inds),[max(runs) 2*accumbins(1)] ,@mean);
for j = 1:size(t,1)
         t1(j,:,:) = accumarray([runs(inds); posd(inds)']',t(j,inds),[max(runs) 2*accumbins(1)] ,@mean);
         t1(j,:,:) = t1(j,:,:)*exp(1i*angle(mean(r1(:,j))));
end
im = superImp(t1,[],0,5);
im = [im(:,accumbins+1:end,:) im(:,accumbins:-1:1,:)];
subplot(223);image(im);set(gca,'fontsize',16,'xtick',[1 100],'xticklabel',{'0','250'});
ylabel('Trial #');xlabel('Position (cm)');title('Position-locked FFPs');
hold all;plot([1 1]*accumbins,[0 1+max(runs)],'w','linewidth',4);
im = superImp(t2,[],0,5);
subplot(224);image(im);set(gca,'fontsize',16,'ytick',[],'xtick',[1 350],...
    'xticklabel',{'0',num2str(round(350*32/1250*10)/10)});
title('Time-locked FFPs');xlabel('Time (s)');
hold all;plot([1 1]*range(offset),[0 1+max(runs)],'w','linewidth',4);