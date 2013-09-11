function plotMultiRun(pos,v,Xf,W,deW)
%% running ica multiple times to validate consistency of peak locations

pos = fixPos(pos);
Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).')];
Xf = Xf(:,~any(isnan(pos')));
pos = pos(~any(isnan(pos')),:);
vel = angVel(pos);
b = nan*ones(size(pos,1),1);
bounds = [.1 .9];
Xf = Xf(:,vel(:,1) >= .05);
pos = pos(vel(:,1) >= .05,1);
pos = pos-min(pos)+eps;
pos = pos/max(pos);
b(pos < bounds(1)) = -1;b(pos > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
w = watershed(b==0);
w = w-1; 
pos(mod(w,2) ==1 ,1) = -pos(mod(w,2) ==1 ,1) + 2*max(pos(:));
pos = ceil(pos*100);
angCol = colormap(hsv(64));
h = zeros(1,200);
if exist('deW','var')
    numLoop = numel(W);
else
    numLoop = size(W,1);
end
figure;for i = numLoop-1:-1:1
    if ~exist('deW','var')
        acts = squeeze(W(i,:,:))*Xf;
    else
        acts = (W{end-i+1}'*deW(end-i:end,:))*Xf;
    end
t1 = zeros(size(acts,1),200);
for j = 1:size(acts,1)
%acts = acts.';acts = acts(:);
%pos = repmat(pos,size(W,2));xs = meshgrid(1:size(W,2),1:size(Xf,2));
t1(j,:) = abs(accumarray(pos,acts(j,:),[],@mean));
end
t1 = t1(max(t1,[],2)>1,:);
t1(t1 < 1) = nan;
[~,peakLoc] = max(t1,[],2);
h(peakLoc) = h(peakLoc) + 1;
c = peakLoc/size(t1,2);%rand(1,numel(peakLoc));%
cc = angCol(max(1,ceil(c*64)),:);
set(gca,'nextPlot','add','ColorOrder',cc,'Color',[0 0 0],'xticklabel',[],'yticklabel',[],'fontsize',16);
hold on;plot(i/2+t1','linewidth',2);drawnow;
end
axis tight;
angC = colormap(hsv(200));
figure;hold on;for i = 1:200
bar(i,h(i),'facecolor',angC(i,:),'edgecolor',angC(i,:));
end
axis tight;set(gca,'color','k');