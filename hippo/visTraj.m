function visTraj(pos,sp,v)

dSkip = 1;
steps = 5;
decay = .99;
Fs = 1250/32;
skip = steps;
if size(v,2) < size(v,1)
    v = v.';
end

col = rand(size(sp,1),3);
angCol = colormap('hsv');
absCol = colormap('jet');
pos(pos == -1) = nan;
nanInds = any(isnan(pos'));
pos = pos(~nanInds,:);
pos = filtLow(pos',Fs,4,8)';
pos = bsxfun(@minus,pos,min(pos));
pos = bsxfun(@rdivide,pos,max(pos));
v = bsxfun(@rdivide,v,std(v,0,2));
[~,fsort] = sort(sum(sp,2),'descend');
%% place fields
% figure;

% pos((end+1):size(sp,2),:) = 0;
% bins{1} = linspace(0,1,20);bins{2} = bins{1};
% h = hist3(pos(:,1:2),bins);
% h1 = hist3(pos(:,3:4),bins);
% for i = 1:size(sp,1)
%     hS(:,:,1) = hist3(pos(sp(inds(i),:) > 0,1:2),bins)./h;
%     hS(:,:,2) = hist3(pos(sp(inds(i),:) > 0,3:4),bins)./h1;
%     hS(:,:,3) = 0;
%     hS = hS/max(hS(:));
%     imagesc(hS);drawnow;pause(1);
%     %hold all;
%     %scatter(pos(sp(inds(i),:)>0,1),pos(sp(inds(i),:)>0,2),col(i),'filled');drawnow;
% end
% return
%% trajectory
figure;  h = axes('Color','k');
im = getframe(h);%  hold all;
for i = 1:skip:size(pos,1)
    image([0 1],[0 1],im.cdata*decay);hold all;%
    inds = i:(i+steps);
dist = 100*sqrt(sum((pos(i,1) - pos(i+dSkip,1)).^2));
dist2 = 100*sqrt(sum((pos(i,3) - pos(i+dSkip,3)).^2));
    scatter(pos(inds,1),pos(inds,2),50*abs(v(2,inds)).^2,absCol(abs2Col(dist),:),'filled');%
    scatter(pos(inds,3),pos(inds,4),50*abs(v(2,inds)).^2,absCol(abs2Col(dist2),:),'filled');%
%    scatter(pos(inds,1),pos(inds,2),50*abs(v(1,inds)).^2,absCol(abs2Col(sum(sp(fsort(3:end),inds)/10)),:),'filled');%
%    scatter(pos(inds,3),pos(inds,4),50*abs(v(2,inds)).^2,angCol(phase2Col(circ_dist(v(1,inds),v(2,inds))),:),'filled');%
%    hold off;
    set(gca,'xlim',[0 1],'ylim',[0 1]);%    title(round(num2str(i/Fs)));
    im = getframe(h);
    drawnow;
hold off;
end

function c = phase2Col(ang)
c = ceil((ang+pi)/(2*pi)*64);

function c = abs2Col(a)
if isnan(a)
    a = 0;
end
c = 1 + floor(64*min(.99,a));