function visScripts(fn,tr)

file = dir('*.mat');
load(file(fn).name,'data');

inds = [6 10:21];
temp = [data{2:end,2}];
f = find(temp == tr);

clear temp;

for j = 1:numel(f)
temp(j,1:numel(data{f(j)+1,inds(1)}),:) = reshape([data{f(j)+1,inds}],numel(data{f(j)+1,inds(1)}),[]);
end
temp(temp == 0) = nan;

pupSize = max(1,(squeeze(temp(:,:,5))-40)*2);
stims = temp(:,:,1);
inds = [2 3 7 8];
temp = temp(:,:,inds);
for i = 1:4
    t = squeeze(temp(:,:,i));
    temp(:,:,i) = temp(:,:,i)-nanmean(t(:));
    temp(:,:,i) = temp(:,:,i)/nanstd(t(:));
end

lims = [-3 3];
step = 1;
figure;
for i = 1:size(temp,1)
    cla%hold off %subplot(4,1,1:3);
    f = find(~isnan(stims(i,:)));
    for j = min(f)-10:step:max(f)+10
        %subplot(4,1,1:3);
        if j < 160
            col = 'b';
        elseif isnan(stims(i,j))
            col = 'g';
        else
            col = 'r';
        end
        biggest = min(size(temp,2),j+step-1);
        scatter(squeeze(temp(i,j:biggest,1)),((j:biggest)-min(f))/20+lims(1),pupSize(i,j:biggest),col,'filled');%squeeze(temp(i,j:biggest,2))
        scatter(squeeze(temp(i,j:biggest,3)),squeeze(temp(i,j:biggest,4)),col,'s');
        plot([stims(i,j) stims(i,j+1)]/5,([j j+1]-min(f))/20+lims(1),'k','linewidth',1);
        set(gca,'xlim',lims,'ylim',lims);
        hold on;
        %subplot(4,1,4);
        %plot([stims(i,1:biggest)/3; temp(i,1:biggest,1)]');
        %set(gca,'xlim',[min(f)-10 max(f)+10],'ylim',lims);
        drawnow;pause(.05);
    end
end