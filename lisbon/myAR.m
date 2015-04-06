function myAR(fn,ILD,ord)

inds = [15 16];

if ~exist('ILD','var')
    ILD = [.5 2 4];
end
%cols = colormap('jet');
[d,~,~,center_ILD] = preProcessRoberto(fn,[inds 6],0,0,[],[-ILD ILD]);
X = [];y = [];
for i = 1:size(d,1)
    f = find(squeeze(d(i,:,end)) ~= 0);
    for j = 1:ord
        temp(numel(inds)*(j-1)+(1:numel(inds)),:) = squeeze(d(i,f-j,1:end-1))';
    end
    X = [X temp];temp = [];
    y = [y squeeze(d(i,f,1:2))'];
end
w = y/X
w1 = pascal(ord+1,1);
w1 = -w1(end,2:end);
xx = [w1; zeros(size(w1))];
yy = xx([2 1],:);
w1 = [xx(:) yy(:)]';
w = w1;
for i = 2:ord
    w((i-1)*numel(inds)+(1:numel(inds)),(i-2)*numel(inds)+(1:numel(inds))) = eye(numel(inds));
end
figure;
for i = 1:size(d,1)
    f = find(squeeze(d(i,:,end)) ~= 0);
    plot(d(i,f,1),d(i,f,2),'k');set(gca,'xlim',[min(squeeze(d(i,f,1))) max(squeeze(d(i,f,1)))],...
        'ylim',[min(squeeze(d(i,f,2))) max(squeeze(d(i,f,2)))]);
    hold on;
    %temp = [];
    %temp(1:ord,:) = squeeze(d(i,f(1:ord),1:2));
    for j = ord+1:numel(f)-1 %ord+1
        temps = [];
        temp = [];
        for k = 1:ord
            temp(numel(inds)*(k-1)+(1:numel(inds))) = squeeze(d(i,f(j)-k,1:end-1))';
        end
        for k = j+(1:10)%numel(f)
            temps(k,:) = temp(1:2);
            temp = temp*w';
        end
        plot(temps(j+1:end,1),temps(j+1:end,2));drawnow;%pause(.5);
    end
    hold off;
    input('');
end