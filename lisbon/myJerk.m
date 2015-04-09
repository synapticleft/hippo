function myJerk(fn,ILD)

inds = [15 16];

if ~exist('ILD','var')
    ILD = [.5 2 4];
end
%cols = colormap('jet');
[d,~,~,center_ILD] = preProcessRoberto(fn,[inds 6],0,0,[],[-ILD ILD]);

xs{2} = -1:.05:2;xs{1} = -3:.05:3;

h = hist3([makeFlat(d(:,:,1)) makeFlat(d(:,:,2))],xs);
h = imfilter(h,fspecial('gaussian',10,3));
ma = max(max(h(xs{1} < 0,xs{2} > 0)));
mb = max(max(h(xs{1} > 0,xs{2} > 0)));
[xa,ya] = find(h == ma);
[xb,yb] = find(h == mb);
ma = [xs{1}(xa) xs{2}(ya)];
mb = [xs{1}(xb) xs{2}(yb)];

figure;
for i = 1:size(d,1)
    f = find(squeeze(d(i,:,end)) ~= 0);
    %plot(d(i,f,1),d(i,f,2),'k');%set(gca,'xlim',[min(squeeze(d(i,f,1))) max(squeeze(d(i,f,1)))],...
    %    'ylim',[min(squeeze(d(i,f,2))) max(squeeze(d(i,f,2)))]);
    dx = squeeze(d(i,f,1));
    dy = squeeze(d(i,f,2));
    dx = filtfilt(gausswin(8),sum(gausswin(8)),dx);
    dy = filtfilt(gausswin(8),sum(gausswin(8)),dy);
    if 0
    mse = zeros(numel(dx),numel(dx),4);cost = mse;costHat = cost;
    for j = 1:numel(dx)
        for k = j+3:numel(dx)
            [msex,costx] = myMin(dx(j:k));%minimumJerk(dx(j),dx(j)-dx(j-1), dx(j)-2*dx(j-1)+dx(j-2),dx(k),dx(k)-dx(k-1), dx(k)-2*dx(k-1)+dx(k-2),k-j);
            [msey,costy] = myMin(dy(j:k));
            mse(j,k,:) = msex+msey;
            cost(j,k,:) = costx + costy;
            %costHat(j,k,:) = costxHat + costyHat;
            %[yc,yj] = myMin(dy(j:k));%minimumJerk(dy(j),dy(j)-dy(j-1), dy(j)-2*dy(j-1)+dy(j-2),dy(k),dy(k)-dy(k-1), dy(k)-2*dy(k-1)+dy(k-2),k-j);
        end
    end
    for j =1:4
        subplot(2,5,j);imagesc(-log(mse(:,:,j)./cost(:,:,j)));
        subplot(2,5,j+5);imagesc(-log(cost(:,:,j)));
%            subplot(3,4,j+8);imagesc(-log(costHat(:,:,j)));
    end
    subplot(2,5,10);plot(dx,dy);
    else
        plot(dx,dy,'k','linewidth',2);
       set(gca,'xlim',[min(xs{1}) max(xs{1})], 'ylim',[min(xs{2}) max(xs{2})]);
    hold on;
    for j = 4:numel(f)-1 %ord+1
        xsa = minimumJerk(dx(j),dx(j)-dx(j-1), dx(j)-2*dx(j-1)+dx(j-2),ma(1),0,0,numel(f)-j);
        xsb = minimumJerk(dx(j),dx(j)-dx(j-1), dx(j)-2*dx(j-1)+dx(j-2),mb(1),0,0,numel(f)-j);
        ysa = minimumJerk(dy(j),dy(j)-dy(j-1), dy(j)-2*dy(j-1)+dy(j-2),ma(2),0,0,numel(f)-j);
        ysb = minimumJerk(dy(j),dy(j)-dy(j-1), dy(j)-2*dy(j-1)+dy(j-2),mb(2),0,0,numel(f)-j);
        plot(xsa,ysa,'b');
        plot(xsb,ysb,'r');drawnow;
    end
    end
    hold off;
    input('');
end

function [mse,cost] = myMin(trajx) %,costHat
    trajxhat = minimumJerk(trajx(3),trajx(3)-trajx(2),trajx(3)-2*trajx(2)+trajx(1),trajx(end),trajx(end)-trajx(end-1),trajx(end)-2*trajx(end-1)+trajx(end-2),numel(trajx)-2);
    mse(1) = sum((trajx(3:end) - trajxhat).^2);
    cost(1) = numel(trajx);%sum(trajx(3:end).^2);
    for i = 1:3
        mse(i+1) = sum((diff(trajx(3:end),i) - diff(trajxhat,i)).^2);
        cost(i+1) = sum(diff(trajx(3:end),i).^2);
        %costHat(i) = sum(diff(trajxhat,i).^2);
    end