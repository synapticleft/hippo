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
    plot(d(i,f,1),d(i,f,2),'k');%set(gca,'xlim',[min(squeeze(d(i,f,1))) max(squeeze(d(i,f,1)))],...
    %    'ylim',[min(squeeze(d(i,f,2))) max(squeeze(d(i,f,2)))]);
    set(gca,'xlim',[min(xs{1}) max(xs{1})], 'ylim',[min(xs{2}) max(xs{2})]);
    hold on;
    dx = squeeze(d(i,f,1));
    dy = squeeze(d(i,f,2));
    dx = filtfilt(gausswin(8),sum(gausswin(8)),dx);
    dy = filtfilt(gausswin(8),sum(gausswin(8)),dy);
    for j = 4:10%numel(f)-1 %ord+1
        xsa = minimumJerk(dx(j),dx(j)-dx(j-1), dx(j)-2*dx(j-1)+dx(j-2),ma(1),0,0,numel(f)-j);
        xsb = minimumJerk(dx(j),dx(j)-dx(j-1), dx(j)-2*dx(j-1)+dx(j-2),mb(1),0,0,numel(f)-j);
        ysa = minimumJerk(dy(j),dy(j)-dy(j-1), dy(j)-2*dy(j-1)+dy(j-2),ma(2),0,0,numel(f)-j);
        ysb = minimumJerk(dy(j),dy(j)-dy(j-1), dy(j)-2*dy(j-1)+dy(j-2),mb(2),0,0,numel(f)-j);
        plot(xsa,ysa,'b');
        plot(xsb,ysb,'r');drawnow;
    end
    hold off;
%    input('');
end