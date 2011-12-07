function [a b inds] = crcns(trace,f1,f2,dims,ratio)
%this program visualizes spatial LFP in two frequency bands and 
%shows the direction in which they move.

global step
global h

Fs = 1250/5;%/32;
figure('Name','Traveling LFP waves in rat hippocampus');
traceT = morFilter(trace,f1,Fs);%hipFilter(trace,f1(1),f1(2),Fs);%
HT = angle(traceT);%hilbert(traceT'))';
traceG = morFilter(trace,f2,Fs);%hipFilter(trace,f2(1),f2(2),Fs);%
HG = angle(traceG);%hilbert(traceG'))';
traceT = real(traceT);traceG = real(traceG);
center = dims([2 1])/2;
scale = 20;
numPast = 20;
scaleTG = 5;
histT = zeros(numPast,2);
histG = zeros(numPast,2);
step = 1;
h = uicontrol('Style','slider','Position', [20 20 100 20],'Value',2,'Min',1,'Max',20, 'Callback',@fixStep);
%tRange = [min(traceT(:)) max(traceT(:))]; gRange = [min(traceG(:)) max(traceG(:))];
tRange = [-1 1]*sqrt(mean(var(traceT(:))))*3; gRange = [-1 1]*sqrt(mean(var(traceG(:))))*3;
colormap gray;
for j = 1:inf
for i = 1:size(trace,2)
    [xt yt] = myGradient(reshape(HT(:,i),dims));
    xt = -xt;yt = -yt;tm = [mean(xt(:)) mean(yt(:))];
    [xg yg] = myGradient(reshape(HG(:,i),dims));
    xg = -xg; yg = -yg; gm = [mean(xg(:)) mean(yg(:))];
    if mod(i,step) == 0
        histG = circshift(histG,[-1 0]);
        if mod(i,step*scaleTG) == 0 
            histT = circshift(histT,[-1 0]);
            histT(numPast,:) = center([2 1]).*ratio+tm*scale;
        end
        histG(numPast,:) = center([2 1]).*ratio+gm*scale;
        %subplot('Position',[0 .8 1 .2]);hold on;
        
        %plot(mod([i-1+tml(1),i+tm(1)],500),[tml(2) tm(2)],'k','LineWidth',2);
        subplot('Position',[0 0 .48 1]);%1,2,1);
        imagesc((1:dims(1))*ratio(1),(1:dims(2))*ratio(2),reshape(traceT(:,i),dims),tRange*.8);hold on;axis image off;
        plot(histT(:,1),histT(:,2),'w','LineWidth',2);
        quiver(xt,yt,'w');
        quiver(center(2)*ratio(1),center(1)*ratio(2),tm(1)*scale,tm(2)*scale,'w','LineWidth',5);
        hold off;
        subplot('Position',[.52 0 .48 1]);%s1);%1,2,2);
        imagesc((1:dims(1))*ratio(1),(1:dims(2))*ratio(2),reshape(traceG(:,i),dims),gRange*.8);hold on;axis image off;
        plot(histG(:,1),histG(:,2),'w','LineWidth',2);
        quiver(xg,yg,'w');
        quiver(center(2)*ratio(1),center(1)*ratio(2),gm(1)*scale,gm(2)*scale,'w','LineWidth',5);
        hold off;
        drawnow;
    end
end
end

function fixStep(a,b)
global step
global h
step = round(get(h,'Value'));

