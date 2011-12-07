function phaseDisplay(data,dims)

sub = 2;trace = 1000;
if nargin < 2
    dims = [8 8];
end

[u s v] = svds(data,sub);
circ_mean(circ_dist(angle(v(1:(end-1),:)),angle(v(2:end,:))))
theta = angle(data);
amp = abs(data);
data1 = u*s*v';
theta1 = angle(data1-data);
amp1 = abs(data1-data);
scale = 3*std(amp(:));
for i = 1:(size(data,2)-trace)
    subplot(231);imagesc(reshape(theta(:,i),dims),[-pi pi]), axis off, colormap hsv
    alpha(min(reshape(amp(:,i)/scale,dims),1));
    subplot(232);imagesc(reshape(theta1(:,i),dims),[-pi pi]), axis off, colormap hsv
    alpha(min(reshape(amp1(:,i)/scale,dims),1));
    subplot(233);
    plot(v(i+(1:trace/100)-1,:));hold all;
    set(gca,'xlim',[-.1 .1]/4,'ylim',[-.1 .1]/4);
    plot(v(i+(1:trace/100)-1,1).*conj(v(i+(1:trace/100)-1,2))*100);
    %plot(v(i+(1:trace/100)-1,2)./exp(1i*angle(v(i+(1:trace/100)-1,1))));%hold on;
    %plot(abs(v(i+(1:trace)-1,1))+.000001*1i,'r');
    hold off;
    subplot(212);
    hold off;
    sPlot([abs(v(i+(1:trace)-1,:)) circ_dist(angle(v(i+(1:trace)-1,1)),angle(v(i+(1:trace)-1,2))) ...
        v(i+(1:trace)-1,1).*conj(v(i+(1:trace)-1,2)) circ_std(theta(:,i+(1:trace)-1))']',[],0,1);
    drawnow;
    
end

%freezeColors