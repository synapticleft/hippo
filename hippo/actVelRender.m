function actVelRender(acts,inds,pos)
%renders activations with arrows pointing in direction of movement
global ColorOrder
vel = diff(pos(:,1:2));
vel = filtLow([0 0;vel]',1250/32,2)';

f1 = figure;
f2 = figure;
dimx = ceil(sqrt(numel(inds)));
dimy = ceil(numel(inds)/dimx);
thresh = 1;
for i= 1:inds
    fi = abs(acts(inds(i),:)) > thresh;
    temp = complex(vel(fi,1),vel(fi,2));
    ColorOrder = hsv2rgb([(angle(temp)/(2*pi) +.5) ones(sum(fi),1) min(1,power(abs(temp),1))]);
    figure(f1);subplot(dimy,dimx,i);arrow3(pos(fi,1:2),pos(min(find(fi)+5,size(pos,1)),1:2)+eps,...
        'o',(abs(acts(inds(i),fi))-thresh/2)/5,(abs(acts(inds(i),fi))-thresh/2)/5);%drawnow
    %figure(f2);subplot(dimy,dimx,i);arrow3(zeros(sum(fi),2),pos(min(find(fi)+5,size(pos,1)),1:2)+eps-pos(fi,1:2),...
        %'o',(abs(acts(inds(i),fi))-thresh/2)/5,(abs(acts(inds(i),fi))-thresh/2)/5);drawnow
        figure(f2);subplot(dimy,dimx,i);hist(angle(temp),100);drawnow;
end