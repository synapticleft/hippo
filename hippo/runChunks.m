function [allX,W,t,samplePos] = runChunks(X,v,pos)

%%Processing of position information
thresh = .05;bounds = [.1 .9];win = [-2 ceil(1250/8/8)];
pos(pos == -1) = nan;
reject = 0;
for i = 1:4
    reject = reject | min([0; diff(pos(:,i))],flipud([0; diff(flipud(pos(:,i)))])) < -20;
end
pos(reject,:) = nan;
% if size(X,2) < size(pos,1)
%     pos = pos(1:size(X,2),:);
% end
for i = 1:4
    nanInds = find(~isnan(pos(:,i)));
    pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
end
nanInds = isnan(pos(:,1)) | isnan(pos(:,3));
%pos = pos(~nanInds,:);Xf = Xf(:,~nanInds); v= v(~nanInds,:);
vel = angVel(pos);vel = vel(:,1);
vel = [0; vel(:,1)];
pos = bsxfun(@minus,pos,nanmean(pos));
[a,~,c] = svd(pos(~nanInds,1:2),'econ');pos = (c\pos(:,1:2)')';%pos = a;pos(nanInds) = nan;
pos = pos(:,1);
for i = 1:size(pos,2)   
    pos(:,i) = pos(:,i) - min(pos(:,i));
    pos(:,i) = pos(:,i)/(max(pos(:,i)));
    pos(:,i) = min(pos(:,i),.9999);
end
pos(nanInds) = 0;
%%THE filtLow function requires some functions from the signal processing
%%toolbox, but is not particularly necessary.
vel = filtLow(vel,1250/32,1);
vel = vel/max(vel);
vel = resample(vel,4,1);pos = resample(pos,4,1);
pos = pos(1:size(X,2),:); vel = vel(1:size(X,2));
inds = vel > thresh;
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
runs = bwlabel(b > 0);
w = watershed(b==0);
w = w-1; %w(w== max(w)) = 0;
allX = zeros(size(X,1)*(range(win)+1),ceil(size(X,2)/range(win)/2));counter = 1;
samplePos = [];
for k = 1:2
    runs1 = b*(-1^k)>0;runs1 = bwlabel(runs1);
    %runs1 = bwlabel(w>0 & mod(w,2) == k-1 & w <=2*max(runs));%b*((-1)^k)>0);
%     for i = 1:max(runs1)
%         sPlot([X(:,runs1 == i);angle(v(runs1 == i)/50000); inds(runs1 == i)'],[],0);pause(1);
%     end
    for i = 1:max(runs1)
        runInds = find(runs1 == i);
        %aV = -angle(v(runInds));
        %[h,pk] = findpeaks(aV,'minpeakheight',pi*.8);
        aV = angle(-v(runInds)*exp(1i*-pi/6));
        [~,pk] = findpeaks(-abs(aV),'minpeakheight',-.3);
        %h(~inds(runInds(pk))) = [];
        pk(~inds(runInds(pk))) = [];
        %plot(aV);hold all;
        %scatter(pk,h,'filled');
        %plot(vel (runInds));hold off;pause(.5);drawnow;
        pk = runInds(1)+pk-1;
        for j = 1:numel(pk)
            temp = X(:,pk(j)+(win(1):win(2)));
            allX(:,counter) = temp(:);
            counter = counter + 1;
        end
        samplePos = [samplePos; [pos(pk) i*ones(numel(pk),1)]];
%        subplot(121);plot(aV);hold all;scatter(pk,pi*ones(1,numel(pk)),'filled');hold off;title(i*10000+counter);drawnow;
%        subplot(122);plot(allX(:,max(1,counter-5):counter-1));drawnow;pause(.05);
    end
end
%[~,sortPos] = sort(samplePos,'ascend');
allX(:,counter:end) = [];
[A,W] = gfastica(allX,'lastEig',size(allX,1),'g','tanh','approach','symm','stabilization','on');
t = W*remmean(allX);
%[u,s,v] = svds(allX,10);
%figure;imagesc(allX(:,sortPos));
% Xm = reshape(mean(allX,2),[size(X,1) size(allX,1)/size(X,1)]);
% figure;imagesc(Xm);
% figure;plot(Xm');