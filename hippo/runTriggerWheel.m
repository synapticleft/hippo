function [posInds,t,tes] = runTriggerWheel(pos,v,Xf,thresh,r,probes,posInds,r2)

%X = getData('ec014.440.h5',97:98,[],[],1);
%whl = decimate(abs(diff(X(1,:))),32);
%whl(2,:) = decimate(abs(diff(X(2,:))),32);
%whl = sum(whl);

ran = round(1250/32)*[-1 15];
if size(v,1) < size(pos,1)
    pos = pos(1:size(v,1),:);
end
% for i = 1:size(vel,2)
%     vel(:,i) = filtLow(vel(:,i),1250/32,1);
% end
inds = bwlabel(pos > thresh);
%inds = bwmorph(inds,'dilate',20);
h = hist(inds,1:max(inds));
inds(ismember(inds,find(h < 1250/32*5))) = 0;
%inds = bwlabel(inds > 0);
indStart = find(diff(inds) > 0);
indEnd = find(diff(inds) < 0);
%figure;plot((1:numel(pos))/1250*32,pos/max(pos));hold all;plot((1:numel(pos))/1250*32,inds>0);return
Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).')];
%Xf = [real(Xf);imag(Xf)];
r1 = pinv(r);
if ~exist('r2','var')
    r2 = r1;
end
if exist('posInds','var') && ~isempty(posInds)
    r = r(posInds,:);r1 = r1(:,posInds);r2 = r2(:,posInds); %% IS THIS RIGHT??
end
Xf = r*Xf;%zscore(Xf,0,2);%
velPlots = zeros(2,numel(indStart),range(ran));
icPlots = zeros(2,size(r,1),numel(indStart),range(ran));
icTime = zeros(size(r,1),numel(indStart),50);
for i = 1:numel(indStart)
    for j = 1:2
        if j == 1
            in = indStart(i)+ran(1):indStart(i)+ran(2)-1;
        else
            in = indEnd(i) -ran(2):indEnd(i) - ran(1) -1;
        end
        velPlots(j,i,:) = pos(in);
        icPlots(j,:,i,:) = Xf(:,in);
    end
    in = indStart(i):indEnd(i);
    for j = 1:size(r,1)
        icTime(j,i,:) = interp1(1:numel(in),Xf(j,in),linspace(1,numel(in),size(icTime,3)));
    end
end
figure;subplot(211);imagesc(squeeze(velPlots(1,:,:)));
subplot(212);imagesc(squeeze(velPlots(2,:,:)));
h1 = figure;h2 = figure;h3 = figure;h4 = figure;
xdim = ceil(sqrt(size(r,1)));ydim = ceil(size(r,1)/xdim);
profile = zeros(2,size(r,1),range(ran));
for i = 1:size(r,1)
    im = imfilter(squeeze(icPlots(1,i,:,:)),fspecial('gaussian',10,2));
    [~,~,profile(1,i,:)] = svds(im,1);
    figure(h1);subplot(xdim,ydim,i);imagesc(complexIm(im*(mean(r1(:,i))),0,1));axis off;
    im = imfilter(squeeze(icPlots(2,i,:,:)),fspecial('gaussian',10,2));
    [~,~,profile(2,i,:)] = svds(im,1);
    figure(h3);subplot(xdim,ydim,i);imagesc(complexIm(im*(mean(r1(:,i))),0,1));axis off;
    im = imfilter(squeeze(icTime(i,:,:)),fspecial('gaussian',5,1));
    figure(h4);subplot(xdim,ydim,i);imagesc(complexIm(im*mean(r1(:,i)),0,1));axis off;
    figure(h2);subplot(xdim,ydim,i);imagesc(complexIm(reshape(r2(:,i)*conj(mean(r1(:,i))),[8 size(r2,1)/8]),0,1));axis off;
end
figure;subplot(211);plot(abs(squeeze(profile(1,:,:)))');
subplot(212);plot(abs(squeeze(profile(2,:,:)))');
icPlots = squeeze(icPlots(1,:,:,:));
icPlots = bsxfun(@times,icPlots,exp(1i*angle(mean(r1))).');
superImpC(icPlots,[],2,prctile(abs(icPlots(:)),99.5));
return

% % %2d stuff
% %
% % for i = 1:size(t,1)
% %    cc(i,:) = xcorr(t(i,:),vel,1000);
% % end
% % sPlot(cc);
% %figure;plot(cc');
% %xdim = ceil(sqrt(size(cc,1)));ydim= ceil(size(cc,1)/xdim);
% %posd = posd(inds,:);veld = veld(inds,:);
% %figure;for i = 1:size(cc,1)
% %    subplot(xdim,ydim,i);imagesc(imfilter(accumarray(veld,t(i,:),accumbins,@mean,0),fspecial('gaussian',5,1)));
% %end
% if ~exist('posInds','var') || isempty(posInds)
%     posInds = 1:size(r1,2);
% end
% xdim = ceil(sqrt(numel(posInds)));ydim = ceil(numel(posInds)/xdim);
% f1 = figure;f2 = figure;
% %[sk,si] = sort(abs(skewness(t,0,2)),'descend');
% %r1 = r1(:,si);t = t(si,:);
% tes = zeros(numel(posInds),accumbins(1),accumbins(2));
% for i = 1:numel(posInds)
%     u = r2(:,i);%r1(1:size(Xf,1)-1,posInds(i));%
%     if exist('probes','var') && ~isempty(probes)
%         up1 = probes;
%         for ii = 1:size(probes,1)
%             for j = 1:size(probes,2)
%                 up1(ii,j) = u(probes(ii,j)+1);%-256
%             end
%         end
%         %    up1 = up1(:,[1:4 6 5 8 7]);
%         up1 = up1(:,[1:12 14 13 16 15]);
%         %up1 = diff(up1);
%         up1 = [up1(:,1:8) zeros(size(up1,1),1) up1(:,9:16)];
%     else
%         up1 = reshape(u*exp(-1i*angle(mean(r1(:,i)))),[8 size(Xf,1)/8]);%[16 6]);%
%     end
%     figure(f1);subplot(xdim,ydim,i);imagesc(complexIm(up1,0,1));axis off;
%     ac = imfilter(accumarray(posd,t(i,:),accumbins,@mean,0),fspecial('gaussian',5,1),'replicate');
%     tes(i,:,:) = ac;
%     figure(f2);subplot(xdim,ydim,i);imagesc(complexIm(ac,0,1,[],[]));axis off;
%     %title(kurtosis(abs(ac(:))))
% end
% [s,f] = sort(max(abs(tes(:,:)),[],2),'descend');
% figure;plot(s);
% s1 = input('cutoff low?: ');
% s2 = input('cutoff high?: ');
% f = f(s1:s2);
% %f = find(max(abs(tes(:,:)),[],2) > 2)
% superImpC(tes,f,1,prctile(abs(tes(:)),99.9));
% sPlot(abs(t));%[10*vel';t]);
% return
% %%FOR 1D TRACK
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
runs = bwlabel(b > 0);
w = watershed(b==0);
w = w-1; 
posd(mod(w,2) ==1 ,1) = posd(mod(w,2) ==1 ,1) + max(posd(:));
pos(mod(w,2) ==1 ,1) = pos(mod(w,2) ==1 ,1) + max(pos(:));
runs1 = round(w/2);
inds = runs1 > 0 & runs1 <= max(runs);
t1 = zeros(size(t,1),max(runs),accumbins(1)*2);
for j = 1:size(t,1)
         t1(j,:,:) = accumarray([runs1(inds); posd(inds,1)']',t(j,inds),[max(runs) 2*accumbins(1)] ,@mean);
end
%h1 = figure;
spatial = randn(size(t1,1),2*accumbins(1));
for i = 1:size(t1,1)
    temp = reshape(t1(i,:),[max(runs) 2*accumbins(1)]);
    [u,s,v1] = svds(temp,1);
    spatial(i,:) = exp(1i*angle(mean(u)))*s*v1';% 
%     if -min(v) > max(v) 
%         v = -v; 
%     end
end

if ~exist('posInds','var') || isempty(posInds)
    posInds = find(max(abs(spatial')) > 10);
else
   posInds = 1:size(r,1);%
end

t = reshape(t(:),size(t));
t1 = reshape(t1(:),size(t1));
figure;plot(abs(spatial)');
spatial = spatial(posInds,:);
[~,peakLoc] = max(abs(spatial)');
[~,indLoc] = sort(peakLoc);
%peakLoc = peakLoc(indLoc);
posInds = posInds(indLoc);
spatial = spatial(indLoc,:);
t = t(posInds,:);t1 = t1(posInds,:,:);
t = bsxfun(@times,t,exp(1i*-angle(v(:,1))).');
% sPlot(t);
% sPlot(morFilter(t,1250/32,8));
% sPlot(abs(t));
h1 = figure;
h2 = figure;
xdim = ceil(sqrt(numel(posInds)));ydim = ceil(numel(posInds)/xdim);
sk = ones(1,numel(posInds));
tes = zeros(numel(posInds),max(runs),2*accumbins(1));
if exist('probes','var') && ~isempty(probes)
    ups = zeros(numel(posInds),size(probes,1),size(probes,2)+1);
else
    ups = zeros(numel(posInds),8,(size(Xf,1))/8);%/2
end

% meanAng = zeros(1,size(t,1));meanAng1 = meanAng;
% for i = 1:size(t,1)
%     meanAng(i) = angle(mean(r1(1:size(Xf,1),posInds(i))));%
%     meanAng1(i) = angle(mean(t(i,posd(:,1) == peakLoc(i))));% < peakLoc(i)+accumbins(1)/10 & posd(:,1) > peakLoc(i)-accumbins(1)/10));
%     %meanAng(i) = angle(meanAng(i));
% %    meanAng(i) = circ_mean(angle(r1(1:size(Xf,1)-1,posInds(i))),abs(r1(1:size(Xf,1)-1,posInds(i))));
% end
% figure;scatter(meanAng,meanAng1);return
% spatial = bsxfun(@times,spatial,exp(1i*(pi/2-meanAng')));
%figure;imagesc(complexIm(spatial,0,1));
%spsh = spatial;
%for i = 1:size(spsh,1)
%    spsh(i,:) = circshift(spsh(i,:),[0 -peakLoc(i)+accumbins(1)]);
%end
%figure;imagesc(complexIm(spsh,0,1));
for i = 1:numel(posInds)
    te = reshape(t1(i,:),[max(runs) 2*accumbins(1)]);
%     if skewness(te(:)) < 0
%         te = -te;
%         sk(i) = -1;
%     end
    u = r1(:,posInds(i));%*exp(1i*(pi/2-meanAng(i)));
    tes(i,:,:) = te*exp(1i*angle(mean(u)));
    figure(h1);subplot(xdim,ydim,i);imagesc(complexIm(imfilter(squeeze(tes(i,:,:)),fspecial('gaussian',5,1)),0,1));axis off;%s(temp(indLoc(i))),[0 max(te(:))].*exp(1i*(pi/2-meanAng(i)))
    %meanAng1(i) = ;
    u = r2(:,posInds(i))*exp(-1i*angle(mean(u)));
    if exist('probes','var') && ~isempty(probes)
        up1 = probes;
        for ii = 1:size(probes,1)
            for j = 1:size(probes,2)
                up1(ii,j) = u(probes(ii,j)+1-min(probes(:)));%-256
            end
        end
        %up1 = up1(:,[1:4 6 5 8 7]);
        up1 = up1(:,[1:12 14 13 16 15]);
        up1 = [up1(:,1:8) zeros(size(up1,1),1) up1(:,9:16)];
    else
        up1 = reshape(u,[8,(size(Xf,1))/8]);
    end
%    ups(i,:,:) = up1;
    figure(h2);subplot(xdim,ydim,i);imagesc(complexIm(up1,0,1));axis off;
end
% h3 = figure;subplot(311);imagesc(reshape(std(r),[8 size(Xf,1)/8]));
% subplot(312);imagesc(reshape(std(r(posInds,:)),[8 size(Xf,1)/8]));
% subplot(313);imagesc(reshape(mean(abs(r(posInds,:))),[8 size(Xf,1)/8]));
% freezeColors(h3);
%sPlot([bsxfun(@times,t,sk'); vel']);
%figure;imagesc(complexIm(corr(ups(:,:)'),0,1));
superImpC(tes,[],1,prctile(abs(tes(:)),99.5));

% figure;
% for i = 1:max(runs1)
%     for j = 1:size(t,1)
%         indUse = runs1' == i & posd(:,1) < peakLoc(j)+accumbins(1)/10 & posd(:,1) > peakLoc(j)-accumbins(1)/10;
% %        plot(pos(indUse,1),posd(indUse,1));
% %        meanAng(j) = angle(mean(t(j,indUse)));
%         scatter(pos(indUse,1),angle(t(j,indUse)*exp(-1i*meanAng(j))),abs(t(j,indUse))*40,'filled');%angle(spatial(i,peakLoc(j)))
%  %       plot(pos(indUse,1),angle(t(j,indUse)*exp(-1i*meanAng(j))),'k');%angle(spatial(i,peakLoc(j)))
%         hold all;
%     end
%     hold off;pause(2);
% end
