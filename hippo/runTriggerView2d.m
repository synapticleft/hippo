function [actMean] = runTriggerView2d(pos,Xf,binSize,v,r) %,v actMean,
%% binning of 2-d data in open maze

if ~exist('v','var') || isempty(v)
    v = Xf(1,:)';
end
if size(v,1) < size(pos,1)
   pos = pos(1:size(v,1),:);
end
% pos = fixPos(pos); %removed for human work
% pos(pos == -1) = nan;
% reject = 1;
% while sum(reject)
% reject = 0;
% for i = 1:4
%     reject = reject | min([0; diff(pos(:,i))],flipud([0; diff(flipud(pos(:,i)))])) < -20;
% end
% pos(reject,:) = nan;
% for i = 1:4
%     nanInds = find(~isnan(pos(:,i)));
%     pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
% end
% end
nanInds = isnan(pos(:,1)) | isnan(pos(:,2));
pos = pos(~nanInds,:);Xf = Xf(:,~nanInds);v = v(~nanInds,:);%sp = sp(:,~nanInds);
%vel = angVel(pos);
%vel = [zeros(1,2); vel(:,1:2)];
%for i = 1:2
%    vel(:,i) = filtLow(vel(:,i),1250/32,1);
%end
%veld = [ vel(:,1:2)];
%vel = vel(:,1);
%vel = vel/max(vel);inds = vel > thresh;
%pos = bsxfun(@minus,pos,mean(pos));
%[a,~,~] = svd(pos(:,1:2),'econ');pos = a;
pos = pos(:,1:2);
for i = 1:2    
    pos(:,i) = pos(:,i) - min(pos(:,i))+eps;
%    pos(:,i) = pos(:,i)/(max(pos(:,i)));
%    pos(:,i) = min(pos(:,i),.9999);
%     veld(:,i) = veld(:,i) - min(veld(:,i));
%     veld(:,i) = veld(:,i)/max(veld(:,i));
%     veld(:,i) = min(veld(:,i),.9999);
     posd(:,i) = ceil(pos(:,i)/binSize);%floor(pos(:,i)*accumbins(min(numel(accumbins),i)))+1;
%     veld(:,i) = floor(veld(:,i)*accumbins(min(numel(accumbins),i)))+1;
end
Xf = bsxfun(@times,Xf,exp(1i*angle(v(:,1))).');
%inds = bwmorph(inds,'dilate',20);
%Xf = Xf(:,inds);posd = posd(inds,:);%veld = veld(inds,:);
%vel = vel(inds);%pos = pos(inds,:);
%if exist('posInds','var') && ~isempty(posInds)
%    r = r(:,posInds); %% IS THIS RIGHT??
%end
if exist('r','var')
    Xf = r*Xf;
end
%Xf = bsxfun(@times,Xf,exp(1i*angle(v(:,1))).');
actMean = zeros(size(Xf,1),max(posd(:,1)),max(posd(:,2)));
% = absMean;
for i = 1:size(Xf,1)
    actMean(i,:,:) = accumarray(posd,Xf(i,:),[],@mean);
%    absMean(i,:,:) = accumarray(posd,abs(Xf(i,:)),[],@mean);
end
%complexAct = 0;
%if ~complexAct
%    t = Xf;%r*zscore(Xf,0,2);
%    clear Xf;
%    sk = sign(skewness(r));%
%%    t = bsxfun(@times,t,sk');
%    r = bsxfun(@times,r,sk);
%else
%    Xf = zscore(Xf,0,2);Xf = complex(Xf(1:end/2,:),Xf(end/2+1:end,:));
%    t = complex(r(:,1:end/2),r(:,end/2+1:end))*conj(Xf);
%end
% %2d stuff
%if ~exist('posInds','var') || isempty(posInds)
%    posInds = 1:size(r,1);
%end
%xdim = ceil(sqrt(numel(posInds)));ydim = ceil(numel(posInds)/xdim);
%f1 = figure;f2 = figure;
%t = bsxfun(@times,t, sign(skewness(t,0,2)));
%t = zscore(t,0,2);
% for i = 1:numel(posInds)
%     u = r(:,i);%complex(r1(1:size(Xf,1)/2-1,i),r1(size(Xf,1)/2+1:end-1,i));%r1(1:size(Xf,1)-1,posInds(i));%
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
%         up1 = reshape(u,[8 8]);
%     end
%     figure(f1);subplot(xdim,ydim,i);imagesc(complexIm(up1,0,1));axis off;
%     temp = accumarray(posd,t(i,:),accumbins,@mean,nan);
%     temp(isnan(temp)) = prctile(temp(:),20);
%     figure(f2);subplot(xdim,ydim,i);imagesc(imfilter(temp,fspecial('gaussian',5,1),'replicate'));axis off;
% end
% sPlot([10*vel';t;abs(v(inds,1)')/1000]);