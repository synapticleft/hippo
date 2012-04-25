function vals = minSpEntArray(pos,u,s,v,spT,spId,spf,Xf,vals)
%% find electrode site that reduces the variance of the phase-position spike plots
if ~exist('Xf','var')
    Xf = u*s*v';
end
bounds = [.1 .9];
pos(pos == -1) = nan;
sV = size(Xf,2);
pos = pos(1:sV,:);
for i = 1:2
    pos(:,i) = interp1(find(~isnan(pos(:,i))),pos(~isnan(pos(:,i)),i),1:size(pos,1));
end
nanInds = isnan(pos(:,1));
pos = pos(~nanInds,:);Xf = Xf(:,~nanInds);
spf = spf(:,1:sV);
spf = spf(:,~nanInds);
nanInds = find(nanInds);
spT = spT - max(nanInds(nanInds < size(pos,1)/2));
inds = spT < 1 | spT > max(size(pos,1)); spT(inds) = [];spId(inds) = [];
pos = bsxfun(@minus,pos,mean(pos));
[pos,~,~] = svd(pos(:,1:2),'econ');
pos(:,1) = pos(:,1)-min(pos(:,1));pos(:,1) = pos(:,1)/max(pos(:,1));
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
%v(:,2) = v(:,2).*exp(1i*-angle(v(:,1)));
binsp{1} = linspace(-pi+pi/60,pi-pi/60,30); binsp{2} = linspace(bounds(1),bounds(2),50);%
if ~exist('vals','var')
figure;
vals = zeros(2,15,64);
sorted = zeros(2,size(vals,2));
for i = 1:size(vals,3)
        temp = bsxfun(@times,spf,exp(1i*-angle(Xf(i,:))));
         for k = 1:2
             runs = bwlabel(b*((-1)^k)>0);
            mem = ismember(floor(spT),find(runs>0)) | ismember(ceil(spT),find(runs>0));
            subTimes = spT(mem);subId = spId(mem);
             [~,spSort] = sort(sum(abs(spf(:,runs > 0)),2),'descend');
             sorted(k,:) = spSort(1:size(vals,2));
            posHist = hist(pos(runs > 0,1),binsp{2});
            for l = 1:size(vals,2)%size(spf,1)
                tempTimes = subTimes(subId == spSort(l));%tempTimes = tempTimes(ismember(floor(tempTimes),find(runs>0)) | ismember(ceil(tempTimes),find(runs>0)));
                hp = hist3([angle(weighted(tempTimes,temp(spSort(l),:).')) weighted(tempTimes,pos(:,1))],binsp);
                hp = bsxfun(@rdivide,hp,posHist);
                vals(k,l,i) = ent(hp);
                %subplot(2,2,(k-1)*2+l);imagesc(hp);title(num2str(vals(k,l,i,j)));drawnow;
            end
         end
end
else
    %vals = vals(:,1:10,:);
    numDisp = 3;
    %bins{1} = linspace(-pi,pi,30);bins{2} = bins{1};%
    for k = 1:2
        figure;
            runs = bwlabel(b*((-1)^k)>0);
            posHist = hist(pos(runs > 0,1),binsp{2});
            mem = ismember(floor(spT),find(runs>0)) | ismember(ceil(spT),find(runs>0));
            subTimes = spT(mem);subId = spId(mem);
             [~,spSort] = sort(sum(abs(spf(:,runs > 0)),2),'descend');
            % vDemod = v(:,2).*exp(1i*-angle(v(:,1)));
            %v2Hist = hist(angle(vDemod(runs > 0)),bins{2});
        for i = 1:size(vals,2)
            subplot(size(vals,2),numDisp,numDisp*(i-1)+1);imagesc(reshape(squeeze(vals(k,i,:)),[8 8]));ylabel(num2str(spSort(i)));set(gca,'ytick',[],'xtick',[]);
            %smooth = imfilter(squeeze(vals(k,i,:,:)),fspecial('gaussian',3,1),'replicate');
            %subplot(size(vals,2),numDisp,numDisp*(i-1)+2);imagesc(smooth);axis off;
            [~,minInd] = min(squeeze(vals(k,i,:)));
            %[x y] = find(smooth == min(smooth(:)));
            %alpha = (range(x) + 1i*range(y));%title(num2str(alpha));
            tempTimes = subTimes(subId == spSort(i));%tempTimes = tempTimes(ismember(floor(tempTimes),find(runs>0)) | ismember(ceil(tempTimes),find(runs>0)));
            temp = spf(spSort(i),:).*exp(1i*-angle(Xf(minInd,:)));%v(:,1).'));
            hp = hist3([angle(weighted(tempTimes,temp.')) weighted(tempTimes,pos(:,1))],binsp);
            hp = bsxfun(@rdivide,hp,posHist);minMax = [min(hp(:)) max(hp(:))];
            subplot(size(vals,2),numDisp,numDisp*(i-1)+2);imagesc(log(hp));axis off;%title(num2str(ent(hp)));
            [~,maxInd] = max(squeeze(vals(k,i,:)));
            temp = spf(spSort(i),:).*exp(1i*-angle(Xf(maxInd,:)));%v(:,1)+alpha*v(:,2)).');
            hp = hist3([angle(weighted(tempTimes,temp.')) weighted(tempTimes,pos(:,1))],binsp);
            hp = bsxfun(@rdivide,hp,posHist);
            subplot(size(vals,2),numDisp,numDisp*(i-1)+3);imagesc(log(hp),log(minMax));axis off;%title(num2str(ent(hp)));
            %h = hist3([angle(weighted(tempTimes,v(:,1))) angle(weighted(tempTimes,vDemod))],bins);
            %h = bsxfun(@rdivide,h,v2Hist);
            %subplot(size(vals,2),numDisp,numDisp*(i-1)+5);imagesc(sqrt(h));axis off;
%             vTemp = ([1 alpha;-conj(alpha) 1]*v.').';
%             vTemp(:,2) = vTemp(:,2).*exp(1i*-angle(vTemp(:,1)));
%             v2TempHist = hist(angle(vTemp(runs > 0,2)),bins{2});%figure;plot(v2TempHist);hold all;plot(v2Hist);return
%             h = hist3([angle(weighted(tempTimes,vTemp(:,1))) angle(weighted(tempTimes,vTemp(:,2)))],bins);
%            % h = bsxfun(@rdivide,h,v2TempHist);
%             temp = abs(u(:,2)./u(:,1) - alpha);
%             subplot(size(vals,2),numDisp,numDisp*(i-1)+6);imagesc(reshape(temp,[8 8]),[0 max(temp)]);axis off;
%             if exist('sh','var')
%                 hold on;scatter(sh(spSort(i),2),4,'filled');
%             end
         end
     end
end

% Xf = u*v';Xf = Xf(:,1:end-1);
% Xf1 = Xf(:,1:end-1).*conj(Xf(:,2:end))./abs(Xf(:,2:end));
% [a b1] = find(vals == min(vals(:)));
% val = range(a)+1i*range(b1)
% figure;imagesc(range,range,vals);
% figure;subplot(311);imagesc(reshape(std(Xf1,0,2)./std(Xf,0,2),[8 8]));
% subplot(312);imagesc(reshape(std(Xf1(:,b>0),0,2)./std(Xf(:,b>0),0,2),[8 8]));
% subplot(313);imagesc(reshape(abs(u(:,2)./u(:,1) - val),[8 8]));

function v = ent(in)
in = in(:)/sum(in(:));
v = -nansum(in.*log2(in));
%v(isnan(v)) = 0;

function out = weighted(inds,vals)
w = mod(inds,1);
out = (1-w).*vals(floor(inds))+w.*vals(ceil(inds));