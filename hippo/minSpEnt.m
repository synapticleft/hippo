function [vals sorted] = minSpEnt(pos,u,s,v,spT,spId,spf)
%% find static superposition of 1st 2 activations that reduces the variance of the demodulated activation
figure;
v = v*s;
bounds = [.1 .9];
pos(pos == -1) = nan;
sV = size(v,1);
pos = pos(1:sV,:);
for i = 1:2
    pos(:,i) = interp1(find(~isnan(pos(:,i))),pos(~isnan(pos(:,i)),i),1:size(pos,1));
end
nanInds = isnan(pos(:,1));
pos = pos(~nanInds,:);v = v(~nanInds,:);
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
vals = zeros(2,2,10,10);
range = linspace(-3,3,size(vals,3));
binsp{1} = linspace(-pi,pi,30); binsp{2} = linspace(bounds(1),bounds(2),50);%
sorted = zeros(2,size(vals,2));
for i = 1:size(vals,3)%range
    for j = 1:size(vals,4)
        alpha = range(i)+1i*range(j);
        vNew = v(:,1) + alpha*v(:,2);
        temp = bsxfun(@times,spf,exp(1i*angle(vNew)).');
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
                vals(k,l,i,j) = ent(hp);
                subplot(2,2,(k-1)*2+l);imagesc(hp);title(num2str(vals(k,l,i,j)));drawnow;
            end
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