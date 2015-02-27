function inds = groupILDs(mat,cILD,dynamic,figsOn)

if nargin < 4
    figsOn = 0;
end

mat = squeeze(mat(:,:,end));
if nargin < 3
    dynamic = 0;
end
mat = mat(:,sum(mat == 0) == 0);

mat = mat-cILD;
neg = mean(mat,2) < 0;%mat(:,1) < 0;
if dynamic
    mat(neg,:) = -mat(neg,:);
    %mat = bsxfun(@minus,mat,mat(:,1));
end

dists = squareform(pdist(mat));

counter = 1;
inds = zeros(1,size(mat,1));
for i = 1:size(mat,1)
    if ~inds(i)
        inds(dists(i,:) < .01) = counter;
        counter = counter + 1;
    end
end
if dynamic 
    inds(neg) = inds(neg) + .5;%-inds(neg);
end
inds = inds *2-1;
[~,s] = sort(inds);
if figsOn
figure;imagesc(dists < 0.01);
figure;imagesc(mat(s,:));
end