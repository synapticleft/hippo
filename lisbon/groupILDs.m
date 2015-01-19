function inds = groupILDs(mat,cILD,dynamic)

mat = squeeze(mat(:,:,end));
if nargin < 3
    dynamic = 0;
end
mat = mat(:,sum(mat == 0) == 0);

mat = mat-cILD;
neg = mat(:,1) < 0;
if dynamic
    mat(neg,:) = -mat(neg,:);
    %mat = bsxfun(@minus,mat,mat(:,1));
end

dists = squareform(pdist(mat));

figure;imagesc(dists < 0.01);

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
figure;imagesc(mat(s,:));