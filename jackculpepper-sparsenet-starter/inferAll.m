% for i = 1:2
%     posd(:,i) = pos(:,i) - min(pos(:,i)) + eps;
% end
% posd = ceil(posd);
% junk = any(isnan(posd'));
% posd(junk,:) = [];
%Xf(:,junk) = [];
posAccum = zeros(M,max(posd(:,1)),max(posd(:,2)));
myInds = meshgrid(1:M,1:Btest);
for i = 1:floor(size(Xf,2)/Btest)
    %close all;
    atest1 = zeros(M,Btest);
    inds = (i-1)*Btest+(1:Btest);
    atest1 = lbfgs(@objfun_a,atest1(:),lb,ub,nb,opts,phi,Xf(:,inds),lambda);
    atest1 = reshape(atest1,M,numel(atest1)/M)';
    imagesc(sqrt(abs(atest1')));drawnow;
    posAccum = posAccum + accumarray([myInds(:) repmat(posd(inds,:),[M 1])],(atest1(:)),[M max(posd)],@sum);
end