sst = sst(:,:)';
sstm = ~isnan(sst(:,1));

for i = 1:648
    inds = (i-1)*100+(1:100);
    sst(inds,:) = morFilter(sst(inds,:),10,120);
    if ~mod(i,50)
        i
    end
end

sst(sstm,:) = [];
sst = bsxfun(@times,sst,exp(-1i*(1:1642)*2*pi/12));
sstAngle = angle(mean(sst,2));
sst = bsxfun(@times,sst,exp(-1i*sstAngle));
%%do complexICA