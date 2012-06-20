% gSlowp.m - slowness penalty on coeff amplitudes gaussian smoothed
%

function as = gSlowp(a,p)
filt = gausswin(p.firstlayer.a_tau_S);
filt = filt/sum(filt);
for i = 1:size(a,1)
    temp = conv(a(i,:),filt);
    temp = conv(temp(end:-1:1),filt);
    as(i,:) = temp(numel(filt):(end-numel(filt)+1));
end
norm = conv(ones(1,size(a,2)),filt);
norm = conv(norm(end:-1:1),filt);
%as = as(:,numel(filt):(end-numel(filt)+1));
norm = norm(numel(filt):(end-numel(filt)+1));
as = as./repmat(norm,[size(a,1) 1]);
as = (a - as);
%S=diff(a,1,2).^2;