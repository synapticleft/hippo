function x = renorm(x,prc)
%% renormalize data

for i = 1:size(x,1)
    if exist('prc','var')
        x(i,:) = x(i,:) - prctile(x(i,:),prc(1));
        x(i,:) = x(i,:) / prctile(x(i,:),prc(2));
    else
        x(i,:) = x(i,:) - min(x(i,:));
        x(i,:) = x(i,:)/max(x(i,:));
    end
end