function X = crop_chunk(F,m,p)

if strcmp(p.data_type,'lfp')
    X = F;
elseif strcmp(p.data_type,'sim')
    X = reshape(F,size(F,1)*size(F,2),size(F,3));
end

% if p.normalize_crop
%     X = X-mean(X(:));
%     X=X/(sqrt(10*var(X(:))));
% end

if p.whiten_patches && isfield(m,'whitenMatrix')
    X = bsxfun(@minus,X,m.imageMean);
    X = m.whitenMatrix*X;
end

if ~p.whiten_patches
    X = bsxfun(@minus,X,m.imageMean);
    X = bsxfun(@rdivide,X,m.imageStd);
    %X = X/p.var;
end