function X = crop_chunk(F,m,p)

X = F;

% if p.normalize_crop
%     X = X-mean(X(:));
%     X=X/(sqrt(10*var(X(:))));
% end

if p.whiten_patches && isfield(m,'whitenMatrix')
%    X = bsxfun(@minus,X,m.imageMean);
    X = m.whitenMatrix*X;
end

if ~p.whiten_patches
%    X = bsxfun(@minus,X,m.imageMean);
X = X/std(X(:));
%    X = bsxfun(@rdivide,X,m.imageStd);
    %X = X/p.var;
end