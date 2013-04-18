function findBestICs(fields,w,thresh,fieldsConv,wConv,whiteFields)

rad = 1;%.75;
fields(:) = zscore(fields(:));
fieldsConv(:) = zscore(fieldsConv(:));
temp = squeeze(abs(mean(fields,2)));
threshInds = find(max(temp,[],2) > thresh);
fields = fields(threshInds,:,:);
w = pinv(w);
w = w(:,threshInds);
fields = bsxfun(@times,fields,exp(1i*angle(mean(w))).');
w = bsxfun(@rdivide,w,sqrt(mean(w.*conj(w))).*exp(1i*angle(mean(w))));
w = bsxfun(@minus,w,mean(w,2));

for i = 1:size(fieldsConv,1)
fieldsConv(i,:,:) = imfilter(squeeze(fieldsConv(i,:,:)),fspecial('gaussian',5,rad));
end
for i = 1:size(fields,1)
fields(i,:,:) = imfilter(squeeze(fields(i,:,:)),fspecial('gaussian',5,rad));
end
cc = corr(abs(fields(:,:))',abs(fieldsConv(:,:))');
figure;
ydim = 8;
xdim = 12;
for i =1:size(fields,1)
    ind = 4*(i-1);
    subtightplot(ydim,xdim,ind+1);imagesc(complexIm(reshape(w(:,i),[8 8])));title(i);
    subtightplot(ydim,xdim,ind+2);imagesc(complexIm(squeeze(fields(i,:,:))));
    subtightplot(ydim,xdim,ind+4);hold off;plot(1,0);set(gca,'nextPlot','add','ColorOrder',squeeze(complexIm(w(:,i))));
    [val,inds(i)] = max(cc(i,:));
    plot(squeeze(wConv(:,inds(i),:))');axis tight;
    subtightplot(ydim,xdim,ind+3);imagesc(complexIm(squeeze(fieldsConv(inds(i),:,:))));title(val);
    %useMe(i) = input('');
end
figure;
chosen = [9:11 22];
ydim = 4;xdim = 4;
for i = 1:4
    subtightplot(xdim,ydim,i);imagesc(complexIm(imfilter(squeeze(whiteFields(i,:,:)),fspecial('gaussian',5,rad))));axis off;
    hold all;plot([1 1]*size(fields,3)/2,[0 size(fields,2)],'w--','linewidth',2);
    im = squeeze(fields(chosen(i),:,:));
    %im = abs(im);
    im = complexIm(1i*im);
    subtightplot(xdim,ydim,ydim+i);imagesc(im);axis off;%complexIm());axis off;
    hold all;plot([1 1]*size(fields,3)/2,[0 size(fields,2)],'w--','linewidth',2);
    subtightplot(xdim,ydim,2*ydim+i);imagesc(complexIm(rot90(rot90(reshape(w(:,chosen(i)),[8 8])))));axis off;
    subtightplot(xdim,ydim,3*ydim+i);hold off;plot(1,0);set(gca,'nextPlot','add','ColorOrder',squeeze(complexIm(w(:,chosen(i)))),'xtick',[],'ytick',[]);
    plot(squeeze(wConv(:,inds(chosen(i)),:))','linewidth',1);axis tight;
end