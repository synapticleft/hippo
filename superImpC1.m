function im = superImpC1(x,rad,thresh,W) %allxc
%%SUPERIMPC combines multiple components into 1 image, assigning each
%%component a different color, for each pixel, choosing the component with
%%the largest magnitude at that location. All components are normalized.
%%It is for complex numbers
%%INPUTS:   c = all 2-d image components
%%          rad = width of gaussian smoothing kernel
%%          maxVal = if you want to normalize all components by a fixed value (default: normalize maximum of each component to 1)
warning off all;
x(:) = zscore(x(:));
%xm = squeeze(abs(mean(x,2)));
%x = x(max(xm,[],2)>thresh,:,:);
x = permute(x,[1 3 2]);
buffer = 24;
f = fspecial('gaussian',5,rad);
tempF = zeros(size(x,1),size(x,2));
for i = 1:size(x,1)
        %filts{i} = [];
        tempFull = squeeze(x(i,:,:));
        %%tempFull = [zeros(size(tempFull,1),buffer) tempFull zeros(size(tempFull,1),buffer)];
        tempFullF = filtfilt(gausswin(12),sum(gausswin(12)),tempFull);
        tempF(i,:) = mean(abs(tempFullF.'));
        [~,locs] = findpeaks(tempF(i,:),'minpeakheight',thresh,'minpeakdistance',buffer/2);
        inds = [];
        for j = 1:numel(locs)
            inds = [inds max(1,locs(j)-buffer/2):min(size(x,2),locs(j)+buffer/2)];
        end
        x(i,~ismember(1:size(x,2),inds),:) = 0;
        tempF(i,~ismember(1:size(x,2),inds)) = 0;
        if exist('rad','var') && rad
             x(i,:,:) = imfilter(squeeze(x(i,:,:)),f);
        end
end
remInds = sum(tempF,2)==0;
%A = pinv(W);
%W(remInds,:) = [];A(:,remInds) = [];
x(remInds,:,:) = [];
tempF(remInds,:) =[];
[a b]= max(abs(x));
upperBound = prctile(abs(a(:)),90);
vs = squeeze(mean(x,3))/upperBound;
a = squeeze(a); b = squeeze(b);
a = max(0,min(1,a/upperBound));
im = zeros(3,size(a,1),size(a,2));
scale = 1;angCol = colormap('hsv');
[~,peakLoc] = max(tempF');
[~,so] = sort(peakLoc,'ascend');
bound = sum(peakLoc < size(x,2)/2);
figure;
c = mod(peakLoc*scale/size(x,2),1);%rand(1,numel(peakLoc));%
cc = angCol(max(1,ceil(c*64)),:);
subtightplot(3,1,2);set(gca,'nextPlot','add','ColorOrder',cc,'Color',[0 0 0],'xticklabel',[],'yticklabel',[],'fontsize',16);
plot(abs(vs)','linewidth',2);hold all;
plot([size(x,2)/2 size(x,2)/2]+1,[0 max(abs(vs(:)))*1.01],'w','linewidth',5);
axis tight;set(gca,'ylim',[0 1]);
subtightplot(3,1,3);hold on;
for i = 1:size(x,1);
    im(1,b == i) = c(i);%i/size(x,1);
    im(2,b == i) = 1;
    im(3,b == i) = a(b == i);
    inds = find(tempF(i,:) > 0);
    s =  abs(vs(i,inds))*100;
    scatter(inds,angle(vs(i,inds)),s,cc(i,:),'filled');%,'ytick',[-3 0 3],'yticklabel',{'-pi','0','pi'}
%        sqrt(min(1,abs(vs(i,:))'/prctile(abs(vs(i,:)),99)))),'filled');
end
plot([size(x,2)/2 size(x,2)/2]+1,[-pi pi],'w','linewidth',5);
set(gca,'Color',[0 0 0],'xtick',[1 size(x,2)/2 size(x,2)],'xticklabel',[0 250 0],'fontsize',16,...
    'xlim',[1 size(x,2)],'ylim',[-pi pi]);
ylabel phase;xlabel('position (cm)');
im = max(0,im);
im = permute(im,[3 2 1]);
im = hsv2rgb(im);
subtightplot(3,1,1);image(im);set(gca,'xtick',[],'fontsize',16);hold all;
plot([size(x,2)/2 size(x,2)/2]+1,[1 size(im,1)],'w','linewidth',5);
ylabel('trial #');
% for i = 1:size(x,1)
%     temp = circshift(x(i,:,:),[0 0 1]);
%     xc(i,i) = corr((temp(:)),(x(i,:)).');
% end

%x = double(abs(x(so,:)) > thresh);
%xc = x*x'/size(x,2);%abs(corr((x(:,:))'));
%figure;imagesc(xc,[0 .1]);hold all;plot([.5 size(xc,1)+.5],[bound bound]+.5,'w','linewidth',3);
%plot([bound bound]+.5,[.5 size(xc,1)+.5],'w','linewidth',3);axis image;

%Ac = abs(corr(A(:,so)));
%Wc = abs(corr(W(so,:)'));
%for i = 1:size(x,1)
%    xc(i,i) = -1;
%end
% figure;subplot(121);imagesc(Ac);axis image;
% subplot(122);imagesc(Wc);axis image;
% %figure;scatter(xc(:),Ac(:),'filled');hold all;scatter(xc(:),Wc(:),'r','filled');
% figure;subplot(121);bar([mean(Ac(xc == 0)) mean(Ac(xc > 0 & xc < .01)) mean(Ac(xc > .01))]);
% subplot(122);bar([mean(Wc(xc == 0)) mean(Wc(xc > 0 & xc < .01)) mean(Wc(xc > .01))]);
% figure;
% for i= 1:5
%     subplot(121);scatter(i,mean(diag(Ac,i)),'filled');hold all
%     subplot(122);scatter(i,mean(diag(Wc,i)),'filled');hold all;
% end
%    axis tight; 
%buff = 5;
%allxc = [];
%for i = -buff:buff
%    dx = diag(xc,i);
%    allxc = [allxc; dx ones(size(dx))*abs(i)+1];
%end
%figure;boxplot(allxc(:,1),allxc(:,2),'plotstyle','compact');
%r = pinv(r).';
%if exist('r','var')
%    r(remInds,:) = [];
%    xcr = abs(corr(r'));
%    figure;scatter(xcr(:),xc(:));
%end