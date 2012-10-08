function im = superImpC(x,frames,rad,maxVal)
%%SUPERIMP combines multiple components into 1 image, assigning each
%%component a different color, for each pixel, choosing the component with
%%the largest magnitude at that location. All components are normalized.
%%INPUTS:   c = all 2-d image components
%%          frames = which components to combine in image (choose ones with well-defined features)
%%          rad = width of gaussian smoothing kernel
%%          maxVal = if you want to normalize all components by a fixed value (default: normalize maximum of each component to 1)
% if exist('frames','var') && ~isempty(frames)
%     x = x(frames(randperm(numel(frames))),:,:);
% else
%     x = x(randperm(size(x,1)),:,:);
% end
vs = zeros(size(x,1),size(x,3));
for i = 1:size(x,1)
    [~,s,v] = svds(squeeze(x(i,:,:)),1);
    vs(i,:) = s*v';
    if exist('rad','var') && rad
        x(i,:,:) = imfilter(squeeze(x(i,:,:)),fspecial('gaussian',5,rad));
    end
    if ~exist('maxVal','var')
        x(i,:,:) = x(i,:,:)/max(max(max(abs(x(i,:,:)))));
    else
        x(i,:,:) =  x(i,:,:)/maxVal;
    end
end
[~,peakLoc] = max(abs(vs)');
x = abs(x);
x = min(1,max(0,x));
[a b]= max(x);
a = squeeze(a); b = squeeze(b);
im = zeros(3,size(a,1),size(a,2));
scale = 6;angCol = colormap('hsv');
figure;

c = mod(peakLoc*scale/size(x,3),1);
cc = angCol(ceil(c*64),:);
subplot(5,1,3);set(gca,'nextPlot','add','ColorOrder',cc,'Color',[0 0 0],'xticklabel',[],'yticklabel',[]);
plot(abs(vs)');axis tight;
subplot(5,1,[4 5]);hold on;
for i = 1:size(x,1);
    im(1,b == i) = c(i);%i/size(x,1);
    im(2,b == i) = 1;
    im(3,b == i) = a(b == i);
%     scatter(1:size(x,3),angle(vs(i,:)*exp(-1i*angle(vs(i,peakLoc(i))))),...
%         abs(vs(i,:))*100/max(abs(vs(i,:))),cc(i,:),'filled');
inds = max(1,peakLoc(i)-10):min(size(x,3),peakLoc(i)+10);
    scatter(inds,angle(vs(i,inds)*exp(-1i*angle(vs(i,peakLoc(i))))),...
        abs(vs(i,inds))*100/max(abs(vs(i,:))),cc(i,:),'filled');
%        sqrt(min(1,abs(vs(i,:))'/prctile(abs(vs(i,:)),99)))),'filled');
end
set(gca,'Color',[0 0 0]);
ylabel phase;xlabel('position (a.u.)');
axis tight;
im = max(0,im);
im = permute(im,[2 3 1]);
im = hsv2rgb(im);
subplot(5,1,[1 2]);image(im);
ylabel('trial #');