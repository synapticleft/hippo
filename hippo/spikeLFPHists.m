function spikeLFPHists(cellNum,LFP,baseLine)

cells = 1:max(cellNum);
h = hist(cellNum,cells);
cells = cells(h > 5000);
xSize = ceil(sqrt(numel(cells)));
ySize = ceil(numel(cells)/xSize);
range = 15000;
if exist('baseLine','var')
    range =10000;
end
figure;
bound{1} = linspace(-range,range,50);
bound{2} = bound{1};
polar{1} = 0:200:15000;
polar{2} = -pi:.2:pi;
if ~exist('baseLine','var')
    baseLine = ones(numel(bound{1}));
end
for i = 1:numel(cells)
    temp = LFP(cellNum == cells(i));
    subplot(xSize,ySize,i);
    pic = log(hist3([real(temp),imag(temp)],bound)./baseLine);
    %pic = max(-12,pic);
    %pic = imfilter(pic,fspecial('gaussian',5,1));
    %imagesc(pic);%,[min(pic(:)) max(min(pic(:))+.1,-5)]
    pic = log(hist3([abs(temp) angle(temp)],polar));
    %imagesc(pic);
    hist(abs(temp),polar{1});
    %hist(angle(temp),-pi:.1:pi);
end