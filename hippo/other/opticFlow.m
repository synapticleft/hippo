function opticFlow(phi)

NumFrames = size(phi,3);
alpha=.01;
iterations=20;
% Width and Height
U=cell(NumFrames-1,1);
V=cell(NumFrames-1,1);
Normal=cell(NumFrames-1,1);
hsize=5;
sigma=1;
h = fspecial('gaussian',hsize,sigma);
height = 8;width = 8;
for i = 1:size(phi,2)
    phiT = squeeze(phi(:,i,:));
    PreviousFrame = reshape(phiT(:,1),[height width]);  
    PreviousFrame= imfilter(PreviousFrame,h,'replicate');
    for j = 2:NumFrames
        I= reshape(phiT(:,j),[height width]);
        CurrentFrame=imfilter(I,h,'replicate');
        [TempU TempV TempNormal]= computeFlow(CurrentFrame,PreviousFrame,alpha,iterations);
        U{j-1}= TempU ;
        V{j-1} = TempV;
        Normal{j-1}= TempNormal;
        PreviousFrame=CurrentFrame;
        plotFlow(TempU,TempV,PreviousFrame);
    end
end

