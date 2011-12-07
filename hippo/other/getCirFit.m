function [a b] = getCirFit(cont,inds)

figure;hold on;
contour(10:10:80,1:8,cont,inds,'LineWidth',2);
%[x y] = meshgrid(1:8,10:10:80);
cont = contourc(10:10:80,1:8,cont,inds);

layer = giveConts(cont);
%layer{4}.vals((end-4):end,:) = [];
for i = 1:length(layer)
    %plot(layer{i}.vals(:,1),layer{i}.vals(:,2)); hold all;
end

col = colormap;
ind = colInd(numel(inds));

figure;hold all;
thetas = -.05:.05:(2*pi);
for i = 1:length(layer)
    scatter(layer{i}.vals(:,1),layer{i}.vals(:,2),'filled','CData',col(ind(i),:));
    [a{i} b{i} c d e] = ls2dcircle(layer{i}.vals,[30 -32]',30,1,1);
    %scatter(a{i}(1),a{i}(2),'filled');
    plot((a{i}(1)+b{i}*cos(thetas)),(a{i}(2) + b{i}*sin(thetas)),'color',col(ind(i),:),'LineWidth',2);
end

function ind = colInd(numInd)
nCol = 63;
inds = (1:numInd)-1;
inds = floor(inds/(numInd-1)*nCol);
ind = inds + 1;


function layer = giveConts(cont)
state = 1;
for i = 1:Inf
    layers{i}.height = cont(1,state);
    numElems = cont(2,state);
    layers{i}.vals = cont(:,(state+1):(state+numElems))';
    %layers{i}.vals(:,1) = layers{i}.vals(:,1) * 10;
    state = state + numElems + 1;
    if state > size(cont,2)
        break
    end
end
state = 1;
i = 1;
while i <= length(layers)
    layer{state}.height = layers{i}.height;
    layer{state}.vals = [];
    while layer{state}.height == layers{i}.height
        layer{state}.vals = [layer{state}.vals; layers{i}.vals];
        i = i+1;
        if i > length(layers)
            break
        end
    end
    state = state + 1;
end