function sPlot(data,xaxis,newFig,norm)
%spaced out plot of 2-d data
data = squeeze(data);
%data = [data; mean(data)];
if ~exist('newFig','var')
    newFig = 1;
end
if newFig
    figure;hold all;
end
if ~iscell(data)
    if exist('norm','var') && norm ~= 0
        data = bsxfun(@rdivide,data,max(1e-20,std(data,0,2)));
    end
    if ~exist('xaxis','var') || isempty(xaxis)
        xaxis = 1:size(data,2);
    end
    spacing = std(data(:));%5*
    if isreal(data)
        plot(xaxis,data'-repmat(spacing*(1:size(data,1))',[1 size(data,2)])','linewidth',2);
    else
        plot(xaxis,bsxfun(@minus,real(data.'),spacing*(1:size(data,1))),'b',...
            xaxis,bsxfun(@minus,imag(data.'),spacing*(1:size(data,1))),'r','linewidth',2);
    end
else
    offSet = 0;
    for i = 1:numel(data)
        if exist('norm','var') && norm ~= 0
            data{i} = bsxfun(@rdivide,data{i},max(1e-20,std(data{i},0,2)));
            data{i} = bsxfun(@minus,data{i},mean(data{i},2));
        end
        spacing = 5*std(data{i}(:));
        plot(linspace(0,1,size(data{i},2)),bsxfun(@minus,real(data{i}.'),offSet+spacing*(1:size(data{i},1))));
        offSet = offSet + spacing*size(data{i},1);
    end
end
axis tight;
%hold all;
%plot([0 0],[-size(data,1) 0],'k--','linewidth',2);
% set(gca,'xlim',[xaxis(1) xaxis(end)],'ylim',spacing*[-size(data,1)-1 .5],'fontsize',20,'ytick',[]);
% xlabel 'time (s)';

% if isreal(data)
%     plot(data'-repmat(spacing*(1:size(data,1))',[1 size(data,2)])');
% else
%     plot(real(data')-repmat(spacing*(1:size(data,1))',[1 size(data,2)])','b');
%     hold on;
%     plot(imag(data')-repmat(spacing*(1:size(data,1))',[1 size(data,2)])','r');
% end