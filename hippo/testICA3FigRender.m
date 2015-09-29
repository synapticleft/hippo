function testICA3FigRender(varExp,widths,cs)
%% this function is fed variables created by testICA3, in order to generate
%% example plots, as shown in a supplementary figure of the paper.

if ndims(widths) == 4
    for i = 1:11
        for j = 1:11
            for k = 1:2
                temp = abs(resample(squeeze(widths(i,j,k,:)),10,1));
                xcInt(i,j,k) = sum(temp/max(temp) > .5);
            end
        end
    end
end


inds = [1 4; 4 4; 10 4; 4 10];
allSNR = linspace(0,.5,11);
allSmooth = max(1,linspace(0,50,11));

figure;
subplot(121);imagesc(allSmooth,allSNR,varExp);
axis square;colorbar;
set(gca,'fontsize',16);
title('Variance explained (r^2)');
xlabel('Place field size (half-width)');
ylabel('Trial-by-trial variation (CV)');
set(gca,'fontsize',16);
subplot(122);
temp = mean(squeeze(xcInt(:,:,1)));
set(gca,'nextPlot','add','ColorOrder',jet(11));
for i= 11:-1:1
indsa = varExp(i,:) > 1/3;
plot(temp(indsa)/10,xcInt(i,indsa,2)/10,'linewidth',2);hold all;
end
%freezeColors;
%subplot(122);freezeColors;colormap(flipud(jet(64)));
%imagesc(allSmooth,allSNR,widths/10);
%set(gca,'fontsize',16);axis square;colorbar;
%title ('FFP size (half-width)');
%xlabel('Place field size (half-width)');
%ylabel('Trial-by-trial variation (CV)');
%set(gca,'fontsize',16);

figure;
col = colormap(hsv(100));
for i = 1:size(inds,1)
    subplot(2,2,i);
    f = find(max(cs{inds(i,1),inds(i,2)}') > 2);
    for j = 1:size(cs{inds(i,1),inds(i,2)},1)
        [~,m] = max(cs{inds(i,1),inds(i,2)}(j,:));
        cc(j,:) = col(m,:);
    end
    set(gca,'nextPlot','add','ColorOrder',cc(f,:));
    plot(cs{inds(i,1),inds(i,2)}(f,:)','linewidth',2);
    set(gca,'fontsize',16,'color',[0 0 0],'ylim',[-1 max(cs{inds(i,1),inds(i,2)}(:))]);
    if ismember(i,[1 2])
        set(gca,'xtick',[]);
    end
    title([num2str(numel(f)) ' FFPs']);
    % (' num2str(round(numel(f)/size(cs{inds(i,1),inds(i,2)},1)*100)) '%)']);
    %, r^2 = ' num2str(round(100*varExp(inds(i,1),inds(i,2)))/100)]);
end
subplot(2,2,3);xlabel 'Position';ylabel 'Activation (z-score)';

%%%%%%%%%%%