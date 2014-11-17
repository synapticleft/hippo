function floatingOrb(Xfd,u,probes1,pos,spacer,fileName)

fs = 4;
gain = .5;
temp = 1:numel(probes1);temp = reshape(temp,size(probes1));
temp = complexIm(u(temp),0,gain);
temp = permute(temp,[3 1 2]);
uc = temp(:,:)';
Xfd = Xfd/std(Xfd(:));
bounds = 2.5;


for i= 1:spacer:450%size(Xfd,2)
    subplot(2,3,[1 2 4 5]);
    scatter(real(Xfd(:,i)),imag(Xfd(:,i)),[],uc,'filled');
    axis image;
    set(gca,'xlim',bounds*[-1 1],'ylim',bounds*[-1 1]);
    title([num2str(i*fs/1250,2) ' s']);
    xlabel 'Real'
    ylabel 'Imaginary'
    set(gca,'Xtick',[],'Ytick',[],'fontsize',16);
    temp = Xfd(:,i);
    subplot(2,3,3);plot(pos(max(1,i-100):i,1),pos(max(1,i-100):i,2),'k','linewidth',2);
    axis square;
    xlabel 'X'
    ylabel 'Y'
    set(gca,'Xtick',[],'Ytick',[],'xlim',[min(pos(:,1))-5 max(pos(:,1))+5],...
        'ylim',[min(pos(:,2))-5 max(pos(:,2))+5],'fontsize',16);
    subplot(2,3,6);imagesc(complexIm(temp(probes1+1),0,gain,1,bounds));
    axis square off;
    drawnow;
    if exist('fileName','var')
        m(ceil(i/spacer)) = getframe(gcf);
    end
end

if exist('fileName','var')
    movie2avi(m,fileName);
end