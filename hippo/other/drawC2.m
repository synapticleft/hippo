function drawC2(u,ve)

[~,ord] = sort(circ_dist(angle(u(:,1)),circ_mean(angle(u(:,1)))),'descend');temp = (exp(1i*(-pi:.1:pi)));
u1 = u(ord,1)/max(abs(u(:,1)));
u2 = u(ord,2)/max(abs(u(:,2)));
t1 = exp(1i*linspace(0,4*pi,50));
t2 = exp(1i*linspace(0,5*pi,50)).*linspace(0,1,50);
c = getCol(size(u,1));
subplot(3,5,3);
plot(real(temp),imag(temp),'k--','linewidth',2);hold all;
scatter(real(u1),imag(u1),30,c,'o','filled');
axis image ;
set(gca,'fontsize',16,'xtick',[],'ytick',[]);
title 'Complex Value'
subplot(3,5,8);
plot(real(temp),imag(temp),'r--','linewidth',2);hold all;
scatter(real(u2),imag(u2),30,c,'^','filled');axis image 
set(gca,'fontsize',16,'xtick',[],'ytick',[],'ylim',[-1 1],'xlim',[-1 1]);
subplot(3,5,13);
plot(real(temp),imag(temp),'k--','linewidth',2);hold all;
plot(real(t2),imag(t2),'r--','linewidth',1.5)
scatter(real(u2),imag(u2),30,c,'^','filled');
scatter(real(u1),imag(u1),30,c,'o','filled');
set(gca,'fontsize',16,'xtick',[],'ytick',[]);axis image ;
xlabel 'Real';ylabel 'Imag';
subplot(3,5,4:5);
one64 = u1*t1;
plot(1,1);hold on;
plot(real(t1),'k--');
set(gca,'ColorOrder',c);plot(real(one64)','linewidth',2);axis tight ;
set(gca,'fontsize',16,'xtick',[],'ytick',[]);title 'Simulated Electrode Time Courses'
subplot(3,5,9:10);
two64 = u2*t1;
plot(1,1);hold on;
plot(real(t1),'r--');
set(gca,'colororder',c);plot(real(two64)','linewidth',2);axis tight ;
set(gca,'fontsize',16,'xtick',[],'ytick',[]);
ylabel 'Real component of Electrode Activations'
subplot(3,5,14:15);
three64 = u1*t1+u2*t2;
plot(1,1);hold on;
plot(real(t1),'k--');plot(real(t2),'r--');
set(gca,'colororder',c);plot(real(three64)','linewidth',2);axis tight ;
set(gca,'fontsize',16,'xtick',[],'ytick',[]);
xlabel 'Time'; 
subplot(3,5,1);imagesc(flipud(reshape(real(u(:,1))/max(abs(u(:,1))),[8 8])),[-1 1]);
colorbar('southoutside');set(gca,'fontsize',16,'xtick',[],'ytick',[]);ylabel 'PC 1';title 'Real Component';
subplot(3,5,6);imagesc(flipud(reshape(real(u(:,2))/max(abs(u(:,2))),[8 8])),[-1 1]);
colorbar('southoutside');set(gca,'fontsize',16,'xtick',[],'ytick',[]);ylabel 'PC 2'
subplot(3,5,11);set(gca,'fontsize',16,'xtick',[],'ytick',[]);ylabel 'PC 1 & 2';
for i = 1:3
    subplot(3,5,2+5*(i-1));
    if i == 2
        imagesc(flipud(reshape(abs(ve(i,:)).^2*100,[8 8])),[0 8]);
    else
        imagesc(flipud(reshape(abs(ve(i,:)).^2*100,[8 8])),[85 100]);
    end
    axis off;colorbar('southoutside');set(gca,'fontsize',16);
    if i == 1
        title '% Variance Explained'
    end
end
drawnow;

function c = getCol(n)
c = repmat(linspace(0,1,n)',[1 3]);
for i = 1:3
c(:,i) = max(0,min((1-abs(c(:,i)-i/4)*2.5)*1.5,1));
end
c = fliplr(c);