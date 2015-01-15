%% scripts for making figures for poster
%% figure for correlation coefficients
xs = ((1:140)-40)/80;figure;plot(xs,b,'linewidth',2);hold all;plot(xs,b1,'linewidth',2);plot(xs,b2,'linewidth',2);plot(xs,b3,'linewidth',2);
plot(xs,[d.ccL],'k--');
plot(xs,[d.ccH],'k--');
plot(xs,[d1.ccH],'k--');
plot(xs,[d1.ccL],'k--');
plot(xs,[d2.ccL],'k--');
plot(xs,[d2.ccH],'k--');
plot(xs,[d3.ccH],'k--');
plot(xs,[d2.ccL],'k--');
plot(xs,[d3.ccL],'k--');
xlim([0 1.2]);
ylim([-.05 .7])
xlabel('Time (s)');
ylabel('Correlation coefficient');
title('Prediction of ILD from Eye and Hand trajectory');
set(gca,'fontsize',16);
xlabel('Time (s)');
ylabel('Correlation coefficient');

%% figure for weights title('Time-varying kernel for predicting ILD');
figure;subplot(2,2,1);
imagesc(xs,1:120,squeeze(a)'./[d.kernH],[-1 1]*4);hold all;
set(gca,'ytick',[1 40 60 100],'yticklabel',[0 .5 0 .5],'fontsize',16);
plot(xs,60*ones(size(xs)),'k--','linewidth',2);
xlim([.15 max(xs)]);
title('Subject 1');
ylabel('Post-stimulus kernel (s)');
%xlabel('Time following stimulus onset (s)');

subplot(2,2,2);
imagesc(xs,1:120,squeeze(a1)'./[d1.kernH],[-1 1]*4);hold all;
set(gca,'ytick',[1 40 60 100],'yticklabel',[],'fontsize',16);
plot(xs,60*ones(size(xs)),'k--','linewidth',2);
xlim([.15 max(xs)]);
title('Subject 2');
%xlabel('Time following stimulus onset (s)');

subplot(2,2,3);
imagesc(xs,1:120,squeeze(a2)'./[d2.kernH],[-1 1]*4);hold all;
set(gca,'ytick',[1 40 60 100],'yticklabel',[],'fontsize',16);
plot(xs,60*ones(size(xs)),'k--','linewidth',2);
xlim([.15 max(xs)]);
ylabel('Muscle              Eye');
title('Subject 3');
xlabel('Time following stimulus onset (s)');

subplot(2,2,4);
imagesc(xs,1:120,squeeze(a3)'./[d3.kernH],[-1 1]*4);hold all;
set(gca,'ytick',[1 40 60 100],'yticklabel',[],'fontsize',16);
plot(xs,60*ones(size(xs)),'k--','linewidth',2);
xlim([.15 max(xs)]);
title('Subject 4');
%xlabel('Time following stimulus onset (s)');

%% figure 3

f = f(:,113:end);
f1 = f1(:,113:end);
f2 = f2(:,113:end);
f3 = f3(:,113:end);
e = e(:,113:end);
e1 = e1(:,113:end);
e2 = e2(:,113:end);
e3 = e3(:,113:end);

es = [e;e1;e2;e3];fs = [f;f1;f2;f3];

esa = bsxfun(@minus,es,mean(es,2));
fsa = bsxfun(@minus,fs,mean(fs,2));

params.Fs = 80;params.tapers = [3 5];
params.trialave = 1;params.err = [2 .05];
[Cs,phis,S12,S1,S2,f,x,y,z] = coherencyc(es',fs',params);

inds = [ones(1,27) 2*ones(1,36) 3*ones(1,36) 4*ones(1,36)];
[C,phi,S12,S1,S2,fs,x,y,z] = coherencyc(esa(inds == 1,:)',fsa(inds == 1,:)',params);
[C1,phi1,S12,S1,S2,fs,x1,y1,z1] = coherencyc(esa(inds == 2,:)',fsa(inds == 2,:)',params);
[C2,phi2,S12,S1,S2,fs,x2,y2,z2] = coherencyc(esa(inds == 3,:)',fsa(inds == 3,:)',params);
[C3,phi3,S12,S1,S2,fs,x3,y3,z3] = coherencyc(esa(inds == 4,:)',fsa(inds == 4,:)',params);

figure;plot(fs,C,'linewidth',2);
hold all;plot(fs,C1,'linewidth',2);
hold all;plot(fs,C2,'linewidth',2);
hold all;plot(fs,C3,'linewidth',2);

%% figure 4

allCCs1 = [diag(corr(e',f')); diag(corr(e1',f1')); diag(corr(e2',f2')); diag(corr(e3',f3'))];
