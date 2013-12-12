
%% Plot Decoding Results...

idx=1:size(f,1);
figure(1)
clf
subplot(2,1,1)
h = fspecial('disk',2);
% imagesc(1:size(f,1),pvec,f')
imagesc(1:size(f,1),pvec,imfilter(f,h)')
% imagesc(f')
%xlim(xl)
hold on
plot(1:length(idx),rpos(idx)*pi,'r.')
[tmp,maxpost]=max(f');
% [tmp,maxpost]=max(imfilter(f,h)');
% baypost = sum(bsxfun(@times,f,1:40),2);
plot(1:length(idx),maxpost/length(pvec)*2*pi,'k.')
% plot(1:length(rpos(1:25:end)),baypost,'k^')
% plot(1:(length(idx)-1),diff(rpos(idx)).^2*250,'y')
hold off
xlabel('Time [bins]')
ylabel('Position')

subplot(2,1,2)
bar(mean(spk(:,idx)),'k')
ylabel('Avg Spikes')
%xlim(xl)
box off

%% Confusion matrix...

edges = [pvec-mean(diff(pvec))/2 max(pvec)+mean(diff(pvec))/2];
[tmp,rpos_bin] = histc(rpos(idx)*pi,edges);

C=[];
for i=1:length(pvec)
    bidx = find(rpos_bin==i);
    C(:,i) = mean(f(bidx,:),1);
end
if all(C(:)>=0)
    C = bsxfun(@rdivide,C,sum(C));
end
figure(2)
clf
imagesc((C))
axis image
colorbar
title('Confusion Matrix')
xlabel('True Position')
ylabel('Predicted Position')

%% Circular error...

err = circ_dist(rpos(idx)*pi,maxpost'/length(pvec)*2*pi);

figure(3)
hist(abs(err),100)
axis tight
xlabel('Absolute Error')
ylabel('Frequency')