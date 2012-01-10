function fftReduce(data,fs,window,rdim)

params.tapers = [1/window window 1];
params.Fs = fs;

[J,t,f] = gmtspecgramc(data',[window window],params);

%J = permute(J,[3 2 1]);
Jm = mean(J,3);
%[u,s,v] = svds(J(:,:),10);
%figure;plot(diag(s));
J = bsxfun(@times,J,exp(-1i*angle(Jm)));
%figure;plot(squeeze(mean(abs(Jm))));
% figure;
% lims = [-1 1]/2;
% for i = 1:size(J,1)
%     scatter(squeeze(real(J(i,5,:))),squeeze(imag(J(i,5,:))));
%     set(gca,'xlim',lims,'ylim',lims);
%     drawnow;pause(.1);
% end
[u,s,v] = svds(J(:,:),rdim);
figure;hold all;
plot(mean(abs(squeeze(Jm)))/4);
for i = 1:size(v,2)
    plot(mean(abs(reshape(v(:,i),[numel(f) size(data,1)]))'));
    %subplot(size(v,2),1,i);
    %imagesc(f,1:size(data,1),abs(reshape(v(:,i),[numel(f) size(data,1)]))');
end

figure;
for i = 1:size(v,2)
    subplot(size(v,2),1,i);
    imagesc(f,1:size(data,1),angle(reshape(v(:,i),[numel(f) size(data,1)]))');
    plot(f,angle(reshape(v(:,i),[numel(f) size(data,1)]))');
end
%[a b c] = ica_lfp_complex(J(:,:).',20);
% [size(a) size(b)]
%hold all;plot(diag(s));
%figure;plot(abs(u));
%figure;imagesc(t,f,squeeze(abs(Jm)));
%figure;imagesc(t,f,squeeze(angle(J(:,9,:))));

function [f g] = minGrad(