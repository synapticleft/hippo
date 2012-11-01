function ims = colorComplex(u,sz)

ims = zeros(size(u,2),sz(1),sz(2),3);

for i = 1:size(u,2)
    temp = reshape(u(:,i)*exp(-1i*circ_mean(angle(u(:,i)))),sz);
    ims(i,:,:,1) = real(temp);
    ims(i,:,:,2) = imag(temp);
    ims(i,:,:,1:2) = ims(i,:,:,1:2)/max(abs(u(:,i)))/2 + .5;
end
figure;
for i = 1:size(u,2)
    subplot(1,size(u,2),i);imagesc(squeeze(ims(i,:,:,:)));
end