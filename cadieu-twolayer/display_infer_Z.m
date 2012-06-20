function display_infer_Z(Z,I,Ih,m,p)
% Z = double(a.*exp(1j*phase));
a = abs(Z);
% Z(2:end,:) = bsxfun(@times,Z(2:end,:),exp(1i*-angle(Z(1,:))));
% Z(1,2:end) = Z(1,2:end).*exp(1i*-angle(Z(1,1:end-1)));
phase = angle(Z);
if p.whiten_patches
    I  = bsxfun(@plus,m.dewhitenMatrix*I, m.imageMean);
    Ih = bsxfun(@plus,m.dewhitenMatrix*Ih,m.imageMean);
end

%phase = phase + -2*pi*sign(phase).*round(abs(phase)./(2*pi));

sfigure(5);
subplot(2,2,2);
hval=max(max(abs(a(:))),.1);
imagesc(a,[0 1]*hval), axis off, colormap gray
title('a')
subplot(2,2,4);
imagesc(phase,[-pi pi]), axis off, colormap hsv
alpha(double(a/max(a(:))));
freezeColors
title('phase')

subplot(2,2,1);
Ival=max(abs(I(:)));
imagesc(real(I),[-1 1]*Ival), axis off, colormap gray
title('I')
subplot(2,2,3);
imagesc(real(Ih),[-1 1]*Ival), axis off, colormap gray
title('Ihat')
drawnow;