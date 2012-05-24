function makeDist(ampRange)

yInt = 1;
slope = 10;%.45;
theta = pi/4;
sig = 1;

spacing = linspace(-ampRange,ampRange,100)';
[x y] = meshgrid(spacing);
c = complex(x,y);
k = yInt + slope*abs(c);
vm = getVM(c,slope,theta);
g = getGauss(c,sig);
% figure;hold all;
% xs = -pi:.1:pi;
% for i = 1:ampRange 
%     plot(xs,getVM(i*exp(xs*1i),slope,theta));
% end
figure;imagesc(spacing,spacing,log(g));
figure;imagesc(spacing,spacing,log(vm));
figure;imagesc(spacing,spacing,log(vm.*g));

%imagesc(exp(abs(c).*cos(angle(c)-pi/3))./besseli(0,abs(c)));
function vm = getVM(dat,slope,theta)
k = slope*abs(dat);
vm = exp(k.*cos(angle(dat)-theta))./(2*pi.*besseli(0,abs(dat)));%abs(dat)

function g = getGauss(dat,sigma)
dat = abs(dat);
g = exp(-(dat.^2)/(2*sigma^2))/sqrt(2*pi*sigma^2);