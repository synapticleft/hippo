function [a b inds] = propIm(trace,f1,f2,t,slices)

%slices = 10;
Fs = 1250;
range = (f1+f2)/2;range = 100;%max(floor(Fs/range/10)
traceOr = trace;
trace = hipFilter(trace,f1,f2);
trace = circshift(trace,[-1 0]);
trace(64,:) = trace(62,:);
H1 = angle(hilbert(trace'))';
ind = t*Fs;

while H1(28,ind) <= H1(28,ind+1)
    ind = ind+1;
end
ind = ind + 1;inds(1) = ind;
a(1,:) = trace(:,ind);
b(1,:) = unwrapped(H1,ind,range);%H1(:,ind);
for i = 1:(slices-1)
    while H1(28,ind) < -pi + i*pi*2/(slices-1) && H1(28,ind) < (H1(28,ind + 1)+pi)
        ind = ind + 1;
    end
    inds(i+1) = ind;
    a(i+1,:) = trace(:,ind);
    b(i+1,:) = unwrapped(H1,ind,range);
end

% range = 500;
%figure;plot((inds(1)-range):(inds(5)+range),trace(5:8:end,(inds(1)-range):(inds(5)+range)));
%figure;plot((inds(1)-range):(inds(5)+range),traceOr(10:10:end,(inds(1)-range):(inds(5)+range)));
%figure;plot((inds(1)-range):(inds(5)+range),H1(5:8:end,(inds(1)-range):(inds(5)+range)));
%hold on;scatter(inds,H1(28,inds),'r');

function b = unwrapped(a,ind,r)

%temp = unwrap(a(:,(ind-r):(ind+r))')';
%temp = temp(:,r+1);
temp = a(:,ind);
temp = reshape(temp,[8 8]);%reshape(a(:,ind),[8 8]);
temp = unwrap(temp')';
temp = unwrap(temp);
b = makeflat(temp);
b = a(:,ind);
b = b - b(28,:);

