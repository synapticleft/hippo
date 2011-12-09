function [shifts, meanTrace] = alignTraces(data,window,step)

maxOffset = 30;
numSteps = floor(size(data,2)/step);
shifts = zeros(numSteps,size(data,1));
meanTrace = zeros(numSteps,window);
figure;
for i = 1:numSteps
    chunk = data(:,(i-1)*step+(1:window));
    subplot(411);plot(chunk');
    subplot(414);plot(mean(chunk));hold all;
    [chunk,shifts(i,:)] = findAlign1(chunk,maxOffset,2,4);
    %[chunk,shifts(i,:,2)] = findAlign1(chunk,maxOffset,3,4);
    meanTrace(i,:) = mean(chunk);
    hold off;drawnow;
    %return
end
shifts(shifts > maxOffset) = window - shifts(shifts > maxOffset);


function [chunk sh] = findAlign1(chunk,offset,sub1,sub2)
cm = mean(chunk);
sh = zeros(1,size(chunk,1));
peaks = zeros(size(chunk,1),offset*2+1);
for i = 1:size(chunk,1)
    peaks(i,:) = xcorr(chunk(i,:),cm,offset,'unbiased');
    [~,sh(i)] = max(peaks(i,:));
    sh(i) = -(sh(i)-offset-1);
    chunk(i,:) = circshift(chunk(i,:),[0 sh(i)]);
end
if nargin > 2
    subplot(sub2,1,sub1);plot(chunk');%peaks');%
    subplot(sub2,1,sub2);plot(mean(chunk));
end

function [chunk sh] = findAlign(chunk,offset,sub1,sub2)
fData = fft(chunk,[],2);
fm = mean(fData);
peak = bsxfun(@times,fData,conj(fm));
peak = fliplr(ifft(peak,[],2));%./abs(peak)
peak(:,(offset+1):(end-offset)) = 0;
[~,sh] = max(peak,[],2);
chunk = shiftAll(chunk, sh);
if nargin > 2
    subplot(sub2,1,sub1);plot(chunk');%peak');axis tight;%;plot(
    subplot(sub2,1,sub2);plot(mean(chunk));
end

function chunk = shiftAll(chunk,theShifts)
for i = 1:size(chunk,1)
    chunk(i,:) = circshift(chunk(i,:),[0 theShifts(i)]);
end