function soundMaria

global slowDown
global pNew
global arDecay
global arRate
global gain
global pixFilt

close all;
pNew = .01;arDecay = .01;arRate = pi/4;gain = 1;pixFilt = 10;
slowDown = 20;
VIDSIZE = [480 854];
%path = 'C:\wiki\';
%allFiles = [path 'JSIms1.avi'];
allFiles = 'Baaa1.avi';
state = ones(VIDSIZE(1),VIDSIZE(2),3);
lineInfo = zeros(1,VIDSIZE(2));
tic;
%%%%%%%%%%%%figure initialization, could be put in another function
f = figure('toolbar','none','menu','none','Name','','KeyPressFcn',@changeProps);
%maximize(f);
subplot('Position',[0 0 1 1]);   %[0 .2 1 .8]
h = imagesc(zeros(VIDSIZE(1),VIDSIZE(2),3));
set(gca,'xtick',[], 'ytick',[],'color','k');
vr = VideoReader(allFiles);%,'preciseFrames',0);
%info = getinfo(vr);
i = 1;
j = slowDown - 1;
while 1
    j = j + 1;
    if j >= slowDown
        %seek(vr,i);
        curFrame = readFrame(vr);%curFrame = getframe(vr);
        j = 0;
%         if i < 1000 %info.numFrames
%             i = i + 1;
%         else
%             i = 1;
%         end
    end
    ra = getRatio(j/slowDown);
    temp = rand(size(lineInfo)) < pNew;
    ar = (1-arDecay)*exp(1i*arRate);
    lineInfo = idealfilterG(temp,pixFilt);
    state = state*ar + bsxfun(@times,double(curFrame)/256,lineInfo/max(.01,max(lineInfo)));
    %state = state*(1-arDecay);% + bsxfun(@times,double(curFrame)/256,max(0,));
    %state(:,lineInfo > 0,:) = bsxfun(@times,double(curFrame(:,lineInfo>0,:)),lineInfo(lineInfo > 0));
    set(h,'CData',max(0,min(1,real(state)/gain)));
    drawnow;
end

function ra = getRatio(frac)
ra = max(0,frac*2-1);

function changeProps(obj,event)
global pNew
global arDecay
global arRate
global gain
global pixFilt
global slowDown

if event.Key == 'z'
    pNew = pNew/2;
elseif event.Key == 'q'
    pNew = pNew*2;
elseif event.Key == 'a'
    pNew = .01;
elseif event.Key == 'x'
    arDecay = arDecay*2;
elseif event.Key == 'w'
    arDecay = arDecay/2;
elseif event.Key == 's'
    arDecay = .01;
elseif event.Key == 'c'
    arRate = arRate/1.2;
elseif event.Key == 'e'
    arRate = min(arRate*1.2,pi);
elseif event.Key == 'd'
    arRate = pi/4;
elseif event.Key == 'v'
    gain = gain*1.1;
elseif event.Key == 'r'
    gain = gain/1.1;
elseif event.Key == 'f'
    gain = 1;
elseif event.Key == 'b'
    pixFilt = pixFilt*2;
elseif event.Key == 't'
    pixFilt = pixFilt/2;
elseif event.Key == 'g'
    pixFilt = 10;
elseif event.Key == 'y'
    slowDown = max(1,slowDown - 1);
elseif event.Key == 'h'
    slowDown = 20;
elseif event.Key == 'n'
    slowDown = slowDown + 1;
end
[pNew arDecay arRate gain pixFilt slowDown]
