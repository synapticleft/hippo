function [ST, y, yHat, kern, xY] = pressTriggered(keyDat,board,sr)
%% predict keypresses from physiological data. Newer than predictKeys (preserves spike-triggered data).

%keyDat(keyDat(:,2) == 0,:) = [];
keysPressed = unique(keyDat(:,2));
keysPressed(keysPressed <= .5) = [];
up = zeros(numel(keysPressed),size(board,2));
down = up;
for i = 1:numel(keysPressed)
    down(i,round(keyDat(keyDat(:,2) == keysPressed(i),1)*sr)) = 1;
    up(i,round(keyDat(keyDat(:,2) == -keysPressed(i),1)*sr)) = 1;
end
%down(end+1,:) = sum(down);
%up(end+1,:) = sum(up);
%press = cumsum(down,2)-cumsum(up,2);

y = [down;up];%];%;press
%board = filtHigh(board,sr,.1);
X = makeToeplitz(board,sr*.2,1);
for i = 1:size(y,1)
    ST{i} = X(:,y(i,:) > 0);
end
return
xY = y*X';
%y = zscore(y,0,2);
ridges = [0 10.^(-4:4)];
[cc, mse , kern] = ridgeCross(y',X',3,ridges);
figure;plot(mse);
figure;plot(cc);
yHat = squeeze(kern(round(end/2)+1,:,:))'*X;
figure;plot(y(6,:));hold all;plot(yHat(6,:));

return
scales = 2.^(1:7);
wave = zeros(size(board,1)*numel(scales),size(board,2));
for i = 1:size(board,1)
    wave((i-1)*numel(scales)+(1:numel(scales)),:) = cwt(board(i,:),scales,'cmor1-1');
end
temp = 0;


for i = 1:size(y,1)
    for j = 1:size(board,1)
        temp = temp + xcov(y(i,:),abs(wave(numel(scales)*(j-1)+3,:)),sr);
    end
end
figure;plot(temp);
X = zscore(abs(wave),0,2);
w = y/X;
yHat = w*X;
figure;plot(y'/100,'b');hold all;plot(yHat','r');
corr(yHat',y')

