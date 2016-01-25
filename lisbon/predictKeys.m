function [y yHat kern xY] = predictKeys(keyDat,board,sr)
%% predict keypresses from physiological data

%keyDat(keyDat(:,2) == 0,:) = [];
keysPressed = unique(keyDat(:,2));
keysPressed(keysPressed <= 0) = [];
up = zeros(numel(keysPressed),size(board,2));
down = up;
for i = 1:numel(keysPressed)
    down(i,round(keyDat(keyDat(:,2) == keysPressed(i),1)*sr)) = 1;
    up(i,round(keyDat(keyDat(:,2) == -keysPressed(i),1)*sr)) = 1;
end
%down(end+1,:) = sum(down);
%up(end+1,:) = sum(up);
%press = cumsum(down,2)-cumsum(up,2);

y = [down;up];%;press
board = filtHigh(board,sr,.1);
X = makeToeplitz(board,sr*1,1);
xY = y*X';
y = zscore(y,0,2);
ridges = [0 10.^(-4:4)];
[cc, mse , kern] = ridgeCross(y',X',3,ridges);
figure;plot(mse);
figure;plot(cc);
yHat = squeeze(kern(round(end/2)+1,:,:))'*X;
figure;plot(y(2,:));hold all;plot(yHat(2,:));

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

function X = makeToeplitz(data,lags,useZ)
if useZ
    data = zscore(data,0,2);
end
X = [];
for i = 1:size(data,1)
    X = [X; toeplitz(data(i,:),zeros(lags,1))'];
end
X = circshift(X,[0 -round(lags/2)]);
X(end+1,:) = 1;