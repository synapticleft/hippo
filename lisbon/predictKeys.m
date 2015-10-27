function predictKeys(keyDat,board,sr)
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
press = cumsum(down,2)-cumsum(up,2);

%clear up down press;

scales = 2.^(1:7);
wave = zeros(size(board,1)*numel(scales),size(board,2));
for i = 1:size(board,1)
    wave((i-1)*numel(scales)+(1:numel(scales)),:) = cwt(board(i,:),scales,'cmor1-1');
end
temp = 0;
down = circshift(down,[0 24]);
y = down;%y = [down;up;press];
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