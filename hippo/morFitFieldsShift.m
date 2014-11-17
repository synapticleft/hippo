function cc = morFitFieldsShift(X,y,pos)
%% from filtered fits of spike trains, make trials x position activity maps

accumbins = 100;
dec = 32;
decs = 4;
levels = 7;

[~,~,isFast] = fixPos(pos);
isFast1 = logical(round(interp([double(isFast); 0],dec/decs)));
%isFast1 = logical(ones(1,size(X,2)));
y = morFilter(y,2^levels,1250/decs);
X = morFilter(X,2^levels,1250/decs);
X = X(:,isFast1);
Xnorm = (X*X')\X;

shifts = -30:30;
for i = 1:numel(shifts)
    temp = circshift(y,[0 shifts(i)]);%y(:,isFast1);
    W = temp(:,isFast1)*Xnorm';
    cc(i,:) = diag(temp(:,isFast1)*(W*X)');
end