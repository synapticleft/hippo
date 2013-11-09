[a,b,c,d] = LoadCluRes('ec014.468',4);
s = LoadSpk('ec014.468.spk.4',8);
%% plot sorted spike waveforms
figure;for i = 1:max(b)
subplot(4,4,i);imagesc(squeeze(mean(s(:,:,b == i),3)));
end
%% plot normalized, sorted spike waveforms
figure;for i = 1:max(b)
temp = squeeze(mean(s(:,:,b == i),3));
subplot(4,4,i);imagesc(bsxfun(@rdivide,temp,std(temp,0,2)));
end
%% accumulate spikes by sorted cell type
a1 = min(size(v,1),max(1,round(a/32/16)));
temp = [];
for i = 1:max(b)
temp(i,:) = accumarray(a1(b ==i),1,[size(v,1) 1],@sum);
end
temp1 = morFilter(temp,8,1250/32);
t = runTriggerView1d(pos,v,temp1,100,.05);
%% load data
X = getData('ec014.468',1,25:32);
X = filtHigh...
[x,y,z] = fastICA(X,'approach','symm');
for i = 1:8
x1(i,:) = decimate(x(i,:).^2,8);
end
x2 = morFilter(x1,8,1250/32);
%% ICA on spike power
[x,y,z] = fastICA(squeeze(std(s,0,2)),'approach','symm');
temp = [];
for i = 1:8
temp(i,:) = accumarray(a1,x(i,:),[size(v,1) 1],@mean);
end
%% fitting to GLM
ss = squeeze(std(double(s),0,2));
spM = zeros(max(b),size(a,1));
for i =3:max(b)
spM(i,b == i) = 1;
end
for i = 3:size(spM,1)
lin(i,:) = glmfit(double(ss'),spM(i,:)','normal');
end
for i = 3:size(spM,1)
yHat(i,:) = glmval(lin(i,:)',ss','identity');
end
for i = 3:size(spM,1)
temp(i,:) = accumarray(a1,yHat(i,:),[size(v,1) 1],@mean);
end
temp1 = morFilter(temp,8,1250/32);
t = runTriggerView1d(pos,v,temp1,100,.05);
%% spatiotemporal fit
 s = permute(s,[3 1 2]);
lin = spM(:,b > 2)/double(s(b > 2,:))';
for i = 3:size(spM,1)
yHat(i,:) = lin(i,:)*double(s(:,:)');
end

%% find most corr spike - hfICA
cc = corr(corr(teSp(:,:)',tes1(:,:)'));
inds = [32 25 34 35];f = [8 21 45 36];
sz = 3;
cc = abs(cc);
figure;
for i = 1:4
    subplot(4,sz,(i-1)*3+1);imagesc(complexIm(teSp(inds(i),:,:)));axis off;
    [~,f(i)] = max(cc(inds(i),:));
    subplot(4,sz,(i-1)*3+3);imagesc(reshape(abs(A(:,f(i))),[8 8]));axis off;
    subplot(4,sz,(i-1)*3+2);imagesc(complexIm(tes1(f(i),:,:)));axis off;
end

spacing = .5;
xdim = 1;
ydim = 1.5;
figure;hold all;
for i = 1:4
    imagesc((i+(-1:0))*ydim+spacing*(i-1),[2 3]*xdim+2*spacing,complexIm(teSp(inds(i),:,:)));
    [~,f(i)] = max(cc(inds(i),:));
    imagesc((i+(-1:0))*ydim+spacing*(i-1),[1 2]*xdim+1*spacing,complexIm(tes1(f(i),:,:)));
    imagesc((i+(-1:0))*ydim+spacing*(i-1),[0 1]*xdim,(reshape(abs(A(:,f(i)))/max(abs(A(:,f(i)))),[8 8])));
end
colormap gray;axis off;
for i = 1:4
    for j = 1:4
        ts(i,j) = cc(inds(j),f(i));
        if ts(i,j) > .2
        plot([(i-.5)*ydim+(i-1)*spacing (j-.5)*ydim+(j-1)*spacing],[xdim xdim+spacing]+xdim+spacing,'k','linewidth',ts(i,j)*10);
        end
    end
end
%%
 figure;for i = 1:4
subplot(4,4,i);showGrid(complex(-squeeze(sp1250H(inds(i)-11,:,17)),eps*ones(1,64)),[8 8],[],[],[],[],[],1);
subplot(4,4,i+4);showGrid(complex(sp20000(inds(i)-11,:,:),eps*ones(1,8,32)),[],[],[],[],[],[],1);
subplot(4,4,i+8);showGrid(complex(sp1250(inds(i)-11,9:16,:),eps*ones(1,8,32)),[],[],[],[],[],[],1);
subplot(4,4,i+12);showGrid(complex(sp80(inds(i)-11,9:16,:),eps*ones(1,8,32)),[],[],[],[],[],[],1);
end
%%
for i = 16:-1:1
sp80Samples((i-1)*8+(1:8),:) = X1(:,max(1,min(size(X1,2),spTimes80-i+8)));
end
%%
inds = [5 21 45 3 35 17 13 59 49 39 43 6];
offSet = [1 7 13 19 3 9 15 21 5 11 17 23];
figure;for i = 1:12
    subplot(4,6,offSet(i));showGrid(tes1(inds(i),:,:));
    subplot(4,6,offSet(i)+1);showGrid(A1(:,inds(i))',[8 8],[],[],[],[],[],1);
end