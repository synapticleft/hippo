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