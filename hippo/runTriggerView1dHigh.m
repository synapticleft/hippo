function t1 = runTriggerView1dHigh(pos,t,accumbins,scale)
%% 1d binning of ICA activations on linear track

[~,posd,~,runs] = fixPos(pos);
posd = max(1,min(2*accumbins,round(resample(posd,scale,1)*accumbins)));
runs = double(runs);
runs = max(1,min(max(runs),round(resample(runs,scale,1))));
runs = runs(1:size(t,2));posd = posd(1:size(t,2));
for j = 1:size(t,1)
         t1(j,:,:) = accumarray([runs' posd],t(j,:)'.^2,[max(runs) 2*accumbins(1)] ,@mean);
end