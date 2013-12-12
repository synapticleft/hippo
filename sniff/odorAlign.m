function [ allSamps ] = odorAlign(pcaLFP, tr, trial)
allSamps = zeros(numel(tr),size(pcaLFP,1),300);
k = 0
for i = [2 4 5 6 7 8 9 10]
    for j = 1:numel(tr)
        if trial(j).odorValve == i
            k = k + 1;
            allSamps(k,:,:) = pcaLFP(:,tr(j)+(-99:200));
        end
    end
end
