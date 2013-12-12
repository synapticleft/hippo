function [ trialLFP ] = onlyTrials(lfp, tr)
trialLFP = zeros(size(lfp,1),numel(tr)*300);
for i = 1:31
    k = 0;
    for trial = 1:numel(tr)
        trialLFP(i,((300*k)+1):(300*(k+1))) = lfp(i,tr(trial)+(-99:200));
        k = k + 1;
    end
end