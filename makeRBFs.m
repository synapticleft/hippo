function a = makeRBFs(lengthBase,numBases,numTrials,useDiff)
%% INPUTS
% lengthbase - length of a single hump (i.e. if sampling rate of sound is
% 10 kHz, and you want .1-s humps, use lengthbase = 10000*.1 = 1000
% numBases - number of humps in stimulus. For the above example, if you
% want a 1 s stimulus, use 1/.1 = 10 numBases
% numTrials - number of trials to generate stimuli for
% useDiff - whether to generate stimuli that have 0-mean fluctuation. NOTE:
% the requirement for stimulus to be zero-sum, and to smoothly transition
% between 'lobes' allows the middle lobes to reach greater height than edge
% ones.
% OUTPUT
% a - numTrials x [lengthbase*numBases] sized matrix of scaling factors.
% this should be added to the appropriate DC component and scaled as needed

a = zeros(numTrials,lengthBase*numBases);

wave= 1-cos((1:lengthBase)/lengthBase*2*pi);
if useDiff
    wave = [wave -wave];
end

for i = 1:numTrials
    if useDiff
        for j = 1:numBases*2-3
            a(i,(j-1)*lengthBase/2+(1:2*lengthBase)) = a(i,(j-1)*lengthBase/2+(1:2*lengthBase))+(ceil(rand*3)-2)*wave;
        end
    else
    for j = 1:numBases*2-1
        a(i,(j-1)*lengthBase/2+(1:lengthBase)) = a(i,(j-1)*lengthBase/2+(1:lengthBase))+(ceil(rand*3)-2)*wave;
    end
    end
end