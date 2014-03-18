function ripChans=rippledet1(X,allchan,detchan)%,minuschan)
% when you have the variable already, not the file
% file: file base name without extension (looks for the eeg file)
% allchan: number of total channels (either 256 or 512)
% detchan: detection channel, counted from 1.
% minuschan: one channel where you don't expect ripples. This will be substracted
% from detch, in order to eliminate systematic noise in the ripple frequency band.
% outputext: extension for the output file, i.e.: '.rpl'.

%rootDir = '/media/storage/hippocampus/';%'/media/Expansion Drive/redwood/';
%file = [rootDir file '/' file '.eeg'];
%parameters
sz = 100000;
samplerate= 1250; %S/s
cutofflowpass=250; %Hz
cutoffhighpass=100; %Hz
order = 8;
window= 20; % ms
ripWin = ceil(20*samplerate/1000);
detectionlevel=7; %SD
%boundarylevel=1; %SD

%load channel
%a=memmapfile(file,'Format','int16');
f1 = fdesign.bandpass('N,F3dB1,F3dB2',order,cutoffhighpass,cutofflowpass,samplerate);
%f1 = fdesign.highpass('n,f3dB',order,cutoffhighpass,samplerate);
d1 = design(f1,'butter');
[B,A]= sos2tf(d1.sosMatrix,d1.ScaleValues);
numSamps = floor(numel(X)/allchan);
ripChans = [];
for i = 1:ceil(numSamps/sz)
    inds = (i-1)*sz+[1 sz];
    inds(2) = min(inds(2),numSamps);
    if exist('detchan','var');
        Filtchan = X(((inds(1):inds(2))-1)*allchan+detchan);
    else
        Filtchan = mean(reshape(X((inds(1)-1)*allchan+1:inds(2)*allchan),[allchan inds(2)-inds(1)+1]));
    end
    Filtchan = filtfilt(B,A,double(Filtchan)')';
    %if exist('minuschan','var') && ~isempty(minuschan)
    %    Filtchan = Filtchan - filtfilt(B,A,double(a.data(((inds(1):inds(2))-1)*allchan+detchan))')';
    %end
    %calculate power in x ms windows
    power = zscore(filtfilt(ones(ceil(window*samplerate/1000/2),1),sum(ones(ceil(window*samplerate/1000/2),1)),Filtchan.^2));
    cross= power > detectionlevel;%power>meanpower+detectionlevel*stdpower;
    cross = bwlabel(cross);
    [max(cross) i ceil(numSamps/sz)]
    for j = 1:max(cross)
        [~,ind] = max(power.*(cross == j));
        ind1 = (inds(1)-1+ind-ripWin)*allchan+1;
        ind2 = (inds(1)+ind+ripWin)*allchan;
        dataChunk = reshape(a.data(ind1:ind2),[allchan 2*ripWin+1]);
        dataChunk = filtfilt(B,A,double(dataChunk'))';
        powInd = ind+(-ripWin:ripWin);
        powInd(powInd > numel(power)) = [];powInd(powInd <1) = [];
        plot(dataChunk');hold all;plot(power(powInd)*10,'r','linewidth',2);hold off;drawnow;
        ripChans = [ripChans mean(dataChunk.^2,2)];
    end
end
end
%     crossbd(:,1)=find(cross-circshift(cross,1)==1);
%     crossbd(:,2)=find(cross-circshift(cross,1)==-1);
%     clear cross
%     ripples=zeros(size(crossbd,1),3);
%     for i=1:size(crossbd,1)
%         ripples(i,2)=find(power(crossbd(i,1):crossbd(i,2))==max(power(crossbd(i,1):crossbd(i,2))))+crossbd(i,1);
%         %detect beginning boundary
%         ripples(i,1)=find(power(1:ripples(i,2))<(meanpower+boundarylevel*stdpower),1,'last');
%         ripples(i,3)=find(power(ripples(i,2):end)<(meanpower+boundarylevel*stdpower),1,'first')+ripples(i,2); 
%         %detect end boundary
%     end
%     ripples=ripples+ceil((window*samplerate/1000)/2); %shift values with half window width
%     
%     %delete duplicate enries
%     for i=size(ripples,1):-1:2
%         if (ripples(i,1)==ripples(i-1,1))
%             ripples(i,:)=[];
%         end
%     end
% 
%     %convert timepoints to samples at 20 kHz
%     ripples=floor((ripples.*1000)./samplerate);
%     hely=find(ripples(:,3)>(30*60*1000),1,'first');
%     %save ripples to evt file
%     fid=fopen([file '.'  outputext '.evt'],'w');
%     for i=hely:size(ripples,1)
%         fprintf (fid,'%i ripple-beg\n',ripples(i,1));
%         fprintf (fid,'%i ripple-peak\n',ripples(i,2));
%         fprintf (fid,'%i ripple-end\n',ripples(i,3));
%     end
% 
%     fclose(fid);