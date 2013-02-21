filename='AB3-58';

           samplerate= 1250; %S/s
            cutofflowpass=250; %Hz
            cutoffhighpass=100; %Hz
    
load('allchannels.mat'); %from 0

%read spw segments boundaries
load('swr.mat');

   %Extract LFP channel where the soma is
   data=memmapfile([filename '.eeg'],'Format','int16');
power = zeros(512,size(swr,2));
%take spw segments
    
    for j=1:size(swr,2)
        disp(['SWR ' int2str(j) ' of ' int2str(size(swr,2))]);
        %j=4;
        beg=swr(2,j)-50; %-100 ms
        fin=swr(2,j)+50; %+100 ms
        LFP=reshape(data.data(floor(beg*samplerate/1000)*512+1:floor(fin*samplerate/1000)*512),512,[]);
        %make a fourier on them
        %whiten
        for i=1:512
 
            %Filter channel
            fNorm = cutoffhighpass / (samplerate /2); 
            [d,e] = butter(4, fNorm, 'high');
            temp= filtfilt(d, e, double(LFP(i,:)));  
            
    
            fNorm = cutofflowpass / (samplerate /2); 
            [d,e] = butter(4, fNorm, 'low');
            Filtchan= filtfilt(d, e, temp);  
            clear temp d e fNorm
    
            %calculate power in x ms windows
            power(i,j)=sqrt(mean(Filtchan.^2));
        
        end
        
    end
clear data
   
powermean=mean(power,2);
powermeanreshaped=powermean(allchannels+1);
imagesc(atan([powermeanreshaped(:,1:8) (zeros(32,1)+min(min(powermeanreshaped))) powermeanreshaped(:,9:16)].*0.001));
imagesc(([powermeanreshaped(:,1:8) (zeros(32,1)+min(min(powermeanreshaped))) powermeanreshaped(:,9:16)]));
    print('-dtiff','-r300','-zbuffer',['Ripplepower' int2str(cutoffhighpass) '-' int2str(cutofflowpass) '.tiff']);
    print('-djpeg','-r300','-zbuffer',['Ripplepower' int2str(cutoffhighpass) '-' int2str(cutofflowpass) '.jpg']);
