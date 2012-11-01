function MakeWhl(fileBase,VideoLensBool,NumWhlLed,invertWhlYBool,invertSyncBool,nSamplesPerScreen,videoRes,ledDist)

% function MakeWhl(fileBase,VideoLensBool,NumWhlLed,invertWhlYBool,invertSyncBool,nSamplesPerScreen,videoRes,ledDist)
%
% Synchronize electrophysiological data, video streams and behavioral events
%
% Inputs
%  A .led file, a .spots file, and a .maz file
%
% Outputs
%  A .whl file , .whlLed and a .evt file
%
% Required parameters
%  fileBase: the (common) base name for the three input files
%
% Optional parameters
%  nSamplesPerScreen: the number of samples to display at once (default = 100)
%
% videoRes(x y): number of pixels in video input (default = [368, 240])
% videoRes(x y): number of pixels in video input (default = [736, 480])
%
% invertSyncBool may be used to correct the polarity of the syncronization pulse
%
% Note: this is an interactive function; follow the instructions on the main window

% This is a modified version of MakeWhlFile_3spots. This uses function CleanWhlFR instead of CleanWhl.
% CleanWhlFR is a modified version of CleanWhl. If the discrimination of Front and 
% Rear LED is imcomplete in MakeWhlFile_3spots program, uncleaned Whl file has 
% -1 -1 in only either whl(:,1:2) or whl(:,3:4). If those Whl file are used in Original 
% version of CleanWhl, the resultant cWhl file use -1 in whl(:,3:4) for interpolation.
% To prevent this problem, the value of whl(:,3:4) is also took into account for the rang of interpolation.

% If the original video file(.m1v) and spot-extracted file(.spots) are inverted in terms of Y-axis, invertWhlYBool should be 1, otherwise [].

% if you have Led to detect the running in the wheel, use "whlLedBool". A .whlLed file will be saved. This file has two columns. The first column has the number of pixels of Wheel Led Right, and the second column has the number of pixels of Wheel Led Left. The timestapms are the same as those of a conventional .whl file. Namely, size(whlLed,1) and size(whl,1) are the same.

% Sometimes the sync signal in the .eeg file is inverted. In that case use "invertSyncBool" should be 1, otherwise [].

% "nSamplesPerScreen" is used to show your data in the figure.

% "videoRes" should be [x y].

% VideolenseBool: 0 if you do not use video lens to magnify the visual range.
%						: 1 if you use video lens to magnify the visual range.
% This Bool is used to call the funtion to conpensate the XY axis ration, distortion caused by the lens, and finally compensate the pixel/centimeter so that it should be 1.

fileBase

if nargin<1
    fprintf('Usage: MakeWhl(fileBase,VideoLensBool,NumWhlLed,invertWhlYBool,invertSyncBool,nSamplesPerScreen,videoRes,ledDist)\n');
    return;
end

if nargin<2 |( VideoLensBool ~= 0 & VideoLensBool ~= 1);
	fprintf('Usage: MakeWhl(fileBase,VideoLensBool,NumWhlLed,invertWhlYBool,invertSyncBool,nSamplesPerScreen,videoRes,ledDist)\n');
	return;
end



if ~exist('nSamplesPerScreen','var') | isempty(nSamplesPerScreen),
  nSamplesPerScreen = 100;
end

%%% Original VideoRes was [368,240]. This value should be changed depending on your camera setting.

%%% The video resolution is [720x480].
if ~exist('videoRes','var') | isempty(videoRes),
  videoRes = [720,480];
end

if  ~exist('invertSyncBool','var') | isempty(invertSyncBool),
	    invertSyncBool = 0;
end

if  ~exist('invertWhlYBool','var') | isempty(invertWhlYBool),
	    invertWhlYBool = 1;
end

if  ~exist('NumWhlLed','var') | isempty(NumWhlLed),
	    NumWhlLed = 0;
end

if  ~exist('ledDist','var') | isempty(ledDist),
	    ledDist = 20;
end



%%% Thresh = 3; I have changed this value because I remove spots(:,2)<=10~40;
Thresh = 60; % to avoid triggers caused by 1-50 spurious pixels. (used in SchmittTrigger for video sync.)
ThreshWhlR = 30; % to avoid triggers caused by 1-50 spurious pixels. (used in SchmittTrigger for whl Led detection)
ThreshWhlL = 30; % to avoid triggers caused by 1-50 spurious pixels. (used in SchmittTrigger for whl Led detection)

% Load .spots, .led and .maz files

fprintf('Loading data ... ');
spots = [];
if exist([fileBase '.spots']),
   spots = load([fileBase '.spots']);
   if invertWhlYBool,
      spots(:,4) = -spots(:,4)+videoRes(2);
   end
   nFrames = max(spots(:,1))+1;
%%% Remove the tiny spots, which are usually not LEDs spots but junk spots.
%%% spots(:,2) represents how many pixels the spot has.
%%% You might want to use the restriction under which spots(:,5) > 0.5, spots(:,6) > 0.5, spots(:,7)is not an integer, spots(:,8) is not an integer, spots(:,9) is not an integer.

  figure(1);
  subplot(1,2,1);
  plot(spots(:,2),spots(:,7),'.','markersize',4,'linestyle','none','color',[0 0 0]);
  xlabel('spots(:,2), number of pixels');
  ylabel('spots(;,7), luminance');
  title('original spots file');
  subplot(1,2,2);
  plot(spots(:,3),spots(:,4),'.','markersize',4,'linestyle','none','color',[0 0 0]);
  xlabel('spots(:,3), x-axis in place');
  ylabel('spots(;,4), y-axis in place');
  title('original spots file');

  InPixThresh = 0;
  InLumThresh = 0;
  oriPixThresh = min(spots(:,2));
  oriLumThresh = min(spots(:,7));
  PixThresh = oriPixThresh;
  LumThresh = oriLumThresh;

  while length(InPixThresh)>0 | length(InLumThresh)>0;

    figure(2);
    subplot(1,2,1);
    plot(spots(:,2),spots(:,7),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(:,2), number of pixels');
    ylabel('spots(;,7), luminance');
    title('current spots file. \n Determine the pixels and luminance threshold');
    subplot(1,2,2);
    plot(spots(:,3),spots(:,4),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(:,3), x-axis in place');
    ylabel('spots(;,4), y-axis in place');
    title('current spots file');
    
    fprintf('\n\n Original threshold is %d for pixels and %d for luminance.\n', oriPixThresh, oriLumThresh);
    fprintf('\n\n Current threshold is %d for pixels and %d for luminance.\n', PixThresh, LumThresh);
    
    InPixThresh = input('\n\n Enter the number of pixels threshold and hit "enter". If the current pixels threshold is OK, hit "enter".','s');
    InLumThresh = input('\n\n Enter the number of luminances threshold and hit "enter". If the current luminance threshold is OK, hit "enter".','s');
    if length(str2num(InPixThresh))>0;
      newPixThresh = str2num(InPixThresh);
      if newPixThresh > PixThresh;
         PixThresh = newPixThresh;
      end
    end
    if length(str2num(InLumThresh))>0;
      newLumThresh = str2num(InLumThresh);
      if newLumThresh > LumThresh;
         LumThresh = newLumThresh;
      end
    end
    spots = spots(find(spots(:,2)>=PixThresh & spots(:,7)>=LumThresh),:);
  end
end

eegAll = bload([fileBase '.led' ], [1 inf]);
if invertSyncBool %invert Sync pulses
  eeg = -eegAll(1,:);
else
  eeg = eegAll(1,:);
end

if exist([fileBase '.maz']),
  events = load([fileBase '.maz']);
end
fprintf('Done\n');

% find upstrokes of eeg square wave - use mean+-1/4s.d. as trigger points
%%% syncEEG is the timestapm at which time the signal goes up beyound the UP-trigger after the signal goes down below the DOWN-trigger.
MeanEeg = mean(eeg); StdEeg = std(eeg);
syncEEG = SchmittTrigger(eeg,MeanEeg+StdEeg/4,MeanEeg-StdEeg/4);
syncEEGDown = SchmittTriggerDown(eeg,MeanEeg+StdEeg/4,MeanEeg-StdEeg/4);
fprintf('The first EEG pulse was detected at %i sec.\n', syncEEG(1)/1250);
fprintf('The second EEG pulse was detected at %i sec.\n', syncEEG(2)/1250);


if ~isempty(spots),
   figure(3);
   
   plot(spots(:,3),spots(:,4),'.','markersize',10);
   title('position of your LEDs');
   zoom on

   remove_more=1;
   while(remove_more==1)
      keyin = input('\n\n Do you want to remove some spots from analysis (yes/no)? ', 's');
        if strcmp(keyin,'yes'),
            input('In the figure window, select the area you want to delete.\n   Then click back in this window and hit ENTER...','s');
            xr = xlim;
            yr = ylim;
            [m, n] = size(spots);
            SpotsToKeep = find(~(spots(:,3)>=xr(1) & spots(:,3)<=xr(2) & spots(:,4)>=yr(1) & spots(:,4)<=yr(2)));
            spots = spots(SpotsToKeep,:);
            plot(spots(:,3),spots(:,4),'.','markersize',10);
        else
            if strcmp(keyin,'no'),
                remove_more=0;
            end
        end
   end
  
    if NumWhlLed == 1 | NumWhlLed == 2;

        fprintf('\nDetection of whl running LED(R)\n------------------------------\n');
        input('   In the figure window, select the area of the whl running LED(R) spots so it fills the axes.\n   Then click back in this window and hit ENTER...','s');
        title('select the area of the whl running Led spots');
        xr = xlim;
        yr = ylim;

        IsWhlLedR = spots(:,3)>=xr(1) & spots(:,3)<=xr(2) & spots(:,4)>=yr(1) & spots(:,4)<=yr(2);
        WhlLedR = find(IsWhlLedR); % entries in spots.txt for our LED
        PixWhlLedR = Accumulate(spots(WhlLedR,1)+1,spots(WhlLedR,2),nFrames); %%% PixWhlLedR is N x 1 matrix, and PixWhlLedR(:,1) has the number of Sync pixels in the same frame of the video. 
     
        syncWhlR = SchmittTrigger(PixWhlLedR,ThreshWhlR,0);
        syncWhlRDown = SchmittTriggerDown(PixWhlLedR,ThreshWhlR,0);
     
        SpotsToKeep = find(~IsWhlLedR);
        spots = spots(SpotsToKeep,:);
        plot(spots(:,3),spots(:,4),'.','markersize',10);

	
        ThreshPlotWhlR = ThreshWhlR.*ones(1,nFrames);
        %ThreshPlotWhlL = ThreshWhlL.*ones(1,nFrames);

        i=1;
        while i
            if strcmp(input('  If you correctly chose whl Led area, ''done''+ENTER to proceed...','s'),'done'),
                i=0;
            end
        end
	
        wLed = PixWhlLedR;

    end

    if NumWhlLed == 2;
   
        fprintf('\nDetection of whl running LED(L)\n------------------------------\n');
        input('   In the figure window, select the area of the whl running LED(L) spots so it fills the axes.\n   Then click back in this window and hit ENTER...','s');
        xr = xlim;
        yr = ylim;

        IsWhlLedL = spots(:,3)>=xr(1) & spots(:,3)<=xr(2) & spots(:,4)>=yr(1) & spots(:,4)<=yr(2);
        WhlLedL = find(IsWhlLedL); % entries in spots.txt for our LED
        PixWhlLedL = Accumulate(spots(WhlLedL,1)+1,spots(WhlLedL,2),nFrames); %%% PixWhlLedR is N x 1 matrix, and PixWhlLedR(:,1) has the number of Sync pixels in the same frame of the video. 
     
        syncWhlL = SchmittTrigger(PixWhlLedL,ThreshWhlL,0);
        syncWhlLDown = SchmittTriggerDown(PixWhlLedL,ThreshWhlL,0);
     
        SpotsToKeep = find(~IsWhlLedL);
        spots = spots(SpotsToKeep,:);
        plot(spots(:,3),spots(:,4),'.','markersize',10);
     
   %ThreshPlotWhlR = ThreshWhlR.*ones(1,nFrames);
        ThreshPlotWhlL = ThreshWhlL.*ones(1,nFrames);

        while i
            if strcmp(input('  If you correctly chose whl Led area, ''done''+ENTER to proceed...','s'),'done'),
                i=0;
            end
        end
   % This function makes a .wled file, which has the same timestamp as a conventional .whl file. The column 1 has the number of pixels of Whl LED right, and the column 2 has the number of pixels of Whl LED left.
     
        wLed = [wLed PixWhlLedL]; % this matrix will be converted into .wled file by using interp1.

    end

    if NumWhlLed == 1;
        figure(4)
        nFigFrames = ceil(nFrames/300);
        i=0;
        while ~strcmp(input('   Hit ENTER to show the next 10 seconds Video samples, or type ''done''+ENTER to proceed...','s'),'done'),
            i = i+1;
            if i > ceil(nFrames/300);
                i = 1;
                title('Repeating');
            end;
            first = (i-1)*300+1;
            last = min(i*300,nFrames);

            plot([first:last],PixWhlLedR(first:last),'color',[0 0 0]);hold on;
            plot([first:last],ThreshPlotWhlR(first:last),'color',[1 0 0]);
            syncWhlRPlot = syncWhlR(find(syncWhlR>=first & syncWhlR<=last));
            syncWhlRDownPlot = syncWhlRDown(find(syncWhlRDown>=first & syncWhlRDown<=last));
            if ~isempty(syncWhlRPlot),
                plot(syncWhlRPlot,ThreshWhlR,'o','color',[0 1 0]);
            end
            if ~isempty(syncWhlRDownPlot),
                plot(syncWhlRDownPlot,ThreshWhlR,'o','color',[0 0 1]);
            end
            xlabel('timestamps in .spots file unit');
            ylabel('number of pixels');
            title(['whl Led Right' ' ' num2str(i) '/' num2str(nFigFrames)]);
            hold off;
        end
    end

    if NumWhlLed == 2;
        figure(4)
        nFigFrames = ceil(nFrames/300);
        i=0;
        while ~strcmp(input('   Hit ENTER to show the next 10 seconds Video samples, or type ''done''+ENTER to proceed...','s'),'done'),
            i = i+1;
            if i > ceil(nFrames/300);
                i = 1;
                title('Repeating');
            end;
                first = (i-1)*300+1;
                last = min(i*300,nFrames);
	
                subplot(2,1,1),
                plot([first:last],PixWhlLedR(first:last),'color',[0 0 0]);hold on;
                plot([first:last],ThreshPlotWhlR(first:last),'color',[1 0 0]);
                syncWhlRPlot = syncWhlR(find(syncWhlR>=first & syncWhlR<=last));
                syncWhlRDownPlot = syncWhlRDown(find(syncWhlRDown>=first & syncWhlRDown<=last));
                if ~isempty(syncWhlRPlot),
                    plot(syncWhlRPlot,ThreshWhlR,'o','color',[0 1 0]);
                end
                if ~isempty(syncWhlRDownPlot),
                    plot(syncWhlRDownPlot,ThreshWhlR,'o','color',[0 0 1]);
                end
                xlabel('timestamps in .spots file unit');
                ylabel('number of pixels');
                title(['whl Led Right' ' ' num2str(i) '/' num2str(nFigFrames)]);
                hold off;
	
                subplot(2,1,2),
                plot([first:last],PixWhlLedL(first:last),'color',[0 0 0]);hold on;
                plot([first:last],ThreshPlotWhlL(first:last),'color',[1 0 0]);
                syncWhlLPlot = syncWhlL(find(syncWhlL>=first & syncWhlL<=last));
                syncWhlLDownPlot = syncWhlLDown(find(syncWhlLDown>=first & syncWhlLDown<=last));
                if ~isempty(syncWhlLPlot),
                    plot(syncWhlLPlot,ThreshWhlL,'o','color',[0 1 0]);
                end
                if ~isempty(syncWhlLDownPlot),
                    plot(syncWhlLDownPlot,ThreshWhlL,'o','color',[0 0 1]);
                end
                xlabel('timestamps in .spots file unit');
                ylabel('number of pixels');
                title(['whl Led Left' ' ' num2str(i) '/' num2str(nFigFrames)]);
                hold off;
	
        end
    end
     
     
    figure(5)  
    plot(spots(:,3),spots(:,4),'.','markersize',10);
    title('select the area of the sync spots');
  
  
    fprintf('\nDetection of synchronizing LED\n------------------------------\n');
    input('   In the figure window, select the area of the sync spot so it fills the axes.\n   Then click back in this window and hit ENTER...','s');
    xr = xlim;
    yr = ylim;

    IsSyncSpot = spots(:,3)>=xr(1) & spots(:,3)<=xr(2) & spots(:,4)>=yr(1) & spots(:,4)<=yr(2);
    SyncSpots = find(IsSyncSpot); % entries in spots.txt for our LED %%% SyncSpots is the index of spots file for the sycn spots.

	figure(6);
	plot(spots(SyncSpots,2),spots(SyncSpots,7),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(SyncSpots,2), number of pixels');
    ylabel('spots(SyncSpots,7), luminance');
    title('SyncSpots');

    PixCnt = Accumulate(spots(SyncSpots,1)+1,spots(SyncSpots,2),nFrames); %%% PixCnt is N x 1 matrix, and PixCnt(:,1) has the number of Sync pixels in the same frame of the video. 

  % now we get down to synchronization
  % - find flash points using "Schmitt trigger"
  % i.e. when it goes from 0 to above threshold
    syncVideo = SchmittTrigger(PixCnt,Thresh,0);
    syncVideoDown = SchmittTriggerDown(PixCnt,Thresh,0);
  
  %%% make sure that your thresholds for eeg and video are OK. 
    MeanEegPlot = MeanEeg.*ones(1,length(eeg));
    StdEegPlotUP = MeanEegPlot+StdEeg/4;
    StdEegPlotDown = MeanEegPlot-StdEeg/4;

    fprintf('   There are %d video flashes Up and %d square wave pulses Up.\n',length(syncVideo),length(syncEEG));

    nFigFrames = ceil(length(eeg)/(125000/2));
    figure(7);
    i=0;
    while ~strcmp(input('   Hit ENTER to show the next 50 seconds EEG samples, or type ''done''+ENTER to proceed...','s'),'done'),
        clf;
        i = i+1;
        if i > ceil(length(eeg)/(125000/2));
            i = 1;
            title('Repeating');
        end;
        first = (i-1)*(125000/2)+1;
        last = min(i*(125000/2),length(eeg));
        plot([first:last],eeg(first:last),'color',[0 0 0]); hold on;
        plot([first:last],MeanEegPlot(first:last),'color',[1 0 0]);
        plot([first:last],StdEegPlotUP(first:last),'color',[0 1 0]);
        plot([first:last],StdEegPlotDown(first:last),'color',[0 1 0]);
        syncEEGPlot = syncEEG(find(syncEEG>=first & syncEEG<=last));
         syncEEGDownPlot = syncEEGDown(find(syncEEGDown>=first & syncEEGDown<=last));
		if length(syncEEGPlot)>0;
            plot(syncEEGPlot, MeanEeg+StdEeg/4,'o','color',[0 0 1]);
		end
		if length(syncEEGDownPlot)>0;
            plot(syncEEGDownPlot, MeanEeg-StdEeg/4,'o','color',[1 0 1]);
		end
		hold off;
        title(['eeg' ' ' num2str(i) '/' num2str(nFigFrames)]);
        xlabel('timestapms in eeg unit');
        ylabel('voltage');
    end
% plot the last part of sync eeg.
	clf;
	i = ceil(length(eeg)/(125000/2));
	 first = (i-1)*(125000/2)+1;
     last = min(i*(125000/2),length(eeg));
     plot([first:last],eeg(first:last),'color',[0 0 0]); hold on;
     plot([first:last],MeanEegPlot(first:last),'color',[1 0 0]);
     plot([first:last],StdEegPlotUP(first:last),'color',[0 1 0]);
     plot([first:last],StdEegPlotDown(first:last),'color',[0 1 0]);
     syncEEGPlot = syncEEG(find(syncEEG>=first & syncEEG<=last));
     syncEEGDownPlot = syncEEGDown(find(syncEEGDown>=first & syncEEGDown<=last));
    if length(syncEEGPlot)>0;
        plot(syncEEGPlot, MeanEeg+StdEeg/4,'o','color',[0 0 1]);
    end
    if length(syncEEGDownPlot)>0;
        plot(syncEEGDownPlot, MeanEeg-StdEeg/4,'o','color',[1 0 1]);
	end
	hold off;
    title(['eeg' ' ' num2str(i) '/' num2str(nFigFrames)]);
    xlabel('timestapms in eeg unit');
    ylabel('voltage');
  
    nFigFrames = ceil(nFrames/1500);
    ThreshPlot = Thresh.*ones(1,nFrames);
    figure(8);
    i=0;
     while ~strcmp(input('   Hit ENTER to show the next 50 seconds Video samples, or type ''done''+ENTER to proceed...','s'),'done'),
        clf;
        i = i+1;
        if i > ceil(nFrames/1500);
            i = 1;
            title('Repeating');
        end;
        first = (i-1)*1500+1;
        last = min(i*1500,nFrames);
        plot([first:last],PixCnt(first:last),'color',[0 0 0]);hold on;
        plot([first:last],ThreshPlot(first:last),'color',[1 0 0]);
        syncVideoPlot = syncVideo(find(syncVideo>=first & syncVideo<=last));
        syncVideoDownPlot = syncVideoDown(find(syncVideoDown>=first & syncVideoDown<=last));
		if length(syncVideoPlot)>0;
            plot(syncVideoPlot,Thresh,'o','color',[0 1 0]);
		end
		if length(syncVideoDownPlot)>0;
            plot(syncVideoDownPlot,Thresh,'o','color',[0 0 1]);
		end
		hold off;
        title(['video sync led' ' ' num2str(i) '/' num2str(nFigFrames)]);
        xlabel('timestapms in .spots unit');
        ylabel('pixels');
    end
% plot the last part of the sync video.
	clf;
	i = ceil(nFrames/1500);
	first = (i-1)*1500+1;
    last = min(i*1500,nFrames);
    plot([first:last],PixCnt(first:last),'color',[0 0 0]);hold on;
    plot([first:last],ThreshPlot(first:last),'color',[1 0 0]);
    syncVideoPlot = syncVideo(find(syncVideo>=first & syncVideo<=last));
    syncVideoDownPlot = syncVideoDown(find(syncVideoDown>=first & syncVideoDown<=last));
	if length(syncVideoPlot)>0;
	   	plot(syncVideoPlot,Thresh,'o','color',[0 1 0]);
	end
	if length(syncVideoDownPlot)>0;
	    plot(syncVideoDownPlot,Thresh,'o','color',[0 0 1]);
	end
	hold off;
    title(['video sync led' ' ' num2str(i) '/' num2str(nFigFrames)]);
    xlabel('timestapms in .spots unit');
    ylabel('pixels');


  % Before manual correction:
  %  syncEEG contains the timestamps of the SYNC upstrokes detected in the EEG
  %  syncVideo contains the timestamps of the LED onsets detected in the video
  % After manual correction:
  %  syncEEG and syncVideo will contain the corrected timestamps
  %  syncVideo will be truncated if it is longer than syncEEG, but not the other way around
  %  Therefore, length(syncVideo) <= length(syncEEG)
    fprintf('\nSynchronization of video and EEG \n--------------------------------\n');
  %%%% syncEEG(length(syncVideo)+1:end) = []; %%%% 061124 tentatively.
    while 1,
        fprintf('   There are %d video flashes Up and %d square wave pulses Up.\n',length(syncVideo),length(syncEEG));
        if length(syncVideo)~=length(syncEEG),
            fprintf('   ***** MISMATCH! *****\n');
            i = input('   To manually correct this, hit ENTER; to drop flashes in excess and continue, type ''done''+ENTER - but you''d better be sure...','s');
        else
            i = input('   To manually edit this, hit ENTER; to continue, type ''done''+ENTER...','s');
        end
        if strcmp(i,'done'),
            if length(syncVideo) > length(syncEEG),
                syncVideo = syncVideo(1:length(syncEEG));
            end
%%% FrameSamples = interp1(syncVideo,syncEEG(1:length(syncVideo)),1:nFrames,'linear',NaN); if the length of syncEEG is longer than syncVideo, syncEEG is truncated.
            break;
        else
            [syncVideo,syncEEG] = DisplaySYNCnFig(syncVideo,syncEEG,nSamplesPerScreen,9);
            title('timestamps of video UP(bottom) and EEG UP(top)');
        end
    end
  
  %%%% syncEEGDown(length(syncVideoDown)+1:end)=[]; %%%% 061124 tentatively.
    while 1,
        fprintf('   There are %d video flashes Down and %d square wave pulses Down.\n',length(syncVideoDown),length(syncEEGDown));
        if length(syncVideoDown)~=length(syncEEGDown),
            fprintf('   ***** MISMATCH! *****\n');
            i = input('   To manually correct this, hit ENTER; to drop flashes in excess and continue, type ''done''+ENTER - but you''d better be sure...','s');
        else
            i = input('   To manually edit this, hit ENTER; to continue, type ''done''+ENTER...','s');
        end
        if strcmp(i,'done'),
            if length(syncVideoDown) > length(syncEEGDown),
                syncVideoDown = syncVideoDown(1:length(syncEEGDown));
            end
            break;
        else
        [syncVideoDown,syncEEGDown] = DisplaySYNCnFig(syncVideoDown,syncEEGDown,nSamplesPerScreen,10);
        title('timestamps of video Down(bottom) and EEG Donw(top)');
        end
    end

  
  	fprintf(' Check the UP and Down of video and eeg. \n')

    nSamplesPerScreenUD = ceil(nSamplesPerScreen/4);
    nSamplesPerScreenUD = min([nSamplesPerScreenUD length(syncVideo) length(syncEEG)]);
    ratioUD = (syncVideo(nSamplesPerScreenUD)-syncVideo(1))/(syncEEG(nSamplesPerScreenUD)-syncEEG(1));
    shiftUD = syncVideo(1) - syncEEG(1)*ratioUD;


    figure(11);hold on;
    xlabel('red is UP, blue is Down');
    ylabel('top is EEG, bottom is video');
    nSync1 = length(syncVideo);
    nSync1Down = length(syncVideoDown);
    nSync2 = length(syncEEG);
    nSync2Down = length(syncEEGDown);
    n = min([nSync1 nSync2]);
    nFigFrames = ceil(n/nSamplesPerScreenUD);
    i = 0;
    while ~strcmp(input('   Hit ENTER to show the next 25 samples, or type ''done''+ENTER to proceed...','s'),'done'),
         i = i+1;
       if i > ceil(n/nSamplesPerScreenUD), i=1; end
       hold on;
       first = (i-1)*nSamplesPerScreenUD+1;
       last1 = min([i*nSamplesPerScreenUD nSync1]);
       last1Down = min([i*nSamplesPerScreenUD nSync1Down]);
       plot(syncVideo(first:last1),zeros(last1-first+1,1),...
              '+','markersize',6,'linestyle','none','color',[1 0 0]);
       plot(syncVideoDown(first:last1Down),zeros(last1Down-first+1,1),...
	      'o','markersize',6,'linestyle','none','color',[0 0 1]);
	  
       last2 = min([i*nSamplesPerScreenUD nSync2]);
       last2Down = min([i*nSamplesPerScreenUD nSync2Down]);
       plot(syncEEG(first:last2)*ratioUD+shiftUD,0.1*ones(last2-first+1,1),'+','markersize',6,'linestyle','none','color',[1 0 0]);
       plot(syncEEGDown(first:last2Down)*ratioUD+shiftUD,0.1*ones(last2Down-first+1,1),'o','markersize',6,'linestyle','none','color',[0 0 1]);
       set(gca,'ylim',[-.025 .125],'xlim',[min([syncVideo(first) syncVideoDown(first) syncEEG(first)*ratioUD+shiftUD syncEEGDown(first)*ratioUD+shiftUD]) max([syncVideo(last1) syncVideoDown(last1Down) syncEEG(last2)*ratioUD+shiftUD syncEEGDown(last2Down)*ratioUD+shiftUD])]);hold off;
       lastUP = min([i*nSamplesPerScreen nSync1 nSync2]);
       lastDown = min([i*nSamplesPerScreen nSync1Down nSync2Down]);
       title(['Make sure that EEG and Video have the same ratio=(duration UP)/(duration Down)' ' ' num2str(i) '/' num2str(nFigFrames)]);
       for j = first:lastUP,
          line([syncVideo(j) syncEEG(j)*ratioUD+shiftUD],[0 0.1],'color',[1 0 0]);
       end
       for j = first:lastDown,
            line([syncVideoDown(j) syncEEGDown(j)*ratioUD+shiftUD],[0 0.1],'color',[0 0 1]);
       end
    end
end

if ~isempty(spots),
    I=1;
    while I==1;
        Keyin = input('\n\n Do you want to continue this analysis (yes/no)? ', 's');
			if strcmp(Keyin,'yes'),
	     		I=0;
            else
                if strcmp(Keyin,'no'),
                return;
			end
    end 
end

% Before manual correction:
%  syncEEG contains the timestamps of the SYNC upstrokes detected in the EEG
%  syncEvents contains the timestamps of the SYNC events in the .maz file
% After manual correction:
%  syncEEG and syncEvents will contain the corrected timestamps
%  syncEvents will be truncated if it is longer than syncEEG, but not the other way around
%  Therefore, length(syncEvents) <= length(syncEEG)
if ~exist('events', 'var')
    events = [];
end

if ~isempty(events),
  fprintf('\nSynchronization of events and EEG\n---------------------------------\n');
  syncEvents = events(find(events(:,2) == 83 & events(:,3) == 89));
  while 1,
      fprintf('   There are %d SYNC events and %d square wave pulses.\n',length(syncEvents),length(syncEEG));
      if length(syncEvents)~=length(syncEEG),
        fprintf('   ***** MISMATCH! *****\n');
        i = input('   To manually correct this, hit ENTER; to drop flashes in excess and continue, type ''done''+ENTER - but you''d better be sure...','s');
      else
        i = input('   To manually edit this, hit ENTER; to continue, type ''done''+ENTER...','s');
      end
      if strcmp(i,'done'),
          if length(syncEvents) > length(syncEEG),
              syncEvents = syncEvents(1:length(syncEEG));
          end
          break;
      else
        [syncEvents,syncEEG] = DisplaySYNC(syncEvents,syncEEG,nSamplesPerScreen);
      end
  end
end

if ~isempty(spots),
  figure(12);
  title('UP flash stat');
  subplot(2,2,1);
  plot(diff(syncVideo),'.','markersize',4);
  ylabel('# video frames between flashes');
  xlabel('UP flash #');

  subplot(2,2,2)
  plot(diff(syncEEG), '.','markersize',4);
  ylabel('# EEG samples between flashes');
  xlabel('UP flash #');

  subplot(2,2,3);
  [b bint r] = regress(syncEEG(1:length(syncVideo)),[syncVideo,ones(size(syncVideo))]);
  plot(r/1.25,'.','markersize',4);
  xlabel('UP Flash #');
  ylabel('deviation from linear fit (ms)');

  subplot(2,2,4);
  plot(diff(syncVideo)./diff(syncEEG(1:length(syncVideo)))*1250,'.','markersize',4);
%  FilterLen = 10;
%  f = filter(ones(FilterLen,1)/FilterLen, 1,diff(syncVideo)./diff(syncEEG(1:length(syncVideo)))*1250);
%  plot(FilterLen:length(f),f(FilterLen:end),'r');
%  ylabel('Frame rate (red is smoothed)');
	 ylabel('Number of video frames per sec');
  xlabel('UP Flash #');
  
  ReportFigKM({[fileBase '_Make_Whl_File_' date]},{[fileBase '_Up_frame']},[],gcf);

  
  figure(13);
  title('Down flash stat');
  subplot(2,2,1);
  plot(diff(syncVideoDown),'.','markersize',4);
  ylabel('# video frames between flashes');
  xlabel('Down flash #');

  subplot(2,2,2)
  plot(diff(syncEEGDown), '.','markersize',4);
  ylabel('# EEG samples between flashes');
  xlabel('Down flash #');

  subplot(2,2,3);
  [b bint r] = regress(syncEEGDown(1:length(syncVideoDown)),[syncVideoDown,ones(size(syncVideoDown))]);
  plot(r/1.25,'.','markersize',4);
  xlabel('Down Flash #');
  ylabel('deviation from linear fit (ms)');

  subplot(2,2,4);
  plot(diff(syncVideoDown)./diff(syncEEGDown(1:length(syncVideoDown)))*1250,'.','markersize',4);
%  FilterLen = 10;
% f = filter(ones(FilterLen,1)/FilterLen, 1,diff(syncVideoDown)./diff(syncEEGDown(1:length(syncVideoDown)))*1250);
%  plot(FilterLen:length(f),f(FilterLen:end),'r');
  ylabel('Number of video frames per sec');
  xlabel('Down Flash #');

   ReportFigKM({[fileBase '_Make_Whl_File_' date]},{[fileBase '_Down_frame']},[],gcf);
  
  
  DurationVideoUP=[];
  DurationVideoDown=[];
  DurationEEGUP=[];
  DurationEEGDown=[];
  N_UP = min([nSync1 nSync1Down nSync2 nSync2Down])-1;
  for i=1:N_UP
    DurationVideoUP=[DurationVideoUP,syncVideoDown(i)-syncVideo(i)];
    DurationVideoDown=[DurationVideoDown,syncVideo(i+1)-syncVideoDown(i)];
    DurationEEGUP=[DurationEEGUP,syncEEGDown(i)-syncEEG(i)];
    DurationEEGDown=[DurationEEGDown,syncEEG(i+1)-syncEEGDown(i)];
  end
  figure(14);
 
  subplot(2,2,1)
  plot((DurationEEGUP./DurationVideoUP).*ratioUD,'markersize',4);
  title('DurationEEGUP./DurationVideoUP');
  
  subplot(2,2,2)
  plot((DurationEEGDown./DurationVideoDown).*ratioUD,'markersize',4);
  title('DurationEEGDown./DurationVideoDown');
  
  subplot(2,2,3)
  plot(DurationEEGUP./DurationEEGDown,'markersize',4);
  title('DurationEEGUP./DurationEEGDown');
  
  subplot(2,2,4)
  plot(DurationVideoUP./DurationVideoDown,'markersize',4);
  title('DurationVidoeUP./DurationVideoDown');

   ReportFigKM({[fileBase '_Make_Whl_File_' date]},{[fileBase '_Up_DOwn_Duration']},[],gcf);

  % now align them - any outside sync range become NaN
  %%% FrameSamples has the EEG index, which correspond to the video index[1:nFrames]. size(FrameSamples) = [nFrames 1] (09/16/05 Kenji)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%% In order to compenstate the video file which has only short-time sync signals.
  FrameSamples = interp1(syncVideo,syncEEG(1:length(syncVideo)),1:nFrames,'linear',NaN);
 %%%%%%%%%A = [syncVideo(1) syncVideo(end)];
 %%%%%%%%%%B = [syncEEG(1) syncEEG(end)];
%%%%%%%%%%%% FrameSamples = interp1(A,B,1:nFrames,'linear','extrap'); % In order to compenstate the video file which has only short-time sync signals.
 %%%%%%%%%%%%FrameSamples(1:300) = NaN;
 %%%%%%%%%%%%FrameSamples(end-300:end) = NaN;
  %%%%%%%%%%%%%%%%%%%%%%%%%%% In order to compenstate the video file which has only short-time sync signals.
  
  % THE FOLLOWING ONLY WORKS IF YOU HAVE 2 LEDS

  % find non-sync spots
  NonSyncSpots = find(~IsSyncSpot);


%%%%%%%%%%%% If you want to use only spots with 2 non=sync spots now, do it. Otherwise do it after separation the spots depending on the RGB color cordinate.

  % find those with 2 non-sync spots per frame
	NonSyncCnt = Accumulate(1+spots(NonSyncSpots,1),1);
	Spots_1 = NonSyncSpots(find(NonSyncCnt(1+spots(NonSyncSpots,1))==1));
	Spots_2 = NonSyncSpots(find(NonSyncCnt(1+spots(NonSyncSpots,1))==2));
	Spots_3 = NonSyncSpots(find(NonSyncCnt(1+spots(NonSyncSpots,1))>=3));

	Spots_1_Ratio = length(Spots_1)./length(NonSyncCnt),
	Spots_2_Ratio = length(Spots_2)./(length(NonSyncCnt).*2),
	Spots_3_Ratio = length(Spots_3)./(length(NonSyncCnt).*3), %% to be more accurate, the length(NonSyncCnt) should be multiplied by the number of the spots in the same frame.

	%% The frames which contain more than two Front and Rear LED should be removed after the color (front vs rear) separation. (060328)
  %%%  GoodSpots is the index for spots file.
  % select front versus rear LED in color space
  
  %% separated clusters of points corresponding to different LED-colors (in green/blue space!!!)
  %% fprintf('\n\n\n Select a cluster that belongs to one color.\n If the color is selected, double-click within the figure.\n');
 
%%%%%%%%%%%% If you want to use only spots with 2 non=sync spots now, do it. Otherwise do it after separation the spots depending on the RGB color cordinate.


	%% rgb = ycbcr2rgb(spots(GoodSpots,7:9));
	rgb_1 = ycbcr2rgb_K(spots(Spots_1,7:9));
	rgb_2 = ycbcr2rgb_K(spots(Spots_2,7:9));
    if ~isempty(Spots_3);
        rgb_3 = ycbcr2rgb_K(spots(Spots_3,7:9));
    end
  
  	RGB_1 = [rgb_1(:,1,1) rgb_1(:,1,2) rgb_1(:,1,3)];
	RGB_2 = [rgb_2(:,1,1) rgb_2(:,1,2) rgb_2(:,1,3)];
    if ~isempty(Spots_3);
        RGB_3 = [rgb_3(:,1,1) rgb_3(:,1,2) rgb_3(:,1,3)];
    end
    
  [U1 S1 V1] = svd(RGB_1,0);
  %[U1 S1 V1] = svd(rgb_1,0);
  PCArgb_1 = U1*S1;

	[U2 S2 V2] = svd(RGB_2,0);
    %[U2 S2 V2] = svd(rgb_2,0);
  	PCArgb_2 = U2*S2;
    
    if ~isempty(Spots_3);
        [U3 S3 V3] = svd(RGB_3,0);
        %[U3 S3 V3] = svd(rgb_3,0);
        PCArgb_3 = U3*S3;
    end
    
  figure(15);
  
  subplot(2,3,1);
  plot(rgb_1(:,1), rgb_1(:,2),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('1 for One-spot in a frame');
  xlabel('Red');
  ylabel('Green');hold off;
  
  subplot(2,3,2);
  plot(rgb_1(:,1), rgb_1(:,3),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('2 for One-spot in a frame');
  xlabel('Red');
  ylabel('Blue');hold off;
  
  subplot(2,3,3);
  plot(rgb_1(:,2), rgb_1(:,3),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('3 for One-spot in a frame');
  xlabel('Green');
  ylabel('Blue');hold off;
  
  subplot(2,3,4);
  plot(PCArgb_1(:,1),PCArgb_1(:,2),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('4 for One-spot in a frame');
  xlabel('PC1');
  ylabel('PC2');
  
  subplot(2,3,5);
  plot(PCArgb_1(:,1),PCArgb_1(:,3),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('5 for One-spot in a frame');
  xlabel('PC1');
  ylabel('PC3');
  
  subplot(2,3,6);
  plot(PCArgb_1(:,2),PCArgb_1(:,3),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('6 for One-spot in a frame');
  xlabel('PC2');
  ylabel('PC3');

	 ReportFigKM({[fileBase '_Make_Whl_File_' date]},{[fileBase '_color separation for 1 spot per frame']},[],gcf);
  
figure(16);
  
  subplot(2,3,1);
  plot(rgb_2(:,1), rgb_2(:,2),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('1 for Two-spots in a frame');
  xlabel('Red');
  ylabel('Green');hold off;
  
  subplot(2,3,2);
  plot(rgb_2(:,1), rgb_2(:,3),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('2 for Two-spots in a frame ');
  xlabel('Red');
  ylabel('Blue');hold off;
  
  subplot(2,3,3);
  plot(rgb_2(:,2), rgb_2(:,3),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('3 for Two-spots in a frame');
  xlabel('Green');
  ylabel('Blue');hold off;
  
  subplot(2,3,4);
  plot(PCArgb_2(:,1),PCArgb_2(:,2),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('4 for Two-spots in a frame');
  xlabel('PC1');
  ylabel('PC2');
  
  subplot(2,3,5);
  plot(PCArgb_2(:,1),PCArgb_2(:,3),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('5 for Two-spots in a frame ');
  xlabel('PC1');
  ylabel('PC3');
  
  subplot(2,3,6);
  plot(PCArgb_2(:,2),PCArgb_2(:,3),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('6 for Two-spots in a frame ');
  xlabel('PC2');
  ylabel('PC3');

	ReportFigKM({[fileBase '_Make_Whl_File_' date]},{[fileBase '_color separation for 2 spots per frame']},[],gcf);

if ~isempty(Spots_3);
figure(17);
  
  subplot(2,3,1);
  plot(rgb_3(:,1), rgb_3(:,2),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('1 for Three or more spots in a frame');
  xlabel('Red');
  ylabel('Green');hold off;
  
  subplot(2,3,2);
  plot(rgb_3(:,1), rgb_3(:,3),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('2 for Three or more spots in a frame');
  xlabel('Red');
  ylabel('Blue');hold off;
  
  subplot(2,3,3);
  plot(rgb_3(:,2), rgb_3(:,3),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('3 for Three or more spots in a frame');
  xlabel('Green');
  ylabel('Blue');hold off;
  
  subplot(2,3,4);
  plot(PCArgb_3(:,1),PCArgb_3(:,2),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('4 for Three or more spots in a frame');
  xlabel('PC1');
  ylabel('PC2');
  
  subplot(2,3,5);
  plot(PCArgb_3(:,1),PCArgb_3(:,3),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('5 for Three or more spots in a frame');
  xlabel('PC1');
  ylabel('PC3');
  
  subplot(2,3,6);
  plot(PCArgb_3(:,2),PCArgb_3(:,3),'.','color',[0 0 0],'markersize',2,'linestyle','none');hold on
  title('6 for Three or more spots in a frame');
  xlabel('PC2');
  ylabel('PC3');

	ReportFigKM({[fileBase '_Make_Whl_File_' date]},{[fileBase '_color separation for 3 or more spots per frame']},[],gcf);
end


  WhichProjection=1;
  while WhichProjection==1
      WP = input('\n\n Which projection do you want to use for discriminationg front and rear LEDs? Choose 1/2/3/4/5/6 from the title of the figure', 's');
  	if strcmp(WP,'1')
		WhichP = 1;
		WhichProjection=0;
	end
	if strcmp(WP,'2')
		WhichP = 2;
		WhichProjection=0;
	end
	if strcmp(WP,'3')
		WhichP = 3;
		WhichProjection=0;
	end
	if strcmp(WP,'4')
		WhichP = 4;
		WhichProjection=0;
	end
	if strcmp(WP,'5')
		WhichP = 5;
		WhichProjection=0;
	end
	if strcmp(WP,'6')
		WhichP = 6;
		WhichProjection=0;
	end
  end
end
    %% plot(rgb(:,2),rgb(:,3),'.','markersize',1);
    %% xlabel('green');
    %% ylabel('Blue');
    %% zoom on;
    %% fprintf('\nDiscrimination of front and rear lights\n---------------------------------------\n');
    %% input('   In the figure window, select the front spot so it fills the axes.\n   Then click back in this window and hit ENTER...','s');
    %% xr = xlim;
    %% yr = ylim;
    %% IsFrontSpot = rgb(:,2)>=xr(1) & rgb(:,2)<=xr(2) & rgb(:,3)>=yr(1) & rgb(:,3)<=yr(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('\n Spots_2. Draw boarders of a cluster by clicking into a figure. Do not click on the edge of the figure.');
  switch WhichP
     case 1

     figure(18);
     
     title('Two-Spots in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_2 = ClusterPoints(rgb_2(:,1:2),1);

     figure(19);
     
     title('Two-Spots in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_2 = ClusterPoints(rgb_2(:,1:2),1);
     case 2
     
     figure(18);
     
     title('Two-Spots in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_2 = ClusterPoints(rgb_2(:,[1 3]),1);

     figure(19);
     
     title('Two-Spots in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_2 = ClusterPoints(rgb_2(:,[1 3]),1);
     
     case 3
     
     figure(18);
     
     title('Two-Spots in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_2 = ClusterPoints(rgb_2(:,2:3),1);

     figure(19);
     
     title('Two-Spots in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_2 = ClusterPoints(rgb_2(:,2:3),1);
     
     case 4
     
     figure(18);
     
     title('Two-Spots in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_2 = ClusterPoints(PCArgb_2(:,1:2),1);

     figure(19);
     
     title('Two-Spots in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_2 = ClusterPoints(PCArgb_2(:,1:2),1);
     
     case 5
     
     figure(18);
     
     title('Two-Spots in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_2 = ClusterPoints(PCArgb_2(:,[1 3]),1);

     figure(19);
     
     title('Two-Spots in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_2 = ClusterPoints(PCArgb_2(:,[1 3]),1);
     
     case 6
     
     figure(18);
     
     title('Two-Spots in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_2 = ClusterPoints(PCArgb_2(:,2:3),1);

     figure(19);
     
     title('Two-Spots in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_2 = ClusterPoints(PCArgb_2(:,2:3),1);
     
  end
   
  figure(20);
  
  subplot(2,3,1);
  plot(rgb_2(:,1), rgb_2(:,2),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(rgb_2(find(IsFrontSpot_2),1),rgb_2(find(IsFrontSpot_2),2),'.','color',[1 0 0],'markersize',3,'linestyle','none');
  plot(rgb_2(find(IsRearSpot_2),1),rgb_2(find(IsRearSpot_2),2),'.','color',[0 0 1],'markersize',3,'linestyle','none');hold off;
  title('1');
  xlabel('Red');
  ylabel('Green');hold off;
  
  subplot(2,3,2);
  plot(rgb_2(:,1), rgb_2(:,3),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(rgb_2(find(IsFrontSpot_2),1),rgb_2(find(IsFrontSpot_2),3),'.','color',[1 0 0],'markersize',3,'linestyle','none');
  plot(rgb_2(find(IsRearSpot_2),1),rgb_2(find(IsRearSpot_2),3),'.','color',[0 0 1],'markersize',3,'linestyle','none');hold off;
  title('2');
  xlabel('Red');
  ylabel('Blue');hold off;
  
  subplot(2,3,3);
  plot(rgb_2(:,2), rgb_2(:,3),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(rgb_2(find(IsFrontSpot_2),2),rgb_2(find(IsFrontSpot_2),3),'.','color',[1 0 0],'markersize',3,'linestyle','none');
  plot(rgb_2(find(IsRearSpot_2),2),rgb_2(find(IsRearSpot_2),3),'.','color',[0 0 1],'markersize',3,'linestyle','none');hold off;
  title('3');
  xlabel('Green');
  ylabel('Blue');hold off;
  
  subplot(2,3,4);
  plot(PCArgb_2(:,1), PCArgb_2(:,2),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(PCArgb_2(find(IsFrontSpot_2),1),PCArgb_2(find(IsFrontSpot_2),2),'.','color',[1 0 0],'markersize',3,'linestyle','none');
  plot(PCArgb_2(find(IsRearSpot_2),1),PCArgb_2(find(IsRearSpot_2),2),'.','color',[0 0 1],'markersize',3,'linestyle','none');hold off;
  title('4');
  xlabel('PC1');
  ylabel('PC2');hold off;
  
  subplot(2,3,5);
  plot(PCArgb_2(:,1), PCArgb_2(:,3),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(PCArgb_2(find(IsFrontSpot_2),1),PCArgb_2(find(IsFrontSpot_2),3),'.','color',[1 0 0],'markersize',3,'linestyle','none');
  plot(PCArgb_2(find(IsRearSpot_2),1),PCArgb_2(find(IsRearSpot_2),3),'.','color',[0 0 1],'markersize',3,'linestyle','none');hold off;
  title('5');
  xlabel('PC1');
  ylabel('PC3');hold off;
  
  subplot(2,3,6);
  plot(PCArgb_2(:,2), PCArgb_2(:,3),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(PCArgb_2(find(IsFrontSpot_2),2),PCArgb_2(find(IsFrontSpot_2),3),'.','color',[1 0 0],'markersize',3,'linestyle','none');
  plot(PCArgb_2(find(IsRearSpot_2),2),PCArgb_2(find(IsRearSpot_2),3),'.','color',[0 0 1],'markersize',3,'linestyle','none');hold off;
  title('6');
  xlabel('PC2');
  ylabel('PC3');hold off;
  
  
  %% Get position information using the index of spot files.
  %% Sometimes the color assinment is opposite.

	OKcolor = input('\n Is the color assinment is OK (or opposite)? y/n','s');
  while 1;
    if strcmp(OKcolor,'y'),
      FrontSpots_2 = Spots_2(find(IsFrontSpot_2 & ~IsRearSpot_2));
      RearSpots_2 = Spots_2(find(IsRearSpot_2 & ~IsFrontSpot_2));
        WeirdSpots_2 = Spots_2(find(~IsFrontSpot_2 & ~IsRearSpot_2));
      break
    end
    if strcmp(OKcolor,'n'),
      FrontSpots_2 = Spots_2(find(IsRearSpot_2 & ~IsFrontSpot_2));
      RearSpots_2 = Spots_2(find(IsFrontSpot_2 & ~IsRearSpot_2));
		WeirdSpots_2 = Spots_2(find(~IsFrontSpot_2 & ~IsRearSpot_2));
      break
    end
  end

%%% NonSyncSpots, FrontSpots, RearSpots are the index of spots file.(060329)
%%% If one frame has two or more Frant or Rear LED spots, ignore the frame.
%%% You might wont to change the algorithm.

	FrontFrames_2 = Accumulate(1+spots(FrontSpots_2,1),1);
	GoodFrontSpots_2 = FrontSpots_2(find(FrontFrames_2(1+spots(FrontSpots_2,1))==1));
	GoodFrontSpotsRatio_2 = length(GoodFrontSpots_2)./length(FrontSpots_2)
	FrontSpots_2 = GoodFrontSpots_2;

	RearFrames_2 = Accumulate(1+spots(RearSpots_2,1),1);
	GoodRearSpots_2 = RearSpots_2(find(RearFrames_2(1+spots(RearSpots_2,1))==1));
	GoodRearSpotsRatio_2 = length(GoodRearSpots_2)./length(RearSpots_2)
	RearSpots_2 = GoodRearSpots_2;

	WiredSpotsRatio_2 = length(WeirdSpots_2)./length(Spots_2)

	figure(21);

    subplot(2,3,1);
    plot(spots(FrontSpots_2,2),spots(FrontSpots_2,7),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(FrontSpots_2,2), number of pixels');
    ylabel('spots(FrontSpots_2,7), luminance');
    title('FrontSpots_2');

	subplot(2,3,2);
    plot(spots(RearSpots_2,2),spots(RearSpots_2,7),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(RearSpots_2,2), number of pixels');
    ylabel('spots(RearSpots_2,7), luminance');
    title('RearSpots_2');

	subplot(2,3,3);
    plot(spots(WeirdSpots_2,2),spots(WeirdSpots_2,7),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(RearSpots_2,2), number of pixels');
    ylabel('spots(RearSpots_2,7), luminance');
    title('WeirdSpots_2');

	subplot(2,3,4);
    plot(spots(FrontSpots_2,3),spots(FrontSpots_2,4),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(FrontSpots_2,3), x');
    ylabel('spots(FrontSpots_2,4), y');
    title('FrontSpots_2');

	subplot(2,3,5);
    plot(spots(RearSpots_2,3),spots(RearSpots_2,4),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(RearSpots_2,3), x');
    ylabel('spots(RearSpots_2,4), y');
    title('WeirdSpots_2');

	subplot(2,3,6);
    plot(spots(WeirdSpots_2,3),spots(WeirdSpots_2,4),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(WeirdSpots_2,3), x');
    ylabel('spots(WeirdSpots_2,4), y');
    title('WeirdSpots_2');

  figure(22);
  plot(spots(FrontSpots_2,3),spots(FrontSpots_2,4),'.','color',[1 0 0],'markersize',10,'linestyle','none');hold on;
  plot(spots(RearSpots_2,3),spots(RearSpots_2,4),'.','color',[0 1 0],'markersize',10,'linestyle','none');hold off;
	title('Two spots in a frame, whl file');

%%%%%%%%%%%%%%%%%%%%%%%

 fprintf('\n Spots_1. Draw boarders of a cluster by clicking into a figure. Do not click on the edge of the figure.');
  switch WhichP
     case 1

     figure(23);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_1 = ClusterPointsMulti(rgb_1(:,1:2),rgb_2(find(IsFrontSpot_2),1:2),rgb_2(find(IsRearSpot_2),1:2),1,2);

     figure(24);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_1 = ClusterPointsMulti(rgb_1(:,1:2),rgb_2(find(IsFrontSpot_2),1:2),rgb_2(find(IsRearSpot_2),1:2),1,2);

     case 2
     
     figure(23);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_1 = ClusterPointsMulti(rgb_1(:,[1 3]),rgb_2(find(IsFrontSpot_2),[1 3]),rgb_2(find(IsRearSpot_2),[1 3]),1,2);

     figure(24);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_1 = ClusterPointsMulti(rgb_1(:,[1 3]),rgb_2(find(IsFrontSpot_2),[1 3]),rgb_2(find(IsRearSpot_2),[1 3]),1,2);
     
     case 3
     
     figure(23);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_1 = ClusterPointsMulti(rgb_1(:,2:3),rgb_2(find(IsFrontSpot_2),2:3),rgb_2(find(IsRearSpot_2),2:3),1,2);

     figure(24);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_1 = ClusterPointsMulti(rgb_1(:,2:3),rgb_2(find(IsFrontSpot_2),2:3),rgb_2(find(IsRearSpot_2),2:3),1,2);
     
     case 4
     
     figure(23);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_1 = ClusterPointsMulti(PCArgb_1(:,1:2),PCArgb_2(find(IsFrontSpot_2),1:2),PCArgb_2(find(IsRearSpot_2),1:2),1,2);


     figure(24);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_1 = ClusterPointsMulti(PCArgb_1(:,1:2),PCArgb_2(find(IsFrontSpot_2),1:2),PCArgb_2(find(IsRearSpot_2),1:2),1,2);
     
     case 5
     
     figure(23);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_1 = ClusterPointsMulti(PCArgb_1(:,[1 3]),PCArgb_2(find(IsFrontSpot_2),[1 3]),PCArgb_2(find(IsRearSpot_2),[1 3]),1,2);

     figure(24);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_1 = ClusterPointsMulti(PCArgb_1(:,[1 3]),PCArgb_2(find(IsFrontSpot_2),[1 3]),PCArgb_2(find(IsRearSpot_2),[1 3]),1,2);
     
     case 6
     
     figure(23);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_1 = ClusterPointsMulti(PCArgb_1(:,2:3),PCArgb_2(find(IsFrontSpot_2),2:3),PCArgb_2(find(IsRearSpot_2),2:3),1,2);

     figure(24);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_1 = ClusterPointsMulti(PCArgb_1(:,2:3),PCArgb_2(find(IsFrontSpot_2),2:3),PCArgb_2(find(IsRearSpot_2),2:3),1,2);

     
  end
   
  figure(25);
  
  subplot(2,3,1);
  plot(rgb_1(:,1), rgb_1(:,2),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(rgb_1(find(IsRearSpot_1),1),rgb_1(find(IsRearSpot_1),2),'x','color',[0 0 1],'markersize',4,'linestyle','none');
	plot(rgb_1(find(IsFrontSpot_1),1),rgb_1(find(IsFrontSpot_1),2),'.','color',[1 0 0],'markersize',5,'linestyle','none');hold off;
  title('1');
  xlabel('Red');
  ylabel('Green');hold off;
  
  subplot(2,3,2);
  plot(rgb_1(:,1), rgb_1(:,3),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(rgb_1(find(IsRearSpot_1),1),rgb_1(find(IsRearSpot_1),3),'x','color',[0 0 1],'markersize',4,'linestyle','none');
	plot(rgb_1(find(IsFrontSpot_1),1),rgb_1(find(IsFrontSpot_1),3),'.','color',[1 0 0],'markersize',5,'linestyle','none');hold off;
  title('2');
  xlabel('Red');
  ylabel('Blue');hold off;
  
  subplot(2,3,3);
  plot(rgb_1(:,2), rgb_1(:,3),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(rgb_1(find(IsRearSpot_1),2),rgb_1(find(IsRearSpot_1),3),'x','color',[0 0 1],'markersize',4,'linestyle','none');
	plot(rgb_1(find(IsFrontSpot_1),2),rgb_1(find(IsFrontSpot_1),3),'.','color',[1 0 0],'markersize',5,'linestyle','none');hold off;
  title('3');
  xlabel('Green');
  ylabel('Blue');hold off;
  
  subplot(2,3,4);
  plot(PCArgb_1(:,1), PCArgb_1(:,2),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(PCArgb_1(find(IsRearSpot_1),1),PCArgb_1(find(IsRearSpot_1),2),'x','color',[0 0 1],'markersize',4,'linestyle','none');
	plot(PCArgb_1(find(IsFrontSpot_1),1),PCArgb_1(find(IsFrontSpot_1),2),'.','color',[1 0 0],'markersize',5,'linestyle','none');hold off;
  title('4');
  xlabel('PC1');
  ylabel('PC2');hold off;
  
  subplot(2,3,5);
  plot(PCArgb_1(:,1), PCArgb_1(:,3),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
   plot(PCArgb_1(find(IsRearSpot_1),1),PCArgb_1(find(IsRearSpot_1),3),'x','color',[0 0 1],'markersize',4,'linestyle','none');
	plot(PCArgb_1(find(IsFrontSpot_1),1),PCArgb_1(find(IsFrontSpot_1),3),'.','color',[1 0 0],'markersize',5,'linestyle','none');hold off;
  title('5');
  xlabel('PC1');
  ylabel('PC3');hold off;
  
  subplot(2,3,6);
  plot(PCArgb_1(:,2), PCArgb_1(:,3),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(PCArgb_1(find(IsRearSpot_1),2),PCArgb_1(find(IsRearSpot_1),3),'x','color',[0 0 1],'markersize',4,'linestyle','none');
	plot(PCArgb_1(find(IsFrontSpot_1),2),PCArgb_1(find(IsFrontSpot_1),3),'.','color',[1 0 0],'markersize',5,'linestyle','none');hold off;
  title('6');
  xlabel('PC2');
  ylabel('PC3');hold off;
  
  
  %% Get position information using the index of spot files.
  %% Sometimes the color assinment is opposite.
  OKcolor = input('\n Is the color assinment is OK (or opposite)? y/n','s');
  while 1;
    if strcmp(OKcolor,'y'),
      FrontSpots_1 = Spots_1(find(IsFrontSpot_1));
      RearSpots_1 = Spots_1(find(IsRearSpot_1));
		WeirdSpots_1 = Spots_1(find(~IsFrontSpot_1 & ~IsRearSpot_1));
      break
    end
    if strcmp(OKcolor,'n'),
      FrontSpots_1 = Spots_1(find(IsRearSpot_1));
      RearSpots_1 = Spots_1(find(IsFrontSpot_1));
		WeirdSpots_1 = Spots_1(find(~IsFrontSpot_1 & ~IsRearSpot_1));
      break
    end
  end
  
if length(FrontSpots_1)~=0;
	FrontFrames_1 = Accumulate(1+spots(FrontSpots_1,1),1);
	GoodFrontSpots_1 = FrontSpots_1(find(FrontFrames_1(1+spots(FrontSpots_1,1))==1));
	GoodFrontSpotsRatio_1 = length(GoodFrontSpots_1)./length(FrontSpots_1)
	FrontSpots_1 = GoodFrontSpots_1;
end
if length(RearSpots_1)~=0;
	RearFrames_1 = Accumulate(1+spots(RearSpots_1,1),1);
	GoodRearSpots_1 = RearSpots_1(find(RearFrames_1(1+spots(RearSpots_1,1))==1));
	GoodRearSpotsRatio_1 = length(GoodRearSpots_1)./length(RearSpots_1)
	RearSpots_1 = GoodRearSpots_1;
end

	WiredSpotsRatio_1 = length(WeirdSpots_1)./length(Spots_1)

	figure(26);

    subplot(2,3,1);
    plot(spots(FrontSpots_1,2),spots(FrontSpots_1,7),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(FrontSpots_1,2), number of pixels');
    ylabel('spots(FrontSpots_1,7), luminance');
    title('FrontSpots_2');

	subplot(2,3,2);
    plot(spots(RearSpots_1,2),spots(RearSpots_1,7),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(RearSpots_1,2), number of pixels');
    ylabel('spots(RearSpots_1,7), luminance');
    title('RearSpots_1');

	subplot(2,3,3);
    plot(spots(WeirdSpots_1,2),spots(WeirdSpots_1,7),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(WeirdSpots_1,2), number of pixels');
    ylabel('spots(WeirdSpots_1,7), luminance');
    title('WeirdSpots_1');

	subplot(2,3,4);
    plot(spots(FrontSpots_1,3),spots(FrontSpots_1,4),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(FrontSpots_1,3), x');
    ylabel('spots(FrontSpots_1,4), y');
    title('FrontSpots_1');

	subplot(2,3,5);
    plot(spots(RearSpots_1,3),spots(RearSpots_1,4),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(RearSpots_1,3), x');
    ylabel('spots(RearSpots_1,4), y');
    title('WeirdSpots_1');

	subplot(2,3,6);
    plot(spots(WeirdSpots_1,3),spots(WeirdSpots_1,4),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(WeirdSpots_1,3), x');
    ylabel('spots(WeirdSpots_1,4), y');
    title('WeirdSpots_1');

  figure(27);
  plot(spots(FrontSpots_1,3),spots(FrontSpots_1,4),'.','color',[1 0 0],'markersize',10,'linestyle','none');hold on;
  plot(spots(RearSpots_1,3),spots(RearSpots_1,4),'x','color',[0 1 0],'markersize',10,'linestyle','none');hold off;
	title('One spot in a frame, whl file');

%%%%%%%%%%%%%%%%
FrontSpots_3 = [];
RearSpots_3 = [];
if ~isempty(Spots_3);
	fprintf('\n Spots_3. Draw boarders of a cluster by clicking into a figure. Do not click on the edge of the figure.');
  switch WhichP
     case 1

     figure(28);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_3 = ClusterPointsMulti(rgb_3(:,1:2),rgb_2(find(IsFrontSpot_2),1:2),rgb_2(find(IsRearSpot_2),1:2),1,2);

     figure(29);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_3 = ClusterPointsMulti(rgb_3(:,1:2),rgb_2(find(IsFrontSpot_2),1:2),rgb_2(find(IsRearSpot_2),1:2),1,2);
     case 2
     
     figure(28);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_3 = ClusterPointsMulti(rgb_3(:,[1 3]),rgb_2(find(IsFrontSpot_2),[1 3]),rgb_2(find(IsRearSpot_2),[1 3]),1,2);

     figure(29);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_3 = ClusterPointsMulti(rgb_3(:,[1 3]),rgb_2(find(IsFrontSpot_2),[1 3]),rgb_2(find(IsRearSpot_2),[1 3]),1,2);
     
     case 3
     
     figure(28);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_3 = ClusterPointsMulti(rgb_3(:,2:3),rgb_2(find(IsFrontSpot_2),2:3),rgb_2(find(IsRearSpot_2),2:3),1,2);

     figure(29);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_3 = ClusterPointsMulti(rgb_3(:,2:3),rgb_2(find(IsFrontSpot_2),2:3),rgb_2(find(IsRearSpot_2),2:3),1,2);
     
     case 4
     
     figure(28);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_3 = ClusterPointsMulti(PCArgb_3(:,1:2),PCArgb_2(find(IsFrontSpot_2),1:2),PCArgb_2(find(IsRearSpot_2),1:2),1,2);


     figure(29);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_3 = ClusterPointsMulti(PCArgb_3(:,1:2),PCArgb_2(find(IsFrontSpot_2),1:2),PCArgb_2(find(IsRearSpot_2),1:2),1,2);
     
     case 5
     
     figure(28);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_3 = ClusterPointsMulti(PCArgb_3(:,[1 3]),PCArgb_2(find(IsFrontSpot_2),[1 3]),PCArgb_2(find(IsRearSpot_2),[1 3]),1,2);

     figure(29);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_3 = ClusterPointsMulti(PCArgb_3(:,[1 3]),PCArgb_2(find(IsFrontSpot_2),[1 3]),PCArgb_2(find(IsRearSpot_2),[1 3]),1,2);
     
     case 6
     
     figure(28);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsFrontSpot_3 = ClusterPointsMulti(PCArgb_3(:,2:3),PCArgb_2(find(IsFrontSpot_2),2:3),PCArgb_2(find(IsRearSpot_2),2:3),1,2);

     figure(29);
     
     title('One-Spot in a frame. Draw boarders of a cluster by clicking into a figure.');
     IsRearSpot_3 = ClusterPointsMulti(PCArgb_3(:,2:3),PCArgb_2(find(IsFrontSpot_2),2:3),PCArgb_2(find(IsRearSpot_2),2:3),1,2);

     
  end
   
  figure(30);
  
  subplot(2,3,1);
  plot(rgb_3(:,1), rgb_3(:,2),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(rgb_3(find(IsFrontSpot_3),1),rgb_3(find(IsFrontSpot_3),2),'.','color',[1 0 0],'markersize',3,'linestyle','none');
  plot(rgb_3(find(IsRearSpot_3),1),rgb_3(find(IsRearSpot_3),2),'.','color',[0 0 1],'markersize',3,'linestyle','none');hold off;
  title('1');
  xlabel('Red');
  ylabel('Green');hold off;
  
  subplot(2,3,2);
  plot(rgb_3(:,1), rgb_3(:,3),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(rgb_3(find(IsFrontSpot_3),1),rgb_3(find(IsFrontSpot_3),3),'.','color',[1 0 0],'markersize',3,'linestyle','none');
  plot(rgb_3(find(IsRearSpot_3),1),rgb_3(find(IsRearSpot_3),3),'.','color',[0 0 1],'markersize',3,'linestyle','none');hold off;
  title('2');
  xlabel('Red');
  ylabel('Blue');hold off;
  
  subplot(2,3,3);
  plot(rgb_3(:,2), rgb_3(:,3),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(rgb_3(find(IsFrontSpot_3),2),rgb_3(find(IsFrontSpot_3),3),'.','color',[1 0 0],'markersize',3,'linestyle','none');
  plot(rgb_3(find(IsRearSpot_3),2),rgb_3(find(IsRearSpot_3),3),'.','color',[0 0 1],'markersize',3,'linestyle','none');hold off;
  title('3');
  xlabel('Green');
  ylabel('Blue');hold off;
  
  subplot(2,3,4);
  plot(PCArgb_3(:,1), PCArgb_3(:,2),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(PCArgb_3(find(IsFrontSpot_3),1),PCArgb_3(find(IsFrontSpot_3),2),'.','color',[1 0 0],'markersize',3,'linestyle','none');
  plot(PCArgb_3(find(IsRearSpot_3),1),PCArgb_3(find(IsRearSpot_3),2),'.','color',[0 0 1],'markersize',3,'linestyle','none');hold off;
  title('4');
  xlabel('PC1');
  ylabel('PC2');hold off;
  
  subplot(2,3,5);
  plot(PCArgb_3(:,1), PCArgb_3(:,3),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(PCArgb_3(find(IsFrontSpot_3),1),PCArgb_3(find(IsFrontSpot_3),3),'.','color',[1 0 0],'markersize',3,'linestyle','none');
  plot(PCArgb_3(find(IsRearSpot_3),1),PCArgb_3(find(IsRearSpot_3),3),'.','color',[0 0 1],'markersize',3,'linestyle','none');hold off;
  title('5');
  xlabel('PC1');
  ylabel('PC3');hold off;
  
  subplot(2,3,6);
  plot(PCArgb_3(:,2), PCArgb_3(:,3),'.','color',[.1 .1 .1],'markersize',3,'linestyle','none');hold on
  plot(PCArgb_3(find(IsFrontSpot_3),2),PCArgb_3(find(IsFrontSpot_3),3),'.','color',[1 0 0],'markersize',3,'linestyle','none');
  plot(PCArgb_3(find(IsRearSpot_3),2),PCArgb_3(find(IsRearSpot_3),3),'.','color',[0 0 1],'markersize',3,'linestyle','none');hold off;
  title('6');
  xlabel('PC2');
  ylabel('PC3');hold off;
  
  
  %% Get position information using the index of spot files.
  %% Sometimes the color assinment is opposite.
  OKcolor = input('\n Is the color assinment is OK (or opposite)? y/n','s');
  while 1;
    if strcmp(OKcolor,'y'),
      FrontSpots_3 = Spots_3(find(IsFrontSpot_3 & ~IsRearSpot_3));
      RearSpots_3 = Spots_3(find(IsRearSpot_3 & ~IsFrontSpot_3));
		WeirdSpots_3 = Spots_3(find(~IsFrontSpot_3 & ~IsRearSpot_3));
      break
    end
    if strcmp(OKcolor,'n'),
      FrontSpots_3 = Spots_3(find(IsRearSpot_3 & ~IsFrontSpot_3));
      RearSpots_3 = Spots_3(find(IsFrontSpot_3 & ~IsRearSpot_3));
		WeirdSpots_3 = Spots_3(find(~IsFrontSpo_3 & ~IsRearSpot_3));
      break
    end
  end

	if ~isempty(FrontSpots_3);
		FrontFrames_3 = Accumulate(1+spots(FrontSpots_3,1),1);
		GoodFrontSpots_3 = FrontSpots_3(find(FrontFrames_3(1+spots(FrontSpots_3,1))==1));
%		MultiFrondSpots_3 = FrontSpots_3(find(FrontFrames_3(1+spots(FrontSpots_3,1))>1));
%%%% you can use the aligorism by which the distance between Front and Rear LEd is used to identify real front or rear led.
		GoodFrontSpotsRatio_3 = length(GoodFrontSpots_3)./length(FrontSpots_3)
		FrontSpots_3 = GoodFrontSpots_3;
	end

	if ~isempty(RearSpots_3);
		RearFrames_3 = Accumulate(1+spots(RearSpots_3,1),1);
		GoodRearSpots_3 = RearSpots_3(find(RearFrames_3(1+spots(RearSpots_3,1))==1));
%		MultiRearSpots_3 = RearSpots_3(find(RearFrames_3(1+spots(RearSpots_3,1))>1));
		GoodRearSpotsRatio_3 = length(GoodRearSpots_3)./length(RearSpots_3)
		RearSpots_3 = GoodRearSpots_3;
	end



	WiredSpotsRatio_3 = length(WeirdSpots_3)./length(Spots_3)

	figure(31);

    subplot(2,3,1);
    plot(spots(FrontSpots_3,2),spots(FrontSpots_3,7),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(FrontSpots_3,2), number of pixels');
    ylabel('spots(FrontSpots_3,7), luminance');
    title('FrontSpots_3');

	subplot(2,3,2);
    plot(spots(RearSpots_3,2),spots(RearSpots_3,7),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(RearSpots_3,2), number of pixels');
    ylabel('spots(RearSpots_3,7), luminance');
    title('RearSpots_3');

	subplot(2,3,3);
    plot(spots(WeirdSpots_3,2),spots(WeirdSpots_3,7),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(RearSpots_3,2), number of pixels');
    ylabel('spots(RearSpots_3,7), luminance');
    title('WeirdSpots_3');

	subplot(2,3,4);
    plot(spots(FrontSpots_3,3),spots(FrontSpots_3,4),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(FrontSpots_3,3), x');
    ylabel('spots(FrontSpots_3,4), y');
    title('FrontSpots_3');

	subplot(2,3,5);
    plot(spots(RearSpots_3,3),spots(RearSpots_3,4),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(RearSpots_3,3), x');
    ylabel('spots(RearSpots_3,4), y');
    title('WeirdSpots_3');

	subplot(2,3,6);
    plot(spots(WeirdSpots_3,3),spots(WeirdSpots_3,4),'.','markersize',4,'linestyle','none','color',[0 0 0]);
    xlabel('spots(WeirdSpots_3,3), x');
    ylabel('spots(WeirdSpots_3,4), y');
    title('WeirdSpots_3');

  figure(32);
  plot(spots(FrontSpots_3,3),spots(FrontSpots_3,4),'.','color',[1 0 0],'markersize',10,'linestyle','none');hold on;
  plot(spots(RearSpots_1,3),spots(RearSpots_1,4),'.','color',[0 1 0],'markersize',10,'linestyle','none');hold off;
	title('One spot in a frame, whl file');

end


%%%%%%%%%%%%%%%

  %% now make trajectory
HeadPos = -1*ones(nFrames,4);
HeadPos(1+spots(FrontSpots_1,1),1:2) = spots(FrontSpots_1,3:4);
HeadPos(1+spots(FrontSpots_2,1),1:2) = spots(FrontSpots_2,3:4);
if ~isempty(FrontSpots_3);
	HeadPos(1+spots(FrontSpots_3,1),1:2) = spots(FrontSpots_3,3:4);
end


HeadPos(1+spots(RearSpots_1,1),3:4) = spots(RearSpots_1,3:4);
HeadPos(1+spots(RearSpots_2,1),3:4) = spots(RearSpots_2,3:4);
if ~isempty(RearSpots_3);
	HeadPos(1+spots(RearSpots_3,1),3:4) = spots(RearSpots_3,3:4);
end


BadFrameInd = find(HeadPos(:,1)==-1 | HeadPos(:,2)==-1 | HeadPos(:,3)==-1 | HeadPos(:,4)==-1);
BadFrontFrameInd = find(HeadPos(:,1)==-1 & HeadPos(:,2)==-1);
BadOnlyFrontFrameInd = find(HeadPos(:,1)==-1 & HeadPos(:,2)==-1 & HeadPos(:,3)~=-1 & HeadPos(:,4)~=-1);
BadRearFrameInd = find(HeadPos(:,3)==-1 & HeadPos(:,4)==-1);
BadOnlyRearFrameInd = find(HeadPos(:,1)~=-1 & HeadPos(:,2)~=-1 & HeadPos(:,3)==-1 & HeadPos(:,4)==-1);

BadFrameRatio = length(BadFrameInd)./length(HeadPos)
BadFrontFrameRatio = length(BadFrontFrameInd)./length(HeadPos)
BadOnlyFrontFrameRatio = length(BadOnlyFrontFrameInd)./length(HeadPos)
BadRearFrameRatio = length(BadRearFrameInd)./length(HeadPos)
BadOnlyRearFrameRatio = length(BadOnlyRearFrameInd)./length(HeadPos)

%%%%%%%%%%%%%%% Conpensate the XY ratio, and the deformity caused by camera lens. Compensate before you use le
if VideoLensBool==0;
	[HeadPos] =  WhlCompensateXYMmaze_060407(HeadPos,fileBase,0);
end
if VideoLensBool==1;
	[HeadPos] = WhlCompensateXYBig_060407(HeadPos,fileBase,0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
FR_Dist = sqrt((HeadPos(:,1)-HeadPos(:,3)).^2 + (HeadPos(:,2)-HeadPos(:,4)).^2);
%%%%%%%%%%%%%% discard the data if the distance between front and rear Led is more than 16cm.
badDist = find(FR_Dist > ledDist & HeadPos(:,1) ~= -1 & HeadPos(:,3) ~= -1);
HeadPos(badDist,:) = -1;

  %% interpolate missing stretches up to 30 frames long
  cHeadPos = CleanWhlForR(HeadPos, 90, 60);
  cHeadPos((find(cHeadPos==-1)))=NaN; % so it doesn't interpolate between 100 and -1 and get 50.
  %% interp1 dose not use NaN for interporate.

  %% now make wheel file by interpolating
%%% If you use too high frequency for sync pulse, it is problemaitc to use interp1( , , ,'linear'), because the fluctuation of the sapmling vias (1/30 Hz = 30ms) is always refrected to the interporated resultant.

  TargetSamples = 0:32:length(eeg);
  GoodFrames = find(isfinite(FrameSamples));
  Whl(:,1:2) = interp1(FrameSamples(GoodFrames),cHeadPos(GoodFrames,1:2),TargetSamples,'linear',-1);
	Whl(:,3:4) = interp1(FrameSamples(GoodFrames),cHeadPos(GoodFrames,3:4),TargetSamples,'linear',-1);
  Whl(find(~isfinite(Whl)))=-1;
  %%% The length of whl file is defined by the length of eeg file.
	%%% In this way the samlping frequency of the video is modified from 29.97 Hz to 1250/32 = 39.0625 Hz.
  
  %%%
  whlLed = [];
  if NumWhlLed == 1 | NumWhlLed ==2;
  whlLed(:,:) = round(interp1(FrameSamples(GoodFrames),wLed(GoodFrames,:),TargetSamples,'linear',-1));
  end
  
  figure(33);hold on;
     Xmax = max(max(Whl(:,1)),max(Whl(:,3)));
     Ymax = max(max(Whl(:,2)),max(Whl(:,4)));
     set(gca,'xlim',[-20 Xmax+20],'ylim',[-20 Ymax+20]);hold on;
     plot(Whl(:,1),Whl(:,2),'.','color',[1 0 0],'markersize',15,'linestyle','none');
     plot(Whl(:,3),Whl(:,4),'.','color',[0 1 0],'markersize',10,'linestyle','none');
     title('Compensated,smoothed whl file');hold off;


  figure(34);hold on;
  DetectionRate_F = [];
  DetectionRateAve_F = length(find(Whl(:,1)~=-1))./ (max(find(Whl(:,1)~=-1)) - min(find(Whl(:,1)~=-1)) +1),
	DetectionRateAve_R = length(find(Whl(:,3)~=-1))./ (max(find(Whl(:,3)~=-1)) - min(find(Whl(:,3)~=-1)) +1),
  k=0;
  fprintf('\nThe Whl file data\n-----------------------------\n');
  nFigFrames = ceil(length(Whl)/100);
  %% The last <100 samples will not be shown.
	%Xmax = max(max(Whl(:,1)),max(Whl(:,3)));
	%Ymax = max(max(Whl(:,2)),max(Whl(:,4)));
  while ~strcmp(input('   Hit ENTER to show the next 100 samples, or type ''done''+ENTER to proceed...','s'),'done'),
     k = k+1;
     if k*100 > length(Whl), break; end
     subplot(2,1,1);hold on;
     cla;
     title('cleaned and interporated .whl file');
     set(gca,'xlim',[-20 Xmax+20],'ylim',[-20 Ymax+20]);
     plot(Whl((k-1)*100+1:k*100,1),Whl((k-1)*100+1:k*100,2),'.','color',[1 0 0],'markersize',15,'linestyle','none');
     plot(Whl((k-1)*100+1:k*100,3),Whl((k-1)*100+1:k*100,4),'.','color',[0 1 0],'markersize',10,'linestyle','none');
     for j=(k-1)*100+1:k*100, line([Whl(j,1) Whl(j,3)],[Whl(j,2) Whl(j,4)],'color',[0 0 0]); end
     subplot(2,1,2);hold on;
     title(['Detection Rate(average and instantaneous)' ' ' num2str(k) '/' num2str(nFigFrames)]);
     DetectionRate_F = [DetectionRate_F length(find(Whl((k-1)*100+1:k*100,1)~=-1))/100];
     plot(DetectionRate_F,'color',[1 0 0]);
     plot([k-1 k],[DetectionRateAve_F DetectionRateAve_F]);
  end
end

%%% k=0;
%%% plotWhl = Whl(find(Whl(:,1)~=-1),:);
%%% fprintf('\nThe Whl file data\n----------------------\n');
%%%% end


if ~isempty(events),
  % Update events: remove all initial SYNC events and replace them with the corrected ones (according to user input)
  nonSYNC = events(events(:,2) ~= 83 | events(:,3) ~= 89,:);
  events = [nonSYNC;[syncEvents 83*ones(size(syncEvents)) 89*ones(size(syncEvents))]];
  events = sortrows(events,1);
  % Throw any events occurring after the last SYNC event
  lastSYNC = find(events(:,2) == 83 & events(:,3) == 89);lastSYNC = lastSYNC(end);
  events = events(1:lastSYNC,:);
  % Synchronize events on electrophysiological data
  timestamps = interp1(syncEvents,syncEEG(1:length(syncEvents))/1250*1000,events(:,1),'linear',-1);
  events(:,1) = timestamps;
  events = double(uint32(events));
end

%%% if ~isempty(spots),
%%% figure(1);
%%% end
while 1,
  i = input('\n\nSave to disk (yes/no)? ', 's');
  if strcmp(i,'yes') | strcmp(i,'no'), break; end
end
if i(1) == 'y'
  if ~isempty(spots),
     fprintf('Saving %s\n', [fileBase '.whl']);
     msave([fileBase '.whl'], Whl);
  end
  if ~isempty(whlLed),
     fprintf('Saving %s\n', [fileBase '.whlLed']);
     msave([fileBase '.whlLed'], whlLed');
  end
  if ~isempty(events),
     fprintf('Saving %s\n', [fileBase '.evt']);
     msave([fileBase '.evt'],events);
  end
end


%%% if VideoLensBool==0;
%	WhlCompensateXYMaze_060407(Whl,fileBase)
%end
%if VideoLensBool==1;
%	WhlCompensateXYBig_060407(Whl,fileBase)
%end

	

if 0 % no longer needed because neuroscope is cool
   while 1,
   i = input('\n\nSave all position information as jpeg (yes/no)? ', 's');
   if strcmp(i,'yes') | strcmp(i,'no'), break; end
   end
   if i(1)=='y'
     figure;
     plot(plotWhl(:,1),plotWhl(:,2), '.','color',[0.5 0.5 0.5],'markersize',10,'linestyle','none' );
     set(gca, 'xlim', [0 videoRes(1)], 'ylim', [0 videoRes(2)], 'Position', [0 0 1 1]);
     set(gca, 'color', [0 0 0]);
     set(gcf, 'inverthardcopy', 'off')
     print(gcf, '-djpeg100', [fileBase '.jpg']);
     eval(['!convert -geometry ' num2str(videoRes(1)) 'x' num2str(videoRes(2)) ' ' fileBase '.jpg' ' ' fileBase '.jpg']);
   end
end
end







