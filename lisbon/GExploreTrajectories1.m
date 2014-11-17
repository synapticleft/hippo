function [xValsCorr,yValsCorr,xValsErr,yValsErr]=GExploreTrajectories1(filename,ses,condition,Nsaturate)

% This function plots trajectories of Roberto's human task. It loads a data
% file for a subject and plots the trajectories from sound onset to target
% arrival both as a line for each trial and as a 2D histogram. Optional
% arguments are ses (specifies sessions to be used as a 1D array with
% elements between 1 and 12, defaults to zero, which selects all sessions),
% condition (specifies condition as an integer between 1 and 6, defaults to
% zero, which plots all conditions) and Nsaturate (positive integer for
% color saturation of the 2D histogram, defaults to 50)
%
% Example ussage
%
% ExploreTrajectories('data_s4_interp_2sPrePost_wRaw'); plots all
% conditions in all sessions.
%
% ExploreTrajectories('data_s4_interp_2sPrePost_wRaw',[4:12],4); plots
% trajectories from condition 4 in sessions 4 to 12.

file = dir('*.mat');
load(file(filename).name,'data');
file(filename).name

%load(filename)

% This loads the cell array with all the data interpolated at 80 Hz (2 sec
% of baseline plus the stimulus plus 2 sec poststim), a scalar center_ILD,
% with the ILD center for the subject (in dB SPL) and the array
% bin_centers, which contains the ILDs that can be used in the x-axis of
% the psychometric curve.


session=data(:,2);  % session
trialN=data(:,3);  % trial number 
ILD_th=data(:,4);  % long term average of the ILD in this trial
ILD_pr=data(:,5);  % mean ILD actually experienced (not including initial constant ILD)
stim=data(:,6);    % time varying stimulus for each trial relative to subject's ILD=0
outcome=data(:,7); % 1=left; 2=right; 3=abort 
cond=data(:,8);    % 1,2,3,4,5,6 and NaN (aborts and super easy)
RT=data(:,9);      % reaction time
pupX=data(:,10);   % x coord of pupil
pupY=data(:,11);   % y coord of pupil
pupMaj=data(:,12); % pupil major axis
pupMin=data(:,13); % pupil minor axis
pupAng=data(:,14); % pupil angle (of the ellipse)
jX=data(:,15);     % joystick x coord
jY=data(:,16);     % joystick y coord
musc1=data(:,17);  % emg1
musc2=data(:,18);  % emg2
accX=data(:,19);   % accelerometer x
accY=data(:,20);   % accelerometer y
accZ=data(:,21);   % accelerometer z
time=data(:,22);   % time within trial
pupXr=data(:,23);   % x coord of pupil raw
pupYr=data(:,24);   % y coord of pupil raw
pupMajr=data(:,25); % pupil major axis raw
pupMinr=data(:,26); % pupil minor axis raw
pupAngr=data(:,27); % pupil angle (of the ellipse) raw
pupT=data(:,28); % pupil raw time

cols = jet(64);


% optional arguments
if nargin<2
    ses=0; %plot all sessions
end

if nargin<3
    condition=0; %plot all conditions
end

if nargin<4
    Nsaturate=50; %color saturation for the 2D histograms
end


if condition>0 % Plot trajectories for a specific condition only
    
    trV=[];
    if sum(ses~=0)>0  % Specify both session and trial type
        disp(' ')
        disp(['Selecting trials from condition = ' num2str(condition) ' from session = ' num2str(ses)])
        disp(' ')        
        tit=['Trajectories from condition = ' num2str(condition) ' from session = ' num2str(ses)];
        sesMat=cell2mat(session(2:end));
        condMat=cell2mat(cond(2:end));
        for ises=1:length(ses)
            curr_ses=ses(ises);
            trV=[trV;find(sesMat==curr_ses & condMat==condition)+1];
        end
        if isempty(trV)
            disp(' ')
            disp('WATCH OUT!! SESSION + TRIAL SPECIFICATION FAILED!!!')
            disp(' ')
            return
        end        
    elseif sum(ses~=0)==0  % Specify trial type and no particular session
        disp(' ')
        disp(['Selecting trials from condition = ' num2str(condition) ' from all sessions.'])
        disp(' ')
        tit=['Trajectories from condition = ' num2str(condition) ' from all sessions'];
        condMat=cell2mat(cond(2:end));
        trV=find(condMat==condition)+1;
        if isempty(trV)
            disp(' ')
            disp('WATCH OUT!! TRIAL SPECIFICATION FAILED!!!')
            disp(' ')
            return
        end         
    end
    
    if isempty(trV) %check that at least some trials were specified
        disp(' ')
        disp('WATCH OUT!! TRIAL SPECIFICATION PROCESS FAILED!!!')
        disp(' ')
        return
    end  

    % Figure Params
    Nrows=2;
    Ncols=2;
    delR=0.02;
    delC=0.02;
    pX=(0.9-(Ncols-1)*delC)/Ncols;
    pY=(0.9-(Nrows-1)*delR)/Nrows;
    
    figure('position',[140         294        1166        1027],'name',tit);
    
    a(1)=axes('position',[0.05 0.95-pY pX pY]);
    a(2)=axes('position',[0.05+pX+delC 0.95-pY pX pY]);
    a(3)=axes('position',[0.05 0.95-2*pY-delR pX pY]);
    a(4)=axes('position',[0.05+pX+delC 0.95-2*pY-delR pX pY]);
    
    xValsCorr=[];
    yValsCorr=[];
    xValsErr=[];
    yValsErr=[];
    for i=1:length(trV)

        tr=trV(i);
        istim=~isnan(stim{tr});

        if condition>3 & outcome{tr}==2 %correct trial
            subplot(a(1));
            plot(jX{tr}(istim),jY{tr}(istim),'b'),hold on
            xValsCorr=[xValsCorr jX{tr}(istim)];
            yValsCorr=[yValsCorr jY{tr}(istim)];
        elseif condition<4 & outcome{tr}==1 %correct trial
            subplot(a(1));
            plot(jX{tr}(istim),jY{tr}(istim),'b'),hold on
            xValsCorr=[xValsCorr jX{tr}(istim)];
            yValsCorr=[yValsCorr jY{tr}(istim)];
        elseif condition<4 & outcome{tr}==2 %error trial
            subplot(a(2))
            plot(jX{tr}(istim),jY{tr}(istim),'r'),hold on
            xValsErr=[xValsErr jX{tr}(istim)];
            yValsErr=[yValsErr jY{tr}(istim)];
        elseif condition>3 & outcome{tr}==1 %error trial
            subplot(a(2))
            plot(jX{tr}(istim),jY{tr}(istim),'r'),hold on
            xValsErr=[xValsErr jX{tr}(istim)];
            yValsErr=[yValsErr jY{tr}(istim)];
        end

    end
    
    subplot(a(1))
    %xlabel('Joystick X coordinate (AU)')
    ylabel('Joystick Y coordinate (AU)')
    xlim([-72 72])
    ylim([-5 150])
    set(gca,'yticklabel',[],'xticklabel',[])
    %legend('Correct',3)
    title('Correct Trials')
    subplot(a(2))
    %xlabel('Joystick X coordinate (AU)')
    set(gca,'yticklabel',[],'xticklabel',[])
    xlim([-70 70])
    ylim([-5 150])
    title('Error Trials')
    %legend('Errors',3)
    
    ctrs{1}=-70:4:70;ctrs{2}=-2.5:5:147.5;

    nCorr=hist3([xValsCorr' yValsCorr'],ctrs);
    subplot(a(3))
    imagesc(ctrs{1},ctrs{2},nCorr(1:end,end:-1:1)');
    set(gca,'yticklabel',[],'xticklabel',[])
    colormap('hot')
    caxis([0 Nsaturate])
    xlabel('Joystick X coordinate (AU)')
    ylabel('Joystick Y coordinate (AU)')

    if ~isempty(xValsErr) %do only if there are error trials
        nErr=hist3([xValsErr' yValsErr'],ctrs);
    else
        nErr=zeros(length(ctrs{1}),length(ctrs{2}));
    end
    subplot(a(4))
    imagesc(ctrs{1},ctrs{2},nErr(1:end,end:-1:1)');
    set(gca,'yticklabel',[],'xticklabel',[])
    caxis([0 Nsaturate])
    xlabel('Joystick X coordinate (AU)')
    ylabel('Joystick Y coordinate (AU)')
        
else %if condition==0, plot all conditions sequentially

    figure('position',[33         117        2141        1204]);
    
    % Figure Params
    Nrows=4;
    Ncols=6;
    delR=0.02;
    delC=0.02;
    pX=(0.9-(Ncols-1)*delC)/Ncols;
    pY=(0.9-(Nrows-1)*delR)/Nrows;
    
    for ic=1:6
    condition=ic;
    
        trV=[];
        if sum(ses~=0)>0  % Specify both session and trial type
            disp(' ')
            disp(['Selecting trials from condition = ' num2str(condition) ' from session = ' num2str(ses)])
            disp(' ')        
            sesMat=cell2mat(session(2:end));
            condMat=cell2mat(cond(2:end));
            for ises=1:length(ses)
                curr_ses=ses(ises);
                trV=[trV;find(sesMat==curr_ses & condMat==condition)+1];
            end
            if isempty(trV)
                disp(' ')
                disp('WATCH OUT!! SESSION + TRIAL SPECIFICATION FAILED!!!')
                disp(' ')
                return
            end        
        elseif sum(ses~=0)==0  % Specify trial type and no particular session
            disp(' ')
            disp(['Selecting trials from condition = ' num2str(condition) ' from all sessions.'])
            disp(' ')
            condMat=cell2mat(cond(2:end));
            trV=find(condMat==condition)+1;
            if isempty(trV)
                disp(' ')
                disp('WATCH OUT!! TRIAL SPECIFICATION FAILED!!!')
                disp(' ')
                return
            end         
        end

        if isempty(trV) %check that at least some trials were specified
            disp(' ')
            disp('WATCH OUT!! TRIAL SPECIFICATION PROCESS FAILED!!!')
            disp(' ')
            return
        end 
        
        a((ic-1)*4+1)=axes('position',[0.05+(pX+delC)*(ic-1) 0.95-pY pX pY]);
        a((ic-1)*4+2)=axes('position',[0.05+(pX+delC)*(ic-1) 0.95-2*pY-delR pX pY]);
        a((ic-1)*4+3)=axes('position',[0.05+(pX+delC)*(ic-1) 0.95-3*pY-2*delR pX pY]);
        a((ic-1)*4+4)=axes('position',[0.05+(pX+delC)*(ic-1) 0.95-4*pY-3*delR pX pY]);
        
        xValsCorr=[];
        yValsCorr=[];
        xValsErr=[];
        yValsErr=[];
        for i=1:length(trV)

            tr=trV(i);
            istim=~isnan(stim{tr});

            if (condition>3 && outcome{tr}==2) || (condition<4 && outcome{tr}==1) %correct trial
                subplot(a((ic-1)*4+1));
                %plot(jX{tr}(istim),jY{tr}(istim),'b'),hold on
                scatter(jX{tr}(istim),jY{tr}(istim),3,cols(min(1:sum(istim),size(cols,1)),:)),hold on;
                xValsCorr=[xValsCorr jX{tr}(istim)];
                yValsCorr=[yValsCorr jY{tr}(istim)];
            elseif (condition<4 && outcome{tr}==2) || (condition>3 && outcome{tr}==1) %error trial
                subplot(a((ic-1)*4+3));
                %plot(jX{tr}(istim),jY{tr}(istim),'r'),hold on
                scatter(jX{tr}(istim),jY{tr}(istim),3,cols(min(1:sum(istim),size(cols,1)),:)),hold on;
                xValsErr=[xValsErr jX{tr}(istim)];
                yValsErr=[yValsErr jY{tr}(istim)];
            end

        end
        subplot(a((ic-1)*4+1))
        if ic==1
            ylabel('Joystick Y coordinate (AU)')
        end
        xlim([-72 72])
        ylim([-5 150])
        set(gca,'yticklabel',[],'xticklabel',[])
        title(['Correct Trials for condition ' num2str(ic)])
        subplot(a((ic-1)*4+3))
        if ic==1
            ylabel('Joystick Y coordinate (AU)')
        end
        set(gca,'yticklabel',[],'xticklabel',[])
        xlim([-70 70])
        ylim([-5 150])
        title('Error Trials')


        ctrs{1}=-70:4:70;ctrs{2}=-2.5:5:147.5;
        nCorr=hist3([xValsCorr' yValsCorr'],ctrs);
        if ~isempty(xValsErr) %do only if there are error trials
            nErr=hist3([xValsErr' yValsErr'],ctrs);
        else
            nErr=zeros(length(ctrs{1}),length(ctrs{2}));
        end
        subplot(a((ic-1)*4+2));
        imagesc(ctrs{1},ctrs{2},nCorr(1:end,end:-1:1)');
        set(gca,'yticklabel',[],'xticklabel',[])
        colormap('hot')
        caxis([0 Nsaturate])
        if ic==1
            ylabel('Joystick Y coordinate (AU)')
        end
        subplot(a((ic-1)*4+4));
        imagesc(ctrs{1},ctrs{2},nErr(1:end,end:-1:1)');
        set(gca,'yticklabel',[],'xticklabel',[])
        caxis([0 Nsaturate])
        xlabel('Joystick X coordinate (AU)')
        if ic==1
            ylabel('Joystick Y coordinate (AU)')
        end
        
    end
 
end











