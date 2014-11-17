function []=TrialExplorerVertical(filename,TrialType)

% This function plots all variables measured by Roberto in version2 of his
% human psychophysics experiment. The user specifies the dataset in
% filename, and can specify trial types and a session (optional, defaults
% to all trials and all sessions in chronological order). The specification
% requires a 2 element array. First element is session (1 to 12), second element is
% trial type (0.5, 2, 4 or 6). If either is equal to zero, only the other
% one is specified. To scroll across trials, press space, to stop, press
% control C. 

load(filename)

% This loads the cell array with all the data interpolated at 80 Hz (2 sec
% of baseline plus the stimulus plus 2 sec poststim), a scalar center_ILD,
% with the ILD center for the subject (in dB SPL) and the array
% bin_centers, which contains the ILDs that can be used in the x-axis of
% the psychometric curve.

% Example usage
%
% TrialExplorerVertical('data_s1_interp_2sPrePost_wRaw',[11 0.5]) ; analyze data in
% from trials with ILD plus or minus 0.5 from session 11 of the
% corresponding experiment.
%
% TrialExplorerVertical('data_s1_interp_2sPrePost_wRaw',[4 0]) ; same for all trial
% types in session 4.
%

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

% Search for Specific trials
if nargin>1 % Select subset of all trials
    if length(TrialType)~=2
        disp(' ')
        disp('WATCH OUT!! SESSION/TRIAL SPECIFICATION NEEDS A 2 ELEMENT ARRAY !!')
        disp(' ')
        return
    end    
    ses=TrialType(1);
    trial=TrialType(2);
    if ses>0 & trial==0  % Only session specified
        disp(' ')
        disp(['Selecting all trials from session = ' num2str(ses)])
        disp(' ')
        session{1}=[];
        sesMat=cell2mat(session);
        trV=find(sesMat==ses)+1;
        if isempty(trV)
            disp(' ')
            disp('WATCH OUT!! SESSION SPECIFICATION FAILED!!!')
            disp(' ')
            return
        end
    elseif ses>0 & trial>0 % Specify both session and trial type
        disp(' ')
        disp(['Selecting trials with abs(ILD) = ' num2str(trial) ' from session = ' num2str(ses)])
        disp(' ')
        session{1}=[];
        ILD_th{1}=[];        
        sesMat=cell2mat(session);
        ILD_thMat=cell2mat(ILD_th);
        trV=find(sesMat==ses & abs(ILD_thMat)==trial)+1;
        if isempty(trV)
            disp(' ')
            disp('WATCH OUT!! SESSION + TRIAL SPECIFICATION FAILED!!!')
            disp(' ')
            return
        end        
    elseif ses==0 & trial>0 % Specify trial type and no particular session
        disp(' ')
        disp(['Selecting trials with abs(ILD) = ' num2str(trial)])
        disp(' ')
        ILD_th{1}=[];    
        ILD_thMat=cell2mat(ILD_th);
        trV=find(abs(ILD_thMat)==trial)+1;
        if isempty(trV)
            disp(' ')
            disp('WATCH OUT!! TRIAL SPECIFICATION FAILED!!!')
            disp(' ')
            return
        end         
    end
else % use all trials
    disp(' ')
    disp('Using all trials and sessions')
    trV=2:trialN{end}; 
end

%Plot results
figure('position',[1214         163         614        1183])

% Figure Params
Nrows=9;
Ncols=1;
delR=0.01;
delC=0.04;
pX=(0.85-(Ncols-1)*delC)/Ncols;
pY=(0.9-(Nrows-1)*delR)/Nrows;
tinit=-0.5;

% Plot 1st trial
i=1;
tr=trV(i); %trial number 
tN=trialN{tr};
sess=session{tr};
ILDth=ILD_th{tr};
outc=outcome{tr};
tit=['Trial ' num2str(tN) ' ; Session ' num2str(sess) ' ; ILD LongTerm = ' num2str(ILDth) ' ; outcome = ' num2str(outc)];


% Plot stim
a(1)=axes('position',[0.1 0.95-pY pX pY]);
istim=~isnan(stim{tr});
t=time{tr};
timeSt=time{tr}(istim);    
p(1)=plot(t,stim{tr}-center_ILD,'b');hold on
p(2)=plot(timeSt,mean(stim{tr}(istim))*ones(1,length(timeSt)),'g');
plot([tinit 10],[0 0],'k')
xlim([tinit ceil(timeSt(end))])
yl=[-12 12];
ylim(yl);
p(3)=plot([timeSt(1) timeSt(1)],yl,'r');
p(4)=plot([timeSt(end) timeSt(end)],yl,'r');
set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])
set(gca,'ytick',[-8 -4 0 4 8])
ylabel('ILD (dB SPL)')
titl(1)=title(tit,'fontsize',14);
grid on

% Plot joystick and accelerometer X
a(4)=axes('position',[0.1 0.95-2*pY-delR pX pY]);    
p(13)=plot(t,jX{tr},'b');hold on
p(14)=plot(t,50*accX{tr},'c');
legend('joystick','accelerometer',2)
xlim([tinit ceil(timeSt(end))])
yl=[-80 80];
ylim(yl)
p(15)=plot([timeSt(1) timeSt(1)],yl,'r');
p(16)=plot([timeSt(end) timeSt(end)],yl,'r');
set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])
set(gca,'ytick',[-60 -30 0 30 60])
ylabel('Joystick and Acc X (AU)')
grid on

% Plot joystick and accelerometer Y
a(5)=axes('position',[0.1 0.95-3*pY-2*delR pX pY]);    
p(17)=plot(t,jY{tr},'b');hold on
p(18)=plot(t,50*accY{tr},'c');
xlim([tinit ceil(timeSt(end))])
yl=[-75 150];
ylim(yl)
p(19)=plot([timeSt(1) timeSt(1)],yl,'r');
p(20)=plot([timeSt(end) timeSt(end)],yl,'r');
set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])
set(gca,'ytick',[-100 -50 0 50 100])
ylabel('Joystick and Acc Y (AU)') 
grid on

% Plot Pupil Center X
a(9)=axes('position',[0.1 0.95-4*pY-3*delR pX pY]);    
p(30)=plot(t,pupX{tr},'c');hold on
p(52)=plot(pupT{tr},pupXr{tr},'bo');
legend('interp/smooth','raw',2)
xlim([tinit ceil(timeSt(end))])
yl=get(gca,'ylim');
ylim(yl)
p(31)=plot([timeSt(1) timeSt(1)],yl,'r');
p(32)=plot([timeSt(end) timeSt(end)],yl,'r');
set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])
ylabel('Pupil Center X (AU)')
grid on

% Plot Pupil Center Y
a(10)=axes('position',[0.1 0.95-5*pY-4*delR pX pY]);    
p(33)=plot(t,pupY{tr},'c');hold on
p(53)=plot(pupT{tr},pupYr{tr},'bo');
xlim([tinit ceil(timeSt(end))])
yl=get(gca,'ylim');
ylim(yl)
p(34)=plot([timeSt(1) timeSt(1)],yl,'r');
p(35)=plot([timeSt(end) timeSt(end)],yl,'r');
set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])
ylabel('Pupil Center Y (AU)')
grid on

% Plot Pupil Major Axis
a(11)=axes('position',[0.1 0.95-6*pY-5*delR pX pY]);    
p(36)=plot(t,pupMaj{tr},'c');hold on
p(54)=plot(pupT{tr},pupMajr{tr},'bo');
xlim([tinit ceil(timeSt(end))])
yl=get(gca,'ylim');
ylim(yl)
p(37)=plot([timeSt(1) timeSt(1)],yl,'r');
p(38)=plot([timeSt(end) timeSt(end)],yl,'r');
set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])
ylabel('Pupil Major Axis (AU)')
grid on

% Plot Pupil Minor Axis
a(12)=axes('position',[0.1 0.95-7*pY-6*delR pX pY]);    
p(39)=plot(t,pupMin{tr},'c');hold on
p(55)=plot(pupT{tr},pupMinr{tr},'bo');
xlim([tinit ceil(timeSt(end))])
yl=get(gca,'ylim');
ylim(yl)
p(40)=plot([timeSt(1) timeSt(1)],yl,'r');
p(41)=plot([timeSt(end) timeSt(end)],yl,'r');
set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])
ylabel('Pupil Minor Axis (AU)') 
grid on

% Plot emg 1 Pronator
a(14)=axes('position',[0.1 0.95-8*pY-7*delR pX pY]);    
p(45)=plot(t,musc1{tr},'b');hold on
xlim([tinit ceil(timeSt(end))])
yl=[0 2.5];
ylim(yl)
p(46)=plot([timeSt(1) timeSt(1)],yl,'r');
p(47)=plot([timeSt(end) timeSt(end)],yl,'r');
set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])
set(gca,'ytick',[0 1 2])
ylabel('EMG Pronator (AU)') 
grid on

% Plot emg 2 Bicep
a(15)=axes('position',[0.1 0.95-9*pY-8*delR pX pY]);    
p(48)=plot(t,musc2{tr},'b');hold on
xlim([tinit ceil(timeSt(end))])
yl=[0.5 4.5];
ylim(yl)
p(49)=plot([timeSt(1) timeSt(1)],yl,'r');
p(50)=plot([timeSt(end) timeSt(end)],yl,'r');
set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[tinit:0.25:ceil(timeSt(end))])
set(gca,'ytick',[1 2 3 4])
ylabel('EMG Bicep (AU)') 
xlabel('Time (s)')
grid on

pause

% Update Plots for subsequent trials
for i=2:length(trV) %loop over trials
    tr=trV(i); %trial number 
    tN=trialN{tr};
    sess=session{tr};
    ILDth=ILD_th{tr};
    outc=outcome{tr};
    tit=['Trial ' num2str(tN) ' ; Session ' num2str(sess) ' ; ILD LongTerm = ' num2str(ILDth) ' ; outcome = ' num2str(outc)];
    istim=~isnan(stim{tr});
    t=time{tr};
    timeSt=time{tr}(istim);
            
    % Plot stim
    subplot(a(1));    
    set(p(1),'xdata',t);
    set(p(1),'ydata',stim{tr}-center_ILD);    
    set(p(2),'xdata',timeSt);
    set(p(2),'ydata',mean(stim{tr}(istim))*ones(1,length(timeSt)));    
    yl=[-12 12];
    xlim([tinit ceil(timeSt(end))])
    set(p(3),'xdata',[timeSt(1) timeSt(1)]);
    set(p(3),'ydata',yl);    
    set(p(4),'xdata',[timeSt(end) timeSt(end)]);
    set(p(4),'ydata',yl);    
    set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])
    set(titl(1),'string',tit);
    
    % Plot joystick and accelerometer X
    subplot(a(4));    
    set(p(13),'xdata',t);
    set(p(13),'ydata',jX{tr});    
    set(p(14),'xdata',t);
    set(p(14),'ydata',50*accX{tr});
    yl=[-80 80];
    xlim([tinit ceil(timeSt(end))])
    set(p(15),'xdata',[timeSt(1) timeSt(1)]);
    set(p(15),'ydata',yl);    
    set(p(16),'xdata',[timeSt(end) timeSt(end)]);
    set(p(16),'ydata',yl);    
    set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])
    
    % Plot joystick and accelerometer Y
    subplot(a(5));    
    set(p(17),'xdata',t);
    set(p(17),'ydata',jY{tr});    
    set(p(18),'xdata',t);
    set(p(18),'ydata',50*accY{tr});
    yl=[-75 150];
    xlim([tinit ceil(timeSt(end))])
    set(p(19),'xdata',[timeSt(1) timeSt(1)]);
    set(p(19),'ydata',yl);    
    set(p(20),'xdata',[timeSt(end) timeSt(end)]);
    set(p(20),'ydata',yl);    
    set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])    

    % Plot Pupil Center X
    subplot(a(9));
    set(p(30),'xdata',t);
    set(p(30),'ydata',pupX{tr}); 
    set(p(52),'xdata',pupT{tr});
    set(p(52),'ydata',pupXr{tr}); 
    xlim([tinit ceil(timeSt(end))])
    yl=[min(pupX{tr})*0.95 max(pupX{tr})*1.05];
    ylim(yl)
    set(p(31),'xdata',[timeSt(1) timeSt(1)]);
    set(p(31),'ydata',yl);    
    set(p(32),'xdata',[timeSt(end) timeSt(end)]);
    set(p(32),'ydata',yl); 
    set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])    

    % Plot Pupil Center Y
    subplot(a(10));    
    set(p(33),'xdata',t);
    set(p(33),'ydata',pupY{tr}); 
    set(p(53),'xdata',pupT{tr});
    set(p(53),'ydata',pupYr{tr});
    xlim([tinit ceil(timeSt(end))])
    yl=[min(pupY{tr})*0.95 max(pupY{tr})*1.05];
    ylim(yl)
    set(p(34),'xdata',[timeSt(1) timeSt(1)]);
    set(p(34),'ydata',yl);    
    set(p(35),'xdata',[timeSt(end) timeSt(end)]);
    set(p(35),'ydata',yl); 
    set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])      

    % Plot Pupil Major Axis
    subplot(a(11));   
    set(p(36),'xdata',t);
    set(p(36),'ydata',pupMaj{tr}); 
    set(p(54),'xdata',pupT{tr});
    set(p(54),'ydata',pupMajr{tr});
    xlim([tinit ceil(timeSt(end))])
    yl=[min(pupMaj{tr})*0.95 max(pupMaj{tr})*1.05];
    ylim(yl)
    set(p(37),'xdata',[timeSt(1) timeSt(1)]);
    set(p(37),'ydata',yl);    
    set(p(38),'xdata',[timeSt(end) timeSt(end)]);
    set(p(38),'ydata',yl); 
    set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])  

    % Plot Pupil Minor Axis
    subplot(a(12)); 
    set(p(39),'xdata',t);
    set(p(39),'ydata',pupMin{tr}); 
    set(p(55),'xdata',pupT{tr});
    set(p(55),'ydata',pupMinr{tr});
    xlim([tinit ceil(timeSt(end))])
    yl=[min(pupMin{tr})*0.95 max(pupMin{tr})*1.05];
    ylim(yl)
    set(p(40),'xdata',[timeSt(1) timeSt(1)]);
    set(p(40),'ydata',yl);    
    set(p(41),'xdata',[timeSt(end) timeSt(end)]);
    set(p(41),'ydata',yl); 
    set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])         

    % Plot emg 1 Pronator
    subplot(a(14));
    set(p(45),'xdata',t);
    set(p(45),'ydata',musc1{tr}); 
    xlim([tinit ceil(timeSt(end))])
    yl=[0 2.5];
    set(p(46),'xdata',[timeSt(1) timeSt(1)]);
    set(p(46),'ydata',yl);    
    set(p(47),'xdata',[timeSt(end) timeSt(end)]);
    set(p(47),'ydata',yl); 
    set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[])  

    % Plot emg 2 Bicep
    subplot(a(15));
    set(p(48),'xdata',t);
    set(p(48),'ydata',musc2{tr}); 
    xlim([tinit ceil(timeSt(end))])
    yl=[0.5 4.5];
    set(p(49),'xdata',[timeSt(1) timeSt(1)]);
    set(p(49),'ydata',yl);    
    set(p(50),'xdata',[timeSt(end) timeSt(end)]);
    set(p(50),'ydata',yl); 
    set(gca,'xtick',[tinit:0.25:ceil(timeSt(end))],'xticklabel',[tinit:0.25:ceil(timeSt(end))])
    
    pause
end



























