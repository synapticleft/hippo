function Data=data_decode
%% Reads the csv files from bonsai and eye tracking data

% 'game_data': data from game
%
% 'bonsai_data': data structure from csv files outputed by bonsai
%     'start_time': time in seconds corresponding to start button release time
%     'action_onset' matrix (binary). Order of columns: right, left, bank, cash, timespan from start_time, wave number
%     'pressure' (analog). same order as above
%     'eye_position': coordinates of eye tracking data
%     'pupil_area': area of pupil
%     
%  'calibration': calibration data for the eye tracking
    

%%

gamedata=getKeyStrokesConstruct;

Data.game_data.name=gamedata.users{1,1}.name;
time=gamedata.users{1,1}.session{1,1}.time;
%time=datestr(719529+time/86400000,'HH:MM:SS.FFF');
%time=hour(time)*3600+minute(time)*60+second(time);
Data.game_data.gametime=time;
Data.game_data.numVertices=gamedata.users{1,1}.session{1,1}.numVertices;
Data.game_data.events=gamedata.users{1,1}.session{1,1}.event;


csv_files=dir('*.csv');

for i=1:length(csv_files)
    
    if strncmp(csv_files(i).name,'acc',3)   %accelerometers. left_up
    fileID = fopen(csv_files(i).name);
    acc = textscan(fileID,'%f %u-%u-%uT%f:%f:%f+00:00 %f %u-%u-%uT%f:%f:%f+00:00');
    fclose(fileID);
    end

    if strncmp(csv_files(i).name,'bank_cash',3) %Bank_cash. value_timestamp_state
    fileID = fopen(csv_files(i).name);
    bank_cash = textscan(fileID,'%f %u-%u-%uT%f:%f:%f+00:00 %s %f %u-%u-%uT%f:%f:%f+00:00 %s');
    fclose(fileID);
    end

    if strncmp(csv_files(i).name,'calib',3) %Calibration start. state_timestamp
    fileID = fopen(csv_files(i).name);
    calib_start = textscan(fileID,'%s %u-%u-%uT%f:%f:%f+00:00');
    fclose(fileID);
    end

    if strncmp(csv_files(i).name,'eye',3) % eye. timestamp
    fileID = fopen(csv_files(i).name);
    eye = textscan(fileID,'%u-%u-%uT%f:%f:%f+00:00');
    fclose(fileID);
    end

    if strncmp(csv_files(i).name,'right_left',3) % Right Left. value_timestamp_state
    fileID = fopen(csv_files(i).name);
    right_left = textscan(fileID,'%f %u-%u-%uT%f:%f:%f+00:00 %s %f %u-%u-%uT%f:%f:%f+00:00 %s');
    fclose(fileID);
    end

    if strncmp(csv_files(i).name,'start',3) %state timestamp start. state_timestamp
    fileID = fopen(csv_files(i).name);
    start = textscan(fileID,'%s %u-%u-%uT%f:%f:%f+00:00');
    fclose(fileID);
    end

    if strncmp(csv_files(i).name,'world',3) %world camera. timestamps
    fileID = fopen(csv_files(i).name);
    world = textscan(fileID,'%u-%u-%uT%f:%f:%f+00:00');
    fclose(fileID);
    end

    if strncmp(csv_files(i).name,'pupil_pos',3) %pupil centroid
    fileID=fopen(csv_files(i).name);
    pupil=textscan(fileID,'%f %f');
    fclose(fileID);
    end

    if strncmp(csv_files(i).name,'marker_pos',3) %calibration marker
    fileID=fopen(csv_files(i).name);
    marker=textscan(fileID,'%f %f');
    fclose(fileID);
    end
end

%%

start_pressed=strcmp(start{1,1}(1:length(start{1,1})),'Pressed');
onsets=diff([start_pressed; 0]==0);
start_point=find(onsets==1);
start_time=start{1,5}(start_point,1)*3600+start{1,6}(start_point,1)*60+start{1,7}(start_point,1);



calib_pressed=strcmp(calib_start{1,1}(1:length(calib_start{1,1})),'Pressed');
onsets=diff([calib_pressed; 0]==0);
calib_point=find(onsets==1);
calib_time=calib_start{1,5}(calib_point,1)*3600+calib_start{1,6}(calib_point,1)*60+calib_start{1,7}(calib_point,1);

if length(calib_point)==1
    calib_time=[calib_time-start_time calib_time+1000-start_time];
else
    calib_time=[calib_time-start_time];
end

action_onset= zeros(length(right_left{1,1}),5);
action_onset(:,1)=strncmp('Pressed',right_left{1,8},5); %right
action_onset(:,2)=strncmp('Pressed',right_left{1,16},5); %left
action_onset(:,3)=strncmp('Pressed',bank_cash{1,8},5); %bank
action_onset(:,4)=strncmp('Pressed',bank_cash{1,16},5); %cash
action_onset(:,5)=(right_left{1,5}*3600+right_left{1,6}*60+right_left{1,7})-repmat(start_time,length(action_onset),1); %Timespan
action_onset(action_onset(:,5)<0,:)=[];
action_onset(action_onset(:,5)>=calib_time(2),:)=[];

pressure=zeros(length(right_left{1,1}),5);
pressure(:,1)=right_left{1,1};
pressure(:,2)=right_left{1,9};
pressure(:,3)=bank_cash{1,1};
pressure(:,4)=bank_cash{1,9};
pressure(:,5)=(right_left{1,5}*3600+right_left{1,6}*60+right_left{1,7})-repmat(start_time,length(pressure),1);
pressure(pressure(:,5)<0,:)=[];
pressure(pressure(:,5)>=calib_time(2),:)=[];

eye_move=zeros(length(eye{1,1}),3);
eye_move(:,3)=eye{1,4}*3600+eye{1,5}*60+eye{1,6}; %timestamps


world_sec=world{1,4}*3600+world{1,5}*60+world{1,6};


Data.bonsai_data.start_time=start_time;
Data.bonsai_data.action_onset=action_onset;
Data.bonsai_data.pressure=pressure;
Data.bonsai_data.eye_position=eye_move;
Data.bonsai_data.pupil_area={};
Data.calibration={};

%% Calibration data

% eye_calib_ind=find(eye_sec>=calib_time & eye_sec<start_time);
% eye_calib_timespan=eye_sec(eye_calib_ind)-repmat(calib_time,length(eye_sec(eye_calib_ind)),1);
% eye_calib=[pupil{1,1}(eye_calib_ind) pupil{1,2}(eye_calib_ind) eye_calib_timespan];
% 
% world_calib_ind=find(world_sec>=calib_time & world_sec<start_time);
% world_calib_timespan=world_sec(world_calib_ind)-repmat(calib_time,length(world_sec(world_calib_ind)),1);
% world_calib=[marker{1,1}(world_calib_ind) marker{1,2}(world_calib_ind) world_calib_timespan];

end



