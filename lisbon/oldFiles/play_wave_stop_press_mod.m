function [DecisionLed, DelayLed, xvalues_for_trial, response, numberOfSecondsElapsed, t1, t2, endPositionSecs, times_for_joystick, yvalues_for_trial, ball_coordinates, randT_prestim] = play_wave_stop_press_mod(data, Fs,stim_duration_sec,ramp_stim_dur_msec, pahandle, duration_of_trial,window, windowRect, black, ydir, dio, condition, current_trial, numberoftrials, ZeroPadding, ILD_vs_time, chunk_duration_msec, center_ILD, xmovement, ymovement, joydelay)


num_of_chunks_max = stim_duration_sec/(chunk_duration_msec/1000);
repetitions=1;

wave_L=data(1,:);
wave_R=data(2,:);

  time_aux=1/length(data):1/length(data):1;
  wave_3=sin(2*pi*30000*time_aux);
  
  wave_3(wave_3(1:end)<=0)=-5;
  wave_3(wave_3(1:end)>0)=5;
  
  sync_pulses =sin(2*pi*1000*time_aux);
  sync_pulse=zeros(1,length(time_aux));
  sync_pulse(1:44)=sync_pulses(1:44);
  sync_pulse(80+1:80+44)=sync_pulses(1:44);
  sync_pulse(440+1:440+44)=sync_pulses(1:44);
  
  
 wave_L=horzcat(ZeroPadding,wave_L);
 wave_R=horzcat(ZeroPadding,wave_R);
 wave_3=horzcat(ZeroPadding,wave_3);
 sync_pulse=horzcat(ZeroPadding,sync_pulse);

    % Read WAV file from filesystem:
    wavedata = [wave_L;wave_R; wave_3; sync_pulse];


% Perform basic initialization of the sound driver (done in main program):
% object pahandle

% Fill the audio playback buffer with the audio data 'wavedata':
PsychPortAudio('FillBuffer', pahandle, wavedata);

%%%we will stay in a loop until goal is reached
goal_reached=0;

%%%joystick action

% Here are the parameters for the ball to be displayed
spotRadius = 12.5; % The radius of the spot.
rotationRadius = 200; % The radius of the rotation.
initialRotationAngle = 3 * pi / 2; % The initial rotation angle in radians.

    % Use the parameters.
    spotDiameter = spotRadius * 1;
    spotRect = [0 0 spotDiameter spotDiameter];
    centeredspotRect = CenterRect(spotRect, windowRect); % Center the spot.
    
    centeredspotRect(2)=1.95*centeredspotRect(2);
    centeredspotRect(4)=centeredspotRect(2)+spotRadius;
    
    % Set up the timer.
    durationInSeconds = duration_of_trial;
    
    % Videogame loop variables
    
    xvalues_for_trial=zeros(1,500*floor(duration_of_trial));
    yvalues_for_trial=zeros(1,500*floor(duration_of_trial));
    times_for_joystick=zeros(1,500*floor(duration_of_trial));
    ball_coordinates = zeros(120*floor(duration_of_trial),4);
%     gradient_coordinates=zeros(120*duration_of_trial,4);
    
    %%%%%% graphics init params above

            boxwidth=50;
            boxseparation=170;
            
            leftBox_X1=(centeredspotRect(1)+centeredspotRect(3))/2-(boxseparation+boxwidth);
            leftBox_X2=leftBox_X1+boxwidth;
            leftBox_Y1=315;
            leftBox_Y2= leftBox_Y1+boxwidth;
            
            rightBox_X1=(centeredspotRect(1)+centeredspotRect(3))/2+boxseparation;
            rightBox_X2= rightBox_X1+boxwidth;
            rightBox_Y1=315;
            rightBox_Y2=rightBox_Y1+boxwidth;
            

numberOfSecondsElapsed=0;
randT_prestim = rand(1)+1.5;
tic

while toc<randT_prestim
            joystick_state=jst;
            ydir=joystick_state(1);
            xdir=joystick_state(2);
            
            if ydir<-.8
                ydir=-.8;
            end
            
                xmovement=xdir;
                ymovement=ydir-3.7;   
            xOffset = rotationRadius*xmovement;
            yOffset = rotationRadius*ymovement;
            offsetCenteredspotRect = OffsetRect(centeredspotRect, xOffset, yOffset);        
            Screen('FillRect', window, [128 128 128], [leftBox_X1 leftBox_Y1 leftBox_X2 leftBox_Y2]);
%             Screen('FillRect', window, [color_fill_2 color_fill_2 color_fill_2], [leftBox_X1 box_y1 leftBox_X2 box_y2]);

            Screen('FillRect', window, [128 128 128], [rightBox_X1 rightBox_Y1 rightBox_X2 rightBox_Y2]);
%             Screen('FillRect', window, [color_fill color_fill color_fill], [rightBox_X1 box_y1 rightBox_X2 box_y2]);

            if abs(xdir)>.04 || ydir<-.02
                ball_color=[255 0 0];
            else
                ball_color=[0 0 127];
            end
            
            Screen('FillOval', window, ball_color, offsetCenteredspotRect);
            Screen('Flip', window);
end


if abs(xdir)>.04 || ydir<.95
  
    
    putvalue(dio.line([1 4 6]),[1 1 1]);    %turn on joystick LED

     Screen('FillRect', window, [255 0 0], [leftBox_X1 leftBox_Y1 leftBox_X2 leftBox_Y2]);
 Screen('FillRect', window, [255 0 0], [rightBox_X1 rightBox_Y1 rightBox_X2 rightBox_Y2]);   
%  Screen('DrawText', window, sprintf('%g seconds remaining...', 0), 380, 500, black);
 Screen('Flip', window);  
[startTime endPositionSecs xruns t2] = PsychPortAudio('Stop', pahandle);

annoying_sound=wavread('Buzz.wav');

annoying_sound=.00002*annoying_sound(1:44100)'/max(annoying_sound(1:44100))*10^(60/20);

wave2=[annoying_sound'; annoying_sound'; zeros(size(annoying_sound))'; zeros(size(annoying_sound))'];

 PsychPortAudio('FillBuffer', pahandle,  wave2);
 PsychPortAudio('Start', pahandle, repetitions, 0, 1);
  pause(1)
 PsychPortAudio('Stop', pahandle);
 
 
DecisionLed = NaN;
DelayLed = NaN;
 response = 3;
 numberOfSecondsElapsed=NaN;
 t1 = NaN;
 t2=NaN;
 endPositionSecs=NaN;

 
else
       
    

t1 = PsychPortAudio('Start', pahandle, repetitions, 0, 1);

%Detect sound onset
value=0;
 while value==0     
 value = getvalue(dio.line9);
 SoundStart=tic;
 end
%Turn on Sync Led 
putvalue(dio.line([1 4 6]),[1 1 1]);    %turn on sync LEDs
DelayLed=toc(SoundStart);

i=1;
    while numberOfSecondsElapsed < durationInSeconds  && goal_reached==0  
     
        
        
%             TimeRemaining = round((durationInSeconds-numberOfSecondsElapsed)*100)/100;
%             Screen('DrawText', window, sprintf('%g seconds remaining...', TimeRemaining), 380, 500, black);

            joystick_state=jst;
            ydir=joystick_state(1);
            xdir=joystick_state(2);
            if ydir<-.8
                ydir=-.8;
            end
            xvalues_for_trial(i)=xdir;
            yvalues_for_trial(i)=ydir;
            
          
            
            if numberOfSecondsElapsed>joydelay && ydir>.90
                goal_reached=1;
            end
            
            if numberOfSecondsElapsed>.5
                
                if yvalues_for_trial(i-15)==yvalues_for_trial(i) && xvalues_for_trial(i-15)==xvalues_for_trial(i)
                    goal_reached=1;
                end
            end
            
                xmovement=xdir;
                ymovement=ydir-3.7;

            times_for_joystick(i)=toc(SoundStart);
%                 if abs(xdir)>.1
%                         xmovement = xdir;  %uncomment to make ball move
%                 end
                             
            xOffset = rotationRadius*xmovement;
            yOffset = rotationRadius*ymovement;
            offsetCenteredspotRect = OffsetRect(centeredspotRect, xOffset, yOffset);        


         if xdir<=-.9 && ydir<=-.7
               numberOfSecondsElapsed = toc(SoundStart);
                
               putvalue(dio.line([1 2 4 6]),[0 1 0 0]);    %turn on joystick LED 
                DecisionLed=toc(SoundStart);

               [startTime endPositionSecs xruns t2] = PsychPortAudio('Stop', pahandle);
               
          
               
                if mean(ILD_vs_time(1:min(floor((endPositionSecs-length(ZeroPadding)/Fs-ramp_stim_dur_msec/1000)/(chunk_duration_msec/1000)),num_of_chunks_max)))<center_ILD

                   t2=mean(ILD_vs_time(1:min(floor((endPositionSecs-length(ZeroPadding)/Fs-ramp_stim_dur_msec/1000)/(chunk_duration_msec/1000)),num_of_chunks_max)));
                    L_color= [0 255 0];
                    response=1;
                    
                    correct_sound=wavread('correct.wav');

correct_sound=.00002*correct_sound'/max(correct_sound)*10^(60/20);

wavecorrect=[correct_sound; correct_sound; zeros(size(correct_sound)); zeros(size(correct_sound))];

 PsychPortAudio('FillBuffer', pahandle,  wavecorrect);
 PsychPortAudio('Start', pahandle, repetitions, 0, 1);
  pause(1)
 PsychPortAudio('Stop', pahandle);
                    
                else
                    t2=mean(ILD_vs_time(1:min(floor((endPositionSecs-length(ZeroPadding)/Fs-ramp_stim_dur_msec/1000)/(chunk_duration_msec/1000)),num_of_chunks_max)));
                    L_color= [255 0 0];
                    response=2;
                    
                                                          wrong_sound=wavread('wrong.wav');

wrong_sound=.00002*wrong_sound'/max(wrong_sound)*10^(60/20);

wavewrong=[wrong_sound; wrong_sound; zeros(size(wrong_sound)); zeros(size(wrong_sound))];

 PsychPortAudio('FillBuffer', pahandle,  wavewrong);
 PsychPortAudio('Start', pahandle, repetitions, 0, 1);
  pause(1)
 PsychPortAudio('Stop', pahandle);
                end
                Screen('FillRect', window, [128 128 128], [rightBox_X1 rightBox_Y1 rightBox_X2 rightBox_Y2]);
%                 Screen('FillRect', window, L_color, [leftBox_X1 leftBox_Y1 leftBox_X2 leftBox_Y2]);
                Screen('FillRect', window, [128 128 128], [leftBox_X1 leftBox_Y1 leftBox_X2 leftBox_Y2]);
               
                Screen('Flip', window);

                
                goal_reached=1;
                
            elseif xdir>=.9 && ydir<=-.7
               numberOfSecondsElapsed = toc(SoundStart);
              
              putvalue(dio.line([1 2 4 6]),[0 1 0 0]);    %turn on joystick LED 
              DecisionLed=toc(SoundStart);
              
              [startTime endPositionSecs xruns t2] = PsychPortAudio('Stop', pahandle);
                if mean(ILD_vs_time(1:min(floor((endPositionSecs-length(ZeroPadding)/Fs-ramp_stim_dur_msec/1000)/(chunk_duration_msec/1000)),num_of_chunks_max)))>center_ILD
                  
                    t2=mean(ILD_vs_time(1:min(floor((endPositionSecs-length(ZeroPadding)/Fs-ramp_stim_dur_msec/1000)/(chunk_duration_msec/1000)),num_of_chunks_max)));
                    R_color=[0 255 0];
                    response=1;
                    
                  correct_sound=wavread('correct.wav');

correct_sound=.00002*correct_sound'/max(correct_sound)*10^(60/20);

wavecorrect=[correct_sound; correct_sound; zeros(size(correct_sound)); zeros(size(correct_sound))];

 PsychPortAudio('FillBuffer', pahandle,  wavecorrect);
 PsychPortAudio('Start', pahandle, repetitions, 0, 1);
  pause(1)
 PsychPortAudio('Stop', pahandle);
                    
                else
                    t2=mean(ILD_vs_time(1:min(floor((endPositionSecs-length(ZeroPadding)/Fs-ramp_stim_dur_msec/1000)/(chunk_duration_msec/1000)),num_of_chunks_max)));
                    R_color=[255 0 0];
                    response=2;
                    
                    wrong_sound=wavread('wrong.wav');

wrong_sound=.00002*wrong_sound'/max(wrong_sound)*10^(60/20);

wavewrong=[wrong_sound; wrong_sound; zeros(size(wrong_sound)); zeros(size(wrong_sound))];

 PsychPortAudio('FillBuffer', pahandle,  wavewrong);
 PsychPortAudio('Start', pahandle, repetitions, 0, 1);
  pause(1)
 PsychPortAudio('Stop', pahandle);
                end 
                Screen('FillRect', window, [128 128 128], [leftBox_X1 leftBox_Y1 leftBox_X2 leftBox_Y2]);
%                 Screen('FillRect', window, R_color, [rightBox_X1 rightBox_Y1 rightBox_X2 rightBox_Y2]);
                Screen('FillRect', window, [128 128 128], [rightBox_X1 rightBox_Y1 rightBox_X2 rightBox_Y2]);
                
                Screen('Flip', window);
 
                 goal_reached=1;
                 
                   
         end 
             
                if goal_reached==1
                    break
                end
                
            numberOfSecondsElapsed=toc(SoundStart);
%             color_fill=128-128*numberOfSecondsElapsed/duration_of_trial;    %gradually change bars color
%             color_fill_2=0+128*numberOfSecondsElapsed/duration_of_trial;    %gradually change bars color
% 
%             box_y1= leftBox_Y2-200*numberOfSecondsElapsed/duration_of_trial; %sand clock
%             box_y2= leftBox_Y2;
            
            Screen('FillRect', window, [128 128 128], [leftBox_X1 leftBox_Y1 leftBox_X2 leftBox_Y2]);
%             Screen('FillRect', window, [color_fill_2 color_fill_2 color_fill_2], [leftBox_X1 box_y1 leftBox_X2 box_y2]);

            Screen('FillRect', window, [128 128 128], [rightBox_X1 rightBox_Y1 rightBox_X2 rightBox_Y2]);
%             Screen('FillRect', window, [color_fill color_fill color_fill], [rightBox_X1 box_y1 rightBox_X2 box_y2]);

            Screen('FillOval', window, [0 0 127], offsetCenteredspotRect);  %change to offsetCenteredspotRect
            
            ball_coordinates(i, :) = offsetCenteredspotRect;  %change to offsetCenteredspotRect if you want the ball to move
            
            Screen('Flip', window);
            
            i=i+1;
    end
    
 if exist('response', 'var')
     %subject responded
 else
         putvalue(dio.line2,1);  %turn on goal reached Led
         response=3;  %trial was an abort/no answer
 end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
     
if response==3;
 DecisionLed=toc(SoundStart);   
 Screen('FillRect', window, [255 0 0], [leftBox_X1 leftBox_Y1 leftBox_X2 leftBox_Y2]);
 Screen('FillRect', window, [255 0 0], [rightBox_X1 rightBox_Y1 rightBox_X2 rightBox_Y2]);   
%  Screen('DrawText', window, sprintf('%g seconds remaining...', 0), 380, 500, black);
 Screen('Flip', window);  
[startTime endPositionSecs xruns t2] = PsychPortAudio('Stop', pahandle);

annoying_sound=wavread('Buzz.wav');

annoying_sound=.00002*annoying_sound(1:44100)'/max(annoying_sound(1:44100))*10^(60/20);

wave2=[annoying_sound'; annoying_sound'; zeros(size(annoying_sound))'; zeros(size(annoying_sound))'];

 PsychPortAudio('FillBuffer', pahandle,  wave2);
 PsychPortAudio('Start', pahandle, repetitions, 0, 1);
  pause(1)
 PsychPortAudio('Stop', pahandle);

else
%do nothing
end

end

 putvalue(dio.line([1 2 4 6]),[0 0 0 0]);  %turn off event leds

% Close the audio device:
PsychPortAudio('DeleteBuffer');



