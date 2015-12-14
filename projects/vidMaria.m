function vidMaria(moviename)
% This demo accepts a pattern for a valid moviename, e.g.,
% moviename=`*.avi`, then it plays all movies in the current working
% directory whose names match the provided pattern, e.g., the `*.mpg`
% pattern would play all MPEG files in the current directory.

% This demo introduces distortions onto an ongoing video stream that are
% scaled by the volume of the sound picked up by the microphone.

% Pressing the Cursor-Up/Down key pauses/unpauses the
% movie and increases/decreases playback rate.
% The right arrow key jumps in 1 seconds steps. SPACE jumps to the
% next movie in the list. ESC ends the demo.
% qwerty toggles different effects.

global accumInd     % keeps track of elapsed time
global state        % stores a video frame and any transforms applied to it
global lastAudio    % memory of audio history, used to detect sudden volume increase (aka beats)
global transform    % which transform is currrently used?
global shrink         % 2-d spatial gradient used to warp image (tunnel effect)
global arDecay
global arRate
global pNew
global pixFilt
global slowDown

transform = 2;
pNew = .01;arDecay = .01;arRate = pi/4;pixFilt = 10;slowDown = 20;

KbName('UnifyKeyNames');
accumInd = 0;
lastAudio = 0;      
gain = 1;           %controls sensitivity to sound
audiodata = 0;      %mean volume of current audio sample
%numPics = 10;
%weights = randn(numPics);
% Setup key mapping:
space=KbName('SPACE');
esc=KbName('ESCAPE');
right=KbName('RightArrow');
left=KbName('LeftArrow');
up=KbName('UpArrow');
down=KbName('DownArrow');
shift=KbName('RightShift');

ts = 'qwerty';
for i = 1:numel(ts)
    tforms(i) = KbName(ts(i));
end
ts = 'azsxdcfvgb';
for i = 1:numel(ts)
    mods(i) = KbName(ts(i));
end


if ~exist('moviename','var') || isempty(moviename)
    moviename = '*.avi';
end

% Initialize with unified keynames and normalized colorspace:
PsychDefaultSetup(2);
InitializePsychSound;

% Open the default audio device [], with mode 2 (== Only audio capture),
% and a required latencyclass of zero 0 == no low-latency mode, as well as
% a frequency of 44100 Hz and 2 sound channels for stereo capture.
% This returns a handle to the audio device:
freq = 44100;
pahandle = PsychPortAudio('Open', [], 2, 0, freq, 1);
% Preallocate an internal audio recording  buffer with a capacity of 10 seconds:
PsychPortAudio('GetAudioData', pahandle, 10);
% Start audio capture immediately and wait for the capture to start.
% We set the number of 'repetitions' to zero,
% i.e. record until recording is manually stopped.
PsychPortAudio('Start', pahandle, 0, 0, 1);

try
    % Open onscreen window with gray background:
    screen = max(Screen('Screens'));
    scrsz = get(screen,'ScreenSize');
    win = PsychImaging('OpenWindow', screen, [0, 0, 0]);

    % Initial display and sync to timestamp:
    Screen('Flip',win);
    iteration = 0;
    abortit = 0;
    
    % Use blocking wait for new frames by default:
    blocking = 1;
    
    % Return full list of movie files from directory+pattern:
    moviefiles=dir(moviename);
    
    for i=1:size(moviefiles,1)
        moviefiles(i).name = [ pwd filesep moviefiles(i).name ];
    end
    moviecount = size(moviefiles,1);
    
    % Playbackrate defaults to 1:
    rate=1;
    
    % Endless loop, runs until ESC key pressed:
    while (abortit<2)
        moviename = moviefiles(mod(iteration, moviecount)+1).name;
        iteration = iteration + 1;
        
        % Open movie file and retrieve basic info about movie:
        [movie, movieduration, fps, imgw, imgh] = Screen('OpenMovie', win, moviename, [], [], 2, [], []);
        
        % reset state and shrink to match the dimensions of new video
        state = zeros(imgh,imgw,3);
        [dx dy] = meshgrid(1:imgw,1:imgh);
        dx = dx-mean(dx(:));dy = dy-mean(dy(:));
        shrink = [];shrink(:,:,1) = -dx/10;shrink(:,:,2) = -dy/10;
        fprintf('Movie: %s  : %f seconds duration, %f fps, w x h = %i x %i...\n', moviename, movieduration, fps, imgw, imgh);
        
        % Start playback of movie. This will start
        % the realtime playback clock and playback of audio tracks, if any.
        % Play 'movie', at a playbackrate = 1, with endless loop=1.
        Screen('PlayMovie', movie, rate, 1);
        tic;ti = Screen('GetMovieTimeIndex', movie);
        % Infinite playback loop: Fetch video frames and display them...
        while 1
            % Check for abortion:
            abortit=0;
            [keyIsDown,secs,keyCode]=KbCheck;
            if (keyIsDown==1 && keyCode(esc))
                % Set the abort-demo flag.
                abortit=2;
                break;
            end;
            
            % Check for skip to next movie:
            if (keyIsDown==1 && keyCode(space))
                % Exit while-loop: This will load the next movie...
                break;
            end;
            for i = 1:numel(tforms)
                if (keyIsDown == 1 && keyCode(tforms(i)))
                    transform = i;
                end
            end
            
            % Only perform video image fetch/drawing if playback is active
            % and the movie actually has a video track (imgw and imgh > 0):
            if ((abs(rate)>0) && (imgw>0) && (imgh>0))
                % Retrieve pending audio data from the drivers internal ringbuffer:
                temp = mean(abs(PsychPortAudio('GetAudioData', pahandle)),2);
                if ~isempty(temp)
                    lastAudio = .95*lastAudio + .05*audiodata;
                    audiodata = temp;
                end
                if toc < slowDown
                    Screen('SetMovieTimeIndex', movie, ti);
                else
                    tic;ti = Screen('GetMovieTimeIndex', movie);
                end
                % Return next frame in movie, in sync with current playback
                % time and sound.
                % tex is either the positive texture handle or zero if no
                % new frame is ready yet in non-blocking mode (blocking == 0).
                % It is -1 if something went wrong and playback needs to be stopped:
                tex = Screen('GetMovieImage', win, movie, blocking);
                % Valid texture returned?
                if tex < 0
                    % No, and there won't be any in the future, due to some
                    % error. Abort playback loop:
                    break;
                end
                
                if tex == 0
                    % No new frame in polling wait (blocking == 0). Just sleep
                    % a bit and then retry.
                    WaitSecs('YieldSecs', 0.005);
                    continue;
                end
                % Draw the new texture immediately to screen:
                if transform == 1
                    jitterback(movie,audiodata,gain,movieduration );
                else
                    I = Screen('GetImage',tex);
                    Screen('Close',tex);
                    if transform == 2
                        I = washCycle(I,gain);% smearUpdate(I,audiodata,gain);%
                    elseif transform == 3
                        I = cosColor(I,audiodata,gain);
                    elseif transform == 4
                        I = shiftss(I,audiodata,gain);
                    elseif transform == 5
                        I = shiftRGB(I,audiodata,gain);
                    elseif transform == 6
                        I = tunnel(I,audiodata,gain);
                    end
                    tex = Screen('MakeTexture',win,I);
                end
                Screen('DrawTexture', win, tex,[],scrsz);
                
                % Update display:
                Screen('Flip', win);
                % Release texture:
                Screen('Close', tex);
            end;
            
            % Further keyboard checks...
            
            if (keyIsDown==1 && keyCode(right))
                % Advance movietime by one second:
                Screen('SetMovieTimeIndex', movie, Screen('GetMovieTimeIndex', movie) + 1);
            end;
            
            if (keyIsDown==1 && keyCode(left))
                gain = 1;pNew = .01;arDecay = .01;arRate = pi/4;pixFilt = 10;slowDown = 20;
            end;
            
            if (keyIsDown==1 && keyCode(up))
                gain = gain * 1.2;
            end;
            
            if (keyIsDown==1 && keyCode(down))
                gain = gain / 1.2;
            end;
            
            if (keyIsDown==1 && keyCode(mods(1)))
                arDecay = arDecay * 1.2;
            end;
            if (keyIsDown==1 && keyCode(mods(2)))
                arDecay = arDecay / 1.2;
            end;
            if (keyIsDown==1 && keyCode(mods(3)))
                arRate = arRate * 1.2;
            end;
            if (keyIsDown==1 && keyCode(mods(4)))
                arRate = arRate / 1.2;
            end;
            if (keyIsDown==1 && keyCode(mods(5)))
                pNew = pNew * 1.2;
            end;
            if (keyIsDown==1 && keyCode(mods(6)))
                pNew = pNew / 1.2;
            end;
            if (keyIsDown==1 && keyCode(mods(7)))
                pixFilt = pixFilt * 1.2;
            end;
            if (keyIsDown==1 && keyCode(mods(8)))
                pixFilt = pixFilt / 1.2;
            end;
            if (keyIsDown==1 && keyCode(mods(9)))
                slowDown = slowDown * 1.2;
            end;
            if (keyIsDown==1 && keyCode(mods(10)))
                slowDown = slowDown / 1.2;
            end;
        end;
        
        Screen('Flip', win);
        KbReleaseWait;
        
        % Done. Stop playback:
        Screen('PlayMovie', movie, 0);
        
        % Close movie object:
        Screen('CloseMovie', movie);
    end;
    
    % Close screens.
    Screen('CloseAll');
    % Stop capture:
    PsychPortAudio('Stop', pahandle);
    % Close the audio device:
    PsychPortAudio('Close', pahandle);
    % Done.
    return;
catch
    % Error handling: Close all windows and movies, release all ressources.
    Screen('CloseAll');
    % Stop capture:
    PsychPortAudio('Stop', pahandle);
    % Close the audio device:
    PsychPortAudio('Close', pahandle);

end

function jitterback(movie,ind,gain,movieduration )
global lastAudio
Screen('SetMovieTimeIndex', movie, mod(Screen('GetMovieTimeIndex', movie) + gain*(ind-lastAudio),movieduration));%-lastAudio


function dat = shiftss(dat,ind,gain)
[x, y] = gradient(double(rgb2gray(dat)));
x(:,:,2) = y;
dat = imwarp(dat,x*ind*gain);
clear x y;

function dat = washCycle(dat,gain)
global state
global accumInd
global arDecay
global arRate
global pNew
global pixFilt

accumInd = accumInd + 1;
%pNew = 2;
%pixFilt = 30;
temp = double(dat(mod(accumInd,size(dat,1))+1,:,1))/255 < pNew;%rand(1,size(dat,2)) < pNew; %
ar = (1-arDecay)*exp(1i*arRate);
lineInfo = idealfilterG(temp,pixFilt);%1;%
state = state*ar + bsxfun(@times,double(dat)/255,lineInfo/max(.01,max(lineInfo)));
dat = max(0,min(1,real(state)*gain));

function dat = shiftRGB(dat,ind,gain)
sh = round(ind*gain);
dat(:,:,1) = circshift(dat(:,:,1),[sh 0 0]);
dat(:,:,2) = circshift(dat(:,:,2),[0 sh 0]);
dat(:,:,3) = circshift(dat(:,:,3),[-sh -sh 0]);

function dat = smearUpdate(dat,ind,gain)
global state
state = state*.9;
indbool = state < ind*gain;
state(indbool) = dat(indbool);
dat = state;

function dat = tunnel(dat,ind,gain)
global state
global shrink
norm = min(1,ind*gain);
%state = max(.9*imwarp(state,norm*shrink),double(dat)/255);
state = (norm)*imwarp(state,shrink) +(1-norm)*double(dat)/255;%imrotate(,norm*360/8,'crop') (other variants)
dat = state;

function dat = cosColor(dat,ind,gain)
global accumInd
accumInd = accumInd + 1;% + ind*gain;
dat = (cos(2*pi*(double(dat)*(ind*gain)+accumInd ))+1)/2;