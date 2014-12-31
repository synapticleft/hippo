function xmas2014a(moviename, pixelFormat, maxThreads)
global accumInd
global state
global lastAudio
global transform

accumInd = 0;
gain = 1;
%VIDSIZE = [720 1280];% [240 320];
%state = zeros(VIDSIZE(1),VIDSIZE(2),3);
% PlayMoviesDemo(moviename [, backgroundMaskOut][, tolerance][, pixelFormat=4][, maxThreads=-1])
%
% This demo accepts a pattern for a valid moviename, e.g.,
% moviename=`*.mpg`, then it plays all movies in the current working
% directory whose names match the provided pattern, e.g., the `*.mpg`
% pattern would play all MPEG files in the current directory.
%
% This demo uses automatic asynchronous playback for synchronized playback
% of video and sound. Each movie plays until end, then rewinds and plays
% again from the start. Pressing the Cursor-Up/Down key pauses/unpauses the
% movie and increases/decreases playback rate.
% The left- right arrow keys jump in 1 seconds steps. SPACE jumps to the
% next movie in the list. ESC ends the demo.
if isempty(moviename)
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
% We retrieve status once to get access to SampleRate:
%s = PsychPortAudio('GetStatus', pahandle);
audiodata = 0;
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

try
    % Open onscreen window with gray background:
    screen = max(Screen('Screens'));
    scrsz = get(screen,'ScreenSize');
    win = PsychImaging('OpenWindow', screen, [0, 0, 0]);
    %PsychImaging('AddTask', 'General', 'UsePanelFitter',VIDSIZE);
     shader = [];
    % Use default pixelFormat if none specified:
    if nargin < 4
        pixelFormat = [];
    end
    
    % Use default maxThreads if none specified:
    if nargin < 5
        maxThreads = [];
    end
    
    % Initial display and sync to timestamp:
    Screen('Flip',win);
    iteration = 0;
    abortit = 0;
    
    % Use blocking wait for new frames by default:
    blocking = 1;
    
    % Default preload setting:
    preloadsecs = [];
    
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
        fprintf('ITER=%i::', iteration);
        
        % Open movie file and retrieve basic info about movie:
        [movie, movieduration, fps, imgw, imgh] = Screen('OpenMovie', win, moviename, [], preloadsecs, [], pixelFormat, maxThreads);
        state = zeros(imgh,imgw,3);
        fprintf('Movie: %s  : %f seconds duration, %f fps, w x h = %i x %i...\n', moviename, movieduration, fps, imgw, imgh);
        
        i=0;
        
        % Start playback of movie. This will start
        % the realtime playback clock and playback of audio tracks, if any.
        % Play 'movie', at a playbackrate = 1, with endless loop=1 and
        % 1.0 == 100% audio volume.
        Screen('PlayMovie', movie, rate, 1, 1.0);
        
        t1 = GetSecs;
        
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
                    lastAudio = audiodata;
                    audiodata = temp;
                end
                %nrsamples = size(audiodata, 2);
                
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
                %tic;
                if transform == 1
                    jitterback(movie,audiodata,gain,movieduration );
                else
                    I = Screen('GetImage',tex);
                    Screen('Close',tex);
                    if transform == 2
                        I = smearUpdate(I,audiodata,gain);
                    elseif transform == 3
                        I = cosColor(I,audiodata,gain);
                    elseif transform == 4
                        I = shiftss(I,audiodata,gain);
                    end
                    tex = Screen('MakeTexture',win,I);
                end
                Screen('DrawTexture', win, tex);%, [], scrsz);%, [], [], [], [], shader);
                Screen('DrawTexture', win, tex,[],scrsz);
                %toc;
                % Update display:
                Screen('Flip', win);
                % Release texture:
                Screen('Close', tex);
                % Framecounter:
                i=i+1;
            end;
            
            % Further keyboard checks...
            
            if (keyIsDown==1 && keyCode(right))
                % Advance movietime by one second:
                Screen('SetMovieTimeIndex', movie, Screen('GetMovieTimeIndex', movie) + 1);
            end;
            
            if (keyIsDown==1 && keyCode(left))
                % Rewind movietime by one second:
                %Screen('SetMovieTimeIndex', movie, Screen('GetMovieTimeIndex', movie) - 1);
                gain = 1;
            end;
            
            if (keyIsDown==1 && keyCode(up))
                gain = gain * 1.2;
            end;
            
            if (keyIsDown==1 && keyCode(down))
                gain = gain / 1.2;
            end;
        end;
        
        telapsed = GetSecs - t1;
        fprintf('Elapsed time %f seconds, for %i frames.\n', telapsed, i);
        
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
catch %#ok<*CTCH>
    % Error handling: Close all windows and movies, release all ressources.
    Screen('CloseAll');
    % Stop capture:
    PsychPortAudio('Stop', pahandle);
    % Close the audio device:
    PsychPortAudio('Close', pahandle);

end

function jitterback(movie,ind,gain,movieduration )
global lastAudio
Screen('SetMovieTimeIndex', movie, mod(Screen('GetMovieTimeIndex', movie) + gain*(ind-lastAudio),movieduration));


function dat = shiftss(dat,ind,gain)
sh = round(ind*gain*100);
dat(:,:,1) = circshift(dat(:,:,1),[sh 0 0]);
dat(:,:,2) = circshift(dat(:,:,2),[0 sh 0]);
dat(:,:,3) = circshift(dat(:,:,3),[-sh -sh 0]);

function dat = smearUpdate(dat,ind,gain)
global state
state = state*.9;%(amt+1)/2;
indbool = state < ind*gain;%*dat %narrower parameter space
state(indbool) = dat(indbool);
dat = state;

function dat = cosColor(dat,ind,gain)
global accumInd
%dat = (cos(2*pi*(double(dat)/256 + ind)) + 1)/2;
accumInd = accumInd + 1;% + ind*gain;
dat = (cos(2*pi*(double(dat)*(ind*gain)+accumInd ))+1)/2;%+ toc*time