function xmas2014a(moviename)
global accumInd
global state
global lastAudio
global transform

accumInd = 0;
gain = 1;

% Initialize with unified keynames and normalized colorspace:
PsychDefaultSetup(2);
InitializePsychSound;

% Open the default audio device [], with mode 2 (== Only audio capture),
% and a required latencyclass of zero 0 == no low-latency mode, as well as
% a frequency of 44100 Hz and 2 sound channels for stereo capture.
% This returns a handle to the audio device:
freq = 44100;
%pahandle = PsychPortAudio('Open', [], 2, 0, freq, 1);
%PsychPortAudio('GetAudioData', pahandle, 10);
%PsychPortAudio('Start', pahandle, 0, 0, 1);
%s = PsychPortAudio('GetStatus', pahandle);
audiodata = 0;

%try
    % Open onscreen window with gray background:
    screen = max(Screen('Screens'));
    scrsz = get(screen,'ScreenSize');
    win = PsychImaging('OpenWindow', screen, [0, 0, 0],[0 0 1280 720]);
     shader = [];
        pixelFormat = [];
        maxThreads = [];
    Screen('Flip',win);
    iteration = 0;
    abortit = 0;
    
    % Default preload setting:
    preloadsecs = [];
    
    % Return full list of movie files from directory+pattern:
    moviefiles=dir(moviename);
        moviename = [ pwd filesep moviefiles(1).name ];
    rate=1;
        % Open movie file and retrieve basic info about movie:
        [movie, movieduration, fps, imgw, imgh] = Screen('OpenMovie', win, moviename, [], preloadsecs, [], pixelFormat, maxThreads);
        state = zeros(imgh,imgw,3);
        fprintf('Movie: %s  : %f seconds duration, %f fps, w x h = %i x %i...\n', moviename, movieduration, fps, imgw, imgh);
        
        i=0;
        
        % Start playback of movie. This will start
        % the realtime playback clock and playback of audio tracks, if any.
        % Play 'movie', at a playbackrate = 1, with endless loop=1 and
        % 1.0 == 100% audio volume.
        Screen('PlayMovie', movie, rate, 1, 0);
        % Infinite playback loop: Fetch video frames and display them...
        tic
        while 1
            toc
                % Retrieve pending audio data from the drivers internal ringbuffer:
%                temp = mean(abs(PsychPortAudio('GetAudioData', pahandle)),2);
%                if ~isempty(temp)
%                    lastAudio = audiodata;
%                    audiodata = temp;
%                end
                %nrsamples = size(audiodata, 2);
                
                % Return next frame in movie, in sync with current playback
                % time and sound.
                % tex is either the positive texture handle or zero if no
                % new frame is ready yet in non-blocking mode (blocking == 0).
                % It is -1 if something went wrong and playback needs to be stopped:
                tex = Screen('GetMovieImage', win, movie, 1)
                if tex < 0
                    break;
                end
%                 if tex == 0
%                     % No new frame in polling wait (blocking == 0). Just sleep
%                     % a bit and then retry.
%                     WaitSecs('YieldSecs', 0.005);
%                     continue;
%                 end
                % Draw the new texture immediately to screen:
                %tic;
%                 if transform == 1
%                     jitterback(movie,audiodata,gain,movieduration );
%                 else
                     I = Screen('GetImage',tex);
                    Screen('Close',tex);
%                        I = smearUpdate(I,audiodata,gain);
                        I = cosColor(I,audiodata,gain);
                    tex = Screen('MakeTexture',win,I);
                Screen('DrawTexture', win, tex);%, [], scrsz);%, [], [], [], [], shader);
                Screen('DrawTexture', win, tex,[],scrsz);
                %toc;
                % Update display:
                Screen('Flip', win);
                % Release texture:
                Screen('Close', tex);
                % Framecounter:
                i=i+1;
        end
        Screen('Flip', win);
        Screen('PlayMovie', movie, 0);
        Screen('CloseMovie', movie);
    Screen('CloseAll');
%     PsychPortAudio('Stop', pahandle);
%     PsychPortAudio('Close', pahandle);
% catch %#ok<*CTCH>
%     % Error handling: Close all windows and movies, release all ressources.
%     Screen('CloseAll');
%     % Stop capture:
% %    PsychPortAudio('Stop', pahandle);
%     % Close the audio device:
%     PsychPortAudio('Close', pahandle);
% end

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