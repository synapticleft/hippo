function rhythMovie(moviename)
global accumInd
global state
global gain
global ranger
global xs
global ys

accumInd = 0;
gain = .2;
ranger = 30^2;

% Initialize with unified keynames and normalized colorspace:
PsychDefaultSetup(2);
InitializePsychSound;
freq = 44100;
audiodata = 0;
space=KbName('SPACE');
esc=KbName('ESCAPE');
right=KbName('RightArrow');
left=KbName('LeftArrow');
up=KbName('UpArrow');
down=KbName('DownArrow');
shift=KbName('RightShift');

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
    %moviefiles=dir(moviename);
        %moviename = [ pwd filesep moviefiles(1).name ];
    rate=1;
        % Open movie file and retrieve basic info about movie:
        [movie, movieduration, fps, imgw, imgh] = Screen('OpenMovie', win, moviename, [], preloadsecs, [], pixelFormat, maxThreads);
        state = zeros(imgh,imgw,3);
        fprintf('Movie: %s  : %f seconds duration, %f fps, w x h = %i x %i...\n', moviename, movieduration, fps, imgw, imgh);
        
        i=0;
        [xs, ys] = meshgrid(1:imgw,1:imgh);
        % Start playback of movie. This will start
        % the realtime playback clock and playback of audio tracks, if any.
        % Play 'movie', at a playbackrate = 1, with endless loop=1 and
        % 1.0 == 100% audio volume.
        Screen('PlayMovie', movie, rate, 1, 1);
        % Infinite playback loop: Fetch video frames and display them...
        tic
        while 1
            toc
                tex = Screen('GetMovieImage', win, movie, 1)
                [keyIsDown,~,keyCode]=KbCheck;
                [x,y] = GetMouse;
                if tex < 0
                    break;
                end
%                 if transform == 1
%                     jitterback(movie,audiodata,gain,movieduration );
%                 else
                     I = Screen('GetImage',tex);
                     I = persist(I);
                     I = spotLight(I,x,y);
                    Screen('Close',tex);
%                        I = smearUpdate(I,audiodata,gain);
%                         I = cosColor(I,audiodata,gain);
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
            if (keyIsDown==1 && keyCode(up))
                gain = gain * 1.2;
            end;
            
            if (keyIsDown==1 && keyCode(down))
                gain = gain / 1.2;
            end;
            if (keyIsDown==1 && keyCode(left))
                ranger = ranger * 1.2;
            end;
            
            if (keyIsDown==1 && keyCode(right))
                ranger = ranger / 1.2;
            end;
        end
        Screen('Flip', win);
        Screen('PlayMovie', movie, 0);
        Screen('CloseMovie', movie);
    Screen('CloseAll');
    
    function dat = spotLight(dat,x,y)
        global ranger
        global xs
        global ys
        x = x - 100;
        y = y - 100;
        inds = exp(-((xs-x).^2 + (ys - y).^2)/ranger);
        dat = dat.*double(repmat(inds,[1 1 3]));
        %temp = zeros(size(dat));
        %temp(x,y) = 1;
        %temp = imfilter(temp,fspecial('gaussian',50,30));

    function dat = persist(dat)
        global state
        global gain
        dat = dat > 128;
        %state = state*gain + double(dat)*(1-gain);
        state = max(circshift(state*(1-exp(-gain)),[1 -3]),dat);%
        dat = state;
        
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