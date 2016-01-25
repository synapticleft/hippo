function halloween2015()
global accumInd
global state
global lastAudio
global transform
global turn

accumInd = 0;
lastAudio = 0;
brighter = 1;
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
KbName('UnifyKeyNames');
% Setup key mapping:
space=KbName('SPACE');
esc=KbName('ESCAPE');
right=KbName('RightArrow');
left=KbName('LeftArrow');
up=KbName('UpArrow');
down=KbName('DownArrow');
shift=KbName('RightShift');

imgh = 480;imgw = 640;

ts = 'qwerty';
for i = 1:numel(ts)
    tforms(i) = KbName(ts(i));
end
jj = KbName('j');
kk = KbName('k');
ll = KbName('l');

try
    % Open onscreen window with gray background:
    screen = max(Screen('Screens'));
    scrsz = get(screen,'ScreenSize');
    win = PsychImaging('OpenWindow', screen, [0, 0, 0]);
    % Initial display and sync to timestamp:
    Screen('Flip',win);
    grabber = Screen('OpenVideoCapture', win, [], [0 0 imgw imgh],[],[],[],[],[],[],8);
    Screen('StartVideoCapture', grabber, 30, 1);
    
    abortit = 0;
    % Playbackrate defaults to 1:
    rate=1;
    audiodata = 1;
    while (abortit<2)
        state = zeros(imgh,imgw,3);
        [dx dy] = meshgrid(1:imgw,1:imgh);
        dx = dx-mean(dx(:));dy = dy-mean(dy(:));
        turn = [];turn(:,:,1) = -dx/10;turn(:,:,2) = -dy/10;
        i=0;
        
        while 1
            % Check for abortion:
            abortit=0;
            [keyIsDown,secs,keyCode]=KbCheck;
            if (keyIsDown==1 && keyCode(esc))
                % Set the abort-demo flag.
                abortit=2;
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
                [tex pts nrdropped,intensity]=Screen('GetCapturedImage', win, grabber, 1);                 
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
%                if transform == 1
%                    jitterback(movie,audiodata,gain,movieduration );
%                else
                    I = Screen('GetImage',tex);
                    I = min(255,max(0,brighter*I));
                    Screen('Close',tex);
                    if transform == 1
                        I = slowEmerge(I,audiodata,gain);
                    elseif transform == 2
                        I = smearUpdate(I,audiodata,gain);
                    elseif transform == 3
                        I = cosColor(I,audiodata,gain);
                    elseif transform == 4
                        I = shiftss(I,audiodata,gain);
                    elseif transform == 5
                        I = tunnel(I,audiodata,gain);
                    end
                    tex = Screen('MakeTexture',win,I);
%                end
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
%                Screen('SetMovieTimeIndex', movie, Screen('GetMovieTimeIndex', movie) + 1);
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
            if (keyIsDown==1 && keyCode(jj))
                % Rewind movietime by one second:
                %Screen('SetMovieTimeIndex', movie, Screen('GetMovieTimeIndex', movie) - 1);
                brighter = 1;
            end;
            
            if (keyIsDown==1 && keyCode(kk))
                brighter = brighter / 1.2;
            end;
            
            if (keyIsDown==1 && keyCode(ll))
                brighter = brighter * 1.2;
            end;
        end;
        
        Screen('Flip', win);
        KbReleaseWait;
    end;
    
    % Close screens.
    Screen('CloseAll');
    Screen('StopVideoCapture', grabber);
    Screen('CloseVideoCapture', grabber);
    return;
catch %#ok<*CTCH>
    % Error handling: Close all windows and movies, release all ressources.
    Screen('CloseAll');
end

function jitterback(movie,ind,gain,movieduration )
global lastAudio
Screen('SetMovieTimeIndex', movie, mod(Screen('GetMovieTimeIndex', movie) + gain*(ind-lastAudio),movieduration));%-lastAudio


function dat = shiftss(dat,ind,gain)
%sh = round(ind*gain*100);
[x, y] = gradient(double(rgb2gray(dat)));
%x = imfilter(x,fspecial('gaussian',10,10));y = imfilter(y,fspecial('gaussian',10,10));
x(:,:,2) = y;
dat = imwarp(dat,x*ind*gain);
clear x y;
%dat(:,:,1) = circshift(dat(:,:,1),[sh 0 0]);
%dat(:,:,2) = circshift(dat(:,:,2),[0 sh 0]);
%dat(:,:,3) = circshift(dat(:,:,3),[-sh -sh 0]);

function dat = slowEmerge(dat,ind,gain)
global state
state = state*(1-ind*gain*.1) + ind*gain*.1*double(dat)/255;
dat = min(1,max(0,state));

function dat = smearUpdate(dat,ind,gain)
global state
state = state*.9;%(amt+1)/2;
indbool = state < ind*gain;%*dat %narrower parameter space
state(indbool) = dat(indbool);
dat = state;

function dat = tunnel(dat,ind,gain)
global state
global turn
norm = min(1,ind*gain);
state = max(.9*imwarp(state,norm*turn),double(dat)/255);
%state = (norm)*imwarp(state,turn) + (1-norm)*double(dat)/255;%imrotate(,norm*360/8,'crop')
dat = state;

function dat = cosColor(dat,ind,gain)
global accumInd
%dat = (cos(2*pi*(double(dat)/256 + ind)) + 1)/2;
accumInd = accumInd + 1/100;% + ind*gain;
dat = (cos(2*pi*(double(dat)*(ind*gain)+accumInd ))+1)/2;
%dat = rgb2hsv(double(dat)/256);
%dat(:,:,1) = (cos(2*pi*(double(dat(:,:,1))*(ind*gain)+accumInd ))+1)/2;%+ toc*time
%dat = hsv2rgb(dat);