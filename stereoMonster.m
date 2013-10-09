function stereoMonster()
% ImagingStereoDemo([stereoMode=8][, usedatapixx = 0][, writeMovie = 0])
%
% Demo on how to use OpenGL-Psychtoolbox to present stereoscopic stimuli
% when the Psychtoolbox imaging pipeline is enabled. Use of the imaging
% pipeline allows for more flexible and high quality stereo display modes,
% but it requires graphics hardware with support for at least framebuffer
% objects and Shadermodel 2.0. See the Psychtoolbox Wiki about gfx-hardware
% recommendations. The demo also shows how to configure the pipeline to
% restrict image processing to some subregion of the display, e.g., to save
% some computation time on low-end hardware.
%
% Press escape key to abort demo, space key to toggle modes of specific
% algorithms.
%
% 'stereoMode' specifies the type of stereo display algorithm to use:
%
% 0 == Mono display - No stereo at all.

%
% 4 == Free fusion (lefteye=left, righteye=right): Left-eye view displayed
% in left half of window, right-eye view displayed in right-half of window.
% Use this for dual-display setups (binocular video goggles, haploscopes,
% polarized stereo setups etc.)

% 'writeMovie' If provided and set to a non-zero value will write a
% Quicktime movie file 'MyTestMovie.mov' into the current working directory
% which captures the full performance of this demo. A setting of 1 will
% only write video, a setting of 2 will also write an audio track with
% a sequence of ten successive beep tones of 1 sec duration.
%

rootDir = '/home/redwoo/Desktop/gsrfilms';%Videos';%;
moviefiles=dir([rootDir filesep '*.avi']);%avi']);
for i=1:size(moviefiles,1)
    moviefiles(i).name = [ rootDir filesep moviefiles(i).name ];
end
moviecount = size(moviefiles,1);

% We start of with non-inverted display:
inverted = 0;
    stereoMode = 4;
    
% This script calls Psychtoolbox commands available only in OpenGL-based
% versions of the Psychtoolbox. (So far, the OS X Psychtoolbox is the
% only OpenGL-base Psychtoolbox.)  The Psychtoolbox command AssertPsychOpenGL will issue
% an error message if someone tries to execute this script on a computer without
% an OpenGL Psychtoolbox
AssertOpenGL;

% Define response key mappings, unify the names of keys across operating
% systems:
KbName('UnifyKeyNames');
space = KbName('space');
escape = KbName('ESCAPE');
right=KbName('RightArrow');
left=KbName('LeftArrow');
up=KbName('UpArrow');
down=KbName('DownArrow');
shift=KbName('RightShift');

%try
% Get the list of Screens and choose the one with the highest screen number.
% Screen 0 is, by definition, the display with the menu bar. Often when
% two monitors are connected the one without the menu bar is used as
% the stimulus display.  Chosing the display with the highest dislay number is
% a best guess about where you want the stimulus displayed.
scrnNum = max(Screen('Screens'));

% Open double-buffered onscreen window with the requested stereo mode,
% setup imaging pipeline for additional on-the-fly processing:

% Prepare pipeline for configuration. This marks the start of a list of
% requirements/tasks to be met/executed in the pipeline:
PsychImaging('PrepareConfiguration');

% Ask to restrict stimulus processing to some subarea (ROI) of the
% display. This will only generate the stimulus in the selected ROI and
% display the background color in all remaining areas, thereby saving
% some computation time for pixel processing: We select the center
% 512x512 pixel area of the screen:
%if ~ismember(stereoMode, [100, 101])
%    PsychImaging('AddTask', 'AllViews', 'RestrictProcessing', CenterRect([0 0 512 512], Screen('Rect', scrnNum)));
%end

% Consolidate the list of requirements (error checking etc.), open a
% suitable onscreen window and configure the imaging pipeline for that
% window according to our specs. The syntax is the same as for
% Screen('OpenWindow'):
[windowPtr, windowRect]=PsychImaging('OpenWindow', scrnNum, 0, [], [], [], stereoMode);

% Initially fill left- and right-eye image buffer with black background
% color:
%Screen('SelectStereoDrawBuffer', windowPtr, 0);
%Screen('FillRect', windowPtr, BlackIndex(scrnNum));
%Screen('SelectStereoDrawBuffer', windowPtr, 1);
%Screen('FillRect', windowPtr, BlackIndex(scrnNum));

% Show cleared start screen:
%Screen('Flip', windowPtr);

% Set up alpha-blending for smooth (anti-aliased) drawing of dots:
Screen('BlendFunction', windowPtr, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%Screen('Flip', windowPtr);

% Maximum number of animation frames to show:
%nmax = 100000;

% Perform a flip to sync us to vbl and take start-timestamp in t:
%t = Screen('Flip', windowPtr);
abortit = 0;iteration = 0;preloadsecs = [];blocking = 1;shader = [];
hz = 0;mode = 1;scale = 1;ctr = 0;
% Run until a key is pressed:
while (abortit<2)%length(t) < nmax
    iteration=iteration + 1;
    moviename=moviefiles(mod(iteration, moviecount)+1).name;
    
    % Open movie file and retrieve basic info about movie:
    [movie movieduration fps imgw imgh] = Screen('OpenMovie', windowPtr, moviename, [], preloadsecs);
    fprintf('Movie: %s  : %f seconds duration, %f fps, w x h = %i x %i...\n', moviename, movieduration, fps, imgw, imgh);
    allImages = zeros(5,imgh,imgw,3);
    % Start playback of movie. This will start
    % the realtime playback clock and playback of audio tracks, if any.
    % Play 'movie', at a playbackrate = 1, with endless loop=1 and
    % 1.0 == 100% audio volume.
    Screen('PlayMovie', movie, 1, 1, 1.0);
    [x y] = meshgrid((1:imgw)-1,(1:imgh)-1);
    x = x-imgw/2;y = y-imgh/2;
    x = exp(1i*angle(fftshift(x+1i*y))*2);
    x = exp(1i*angle(tanh(real(x)*10) + 1i*tanh(imag(x)*10)));
    while 1
        abortit=0;
        [keyIsDown,~,keyCode]=KbCheck;
        KbReleaseWait;
        tex = Screen('GetMovieImage', windowPtr, movie, blocking);
         % Valid texture returned?
        if tex < 0
            break;
        end
        if tex == 0
            % No new frame in polling wait (blocking == 0). Just sleep
            % a bit and then retry.
            WaitSecs('YieldSecs', 0.005);
              continue;
        end
        myImage = Screen('GetImage', tex);
         if mode == 0
            texL = Screen('MakeTexture',windowPtr,myImage*((1+cos(GetSecs*hz))/2));
            texR = Screen('MakeTexture',windowPtr,myImage*((1-cos(GetSecs*hz))/2));
         elseif mode == 1
             myImage = 256*im2double(myImage);
            y = (1+real(exp(1i*GetSecs*hz)*x))/2;
            y(1,1) = .5;
            y = repmat(y,[1 1 3]);
              y = max(0,min(255,real(ifft2(fft2(myImage).*y))));
             texL = Screen('MakeTexture',windowPtr,y);
             y = max(0,min(255,myImage-y));
             texR = Screen('MakeTexture',windowPtr,y);%real(ifft2(fftshift(myImage.*(1-y))))
        elseif mode == 2
            myImage = 256*im2double(myImage);
            %im1 = 256*(1+cos(double(myImage*2*pi + GetSecs*hz))/2);
            y = mod(myImage*scale+GetSecs*hz,1);
            y = abs(2*y-1);
            texL = Screen('MakeTexture',windowPtr,myImage.*y);
            texR = Screen('MakeTexture',windowPtr,myImage.*(1-y));
         elseif mode == 3
            texL = Screen('MakeTexture',windowPtr,myImage);
            texR = Screen('MakeTexture',windowPtr,squeeze(allImages(ctr+1,:,:,:)));
         elseif mode == 4
            myImage = 256*im2double(myImage);
            y = repmat((1+cos((1:imgh)*scale*2*pi/imgh + GetSecs*hz)')/2,[1 imgw 3]);
            texL = Screen('MakeTexture',windowPtr,myImage.*y);
            texR = Screen('MakeTexture',windowPtr,myImage.*(1-y));
        elseif mode == 5
            myImage = 256*im2double(myImage);
            y = exp(1i*((1:imgh)*scale*2*pi/imgh + GetSecs*hz)');
            y = repmat(y,[1 imgw 3]);
            y(:,:,2) = y(:,:,2)*exp(1i*2*pi/3);
            y(:,:,3) = y(:,:,3)*exp(1i*4*pi/3);
            y = (1+real(y))/2;
            texL = Screen('MakeTexture',windowPtr,myImage.*y);
            texR = Screen('MakeTexture',windowPtr,myImage.*(1-y));
        end
        %allImages = circshift(allImages,[1 0 0 0]);
        %allImages(1,:,:,:) = myImage;
        allImages(ctr+1,:,:,:) = myImage;%mod(ctr-round(scale)+1,size(allImages,1))+1,:,:,:) = myImage;%ctr+1
        ctr = mod(ctr + 1,size(allImages,1));
        % Select left-eye image buffer for drawing:
        Screen('SelectStereoDrawBuffer', windowPtr, 0);
        Screen('DrawTexture', windowPtr, texL, [], windowRect, [], [], [], [], shader);
        % Select right-eye image buffer for drawing:
        Screen('SelectStereoDrawBuffer', windowPtr, 1);
        Screen('DrawTexture', windowPtr, texR, [], windowRect, [], [], [], [], shader);
        % Tell PTB drawing is finished for this frame:
        Screen('DrawingFinished', windowPtr);
        
        % Flip stim to display and take timestamp of stimulus-onset after
        % displaying the new stimulus and record it in vector t:
        Screen('Flip', windowPtr);
        Screen('Close',tex);
        Screen('Close',texL);
        Screen('Close',texR);
        if (keyIsDown==1 && keyCode(escape))
            % Set the abort-demo flag.
            abortit=2;
            break;
        elseif (keyIsDown==1 && keyCode(space))
            % Exit while-loop: This will load the next movie...
            break;
        elseif (keyIsDown==1 && keyCode(up))
            hz = min(20*2*pi,max(hz*1.5,.1));
        elseif (keyIsDown==1 && keyCode(down))
            hz = max(.1,hz/1.5)-.1;
        elseif (keyIsDown==1 && keyCode(right))
            scale = min(30,scale*2);
        elseif (keyIsDown==1 && keyCode(left))
            scale = max(1,scale/2);
        elseif (keyIsDown==1 && keyCode(shift))
            mode = mod(mode+1,6);
        end
    end
    Screen('Flip', windowPtr);
    % Done. Stop playback:
    Screen('PlayMovie', movie, 0);
    % Close movie object:
    Screen('CloseMovie', movie);
end

% Done. Close the onscreen window:
Screen('CloseAll')