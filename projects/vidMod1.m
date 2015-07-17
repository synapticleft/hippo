function vidMod1(fn,sndFile)
%Take in video and transform it using matlab's builtin functions
[y,fs] = audioread(sndFile);
tic;
ind = 0;
videoFReader = vision.VideoFileReader(fn);%,'AudioOutputPort',1
videoFWriter = vision.VideoFileWriter('test2.avi','FrameRate',videoFReader.info.VideoFrameRate,'AudioInputPort',true);
%%videoPlayer = vision.VideoPlayer;
gain = 1;
ind = 0;
frame = 0;yState = 0;
samplingRatio = round(fs/videoFReader.info.VideoFrameRate);
while ~isDone(videoFReader)
  [im] = step(videoFReader);
  im = cosColor(im,exp(max(0,ind-10)/8)-1,gain);
  image(im);drawnow;
  ind = ind + .05;
  echo = min(ind/40,1);
  ySamp = y(frame+(1:samplingRatio),1);
  yState = (1-echo)*ySamp + echo*yState;
  %yf = fft(ySamp);
  %yf = yf*(cos(((1:length(yf))-length(yf)/2)/100)/2+1);
  step(videoFWriter,im,yState);
  frame = frame + samplingRatio;
%  step(videoPlayer, videoFrame);
end
%release(videoPlayer);
release(videoFWriter);
release(videoFReader);


function dat = cosColor(dat,ind,gain)
%dat = (cos(2*pi*(double(dat)/256 + ind)) + 1)/2;
%accumInd = accumInd + 1;% + ind*gain;
%dat = (cos(2*pi*(double(dat)*(ind*gain)+tim ))+1)/2;%+ toc*time
dat = (sin((dat*(1+ind)-1/2+ind)*pi)+1)/2;