function picsToMovie

cycle = 64;
%pathIn = ['I:\Downloaded Albums\pics\'];
%pathOut = ['C:\wiki\'];
file = 'Maria.avi';
%colormap gray;
counter = 0;
v.width = 640;v.height = 480;
    tic;
%    folders = getDirs(pwd);
wr = VideoWriter(file,'Motion JPEG AVI');
open(wr);

%for i = 1:numel(folders)
%    cd folders{i}
%    d = dir([path num2str(i) '_*.jpg']);
d = rdir(['**', '/*.*']);
d = {d.name}';
numel(d)
for j = 1:numel(d)
    if ismember(d{j}(end-2:end),{'JPG','jpg'})
        try
        im = imread(d{j});
        im = imresize(im,min(v.height/size(im,1),v.width/size(im,2)),'bicubic');
        im = im(1:(2*floor(size(im,1)/2)),1:(2*floor(size(im,2)/2)),:);
        if ndims(im) == 2
            im = double(im);
            im = uint8(im/max(im(:))*255);
            im = repmat(im,[1 1 3]);
        end
        im = padarray(im,[(v.height-size(im,1))/2 (v.width-size(im,2))/2 0],0);
        writeVideo(wr,im);
         catch
             d{j}
         end
        
%         counter = counter + 1;
%         v.frames(mod(counter-1,cycle)+1).cdata = im;
%         v.times(mod(counter-1,cycle)+1) = counter /10;
%         if ~mod(counter,cycle)
%             toc
%             if counter / cycle == 1
%                 mmwrite([pathOut file],v,'Continue');
%             elseif counter/cycle < numel(d)/cycle
%                 mmwrite([pathOut file],v,'Continue','Initialized');
%             else
%                 mmwrite([pathOut file],v,'Initialized');
%                 return
%             end
%             if ~mod(counter/cycle,10)
%                 counter/cycle
%             end
%             toc
%         end
    end
    if ~mod(j,100)
        j
    end
end
close(wr);

% function nameFolds = getDirs(fol)
% d = dir(fol);
% isub = [d(:).isdir]; %# returns logical vector
% nameFolds = {d(isub).name}';
% nameFolds(ismember(nameFolds,{'.','..'})) = [];