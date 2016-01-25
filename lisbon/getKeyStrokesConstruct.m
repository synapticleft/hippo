function fullData = getKeyStrokesConstruct %[dataOut,pos]
%reads all of the files in a directory and organizes the data into a struct
%according to the identity of the player
% 'Pos' stores all vertices visited and their times
% 'Strokes' stores all keystrokes
% 'Lose' stores all losses
% 'Pop' stores all bubble loss events
% 'Date' stores identifying date (first play logged by browser)
% 'Name' stores player name
% 'CPU' stores device type
% 'gametype' stores number of vertices

%Time
%event
%wave
%score
%bank
%angle
%graphicScore
%graphicBank
%waveTime
%L = lose, V = vertex, B = bubble
num = [-40,-39,-38,-37,0,1,37,38,39,40];
files = dir('*.txt');
%cols = ['b','r'];
fullData.users = {};
for k = 1:numel(files)
    data = loadjson(files(k).name);%file1);
    data = data.data;
    eq = [];
    for l = 1:numel(fullData.users)
        eq(l) = streq(data{1}{1},fullData.users{l}{1}) && streq(data{1}{3},fullData.users{l}{2});
    end
    if ~sum(eq)%~ismember(data{1,1}{1},names)
        l = l+1;
        fullData.users{l}.session = [];
    else
        l = find(eq);
    end
    fullData.users{l}.name = data{1}{1};
    fullData.users{l}.system = data{1}{3};
    fullData.users{l}.session{end+1}.time = data{1}{2};
    fullData.users{l}.session{end}.numVertices = data{1}{4};
    fullData.users{l}.session{end}.keyPress = [];
    fullData.users{l}.session{end}.keyRelease = [];
    fullData.users{l}.session{end}.keyposition = [];
    fullData.users{l}.session{end}.lose = [];
    fullData.users{l}.session{end}.pop = [];
    %startFile1 = data{1}{2}/86400/1000 + datenum(1970,1,1);
    
    for i = 3:length(data)-1
        if ~iscell(data{i})
            if data{i}(2) > 0
                fullData.users{l}.session{end}.keyPress(end+1) = data{i};
            else
                fullData.users{l}.session{end}.keyRelease(end+1) = data{i};
            end
            %        dataOut(i-2,:) = data{i};
        else
            if data{i}{2}{1} == 'V'
                fullData.users{l}.session{end}.position(end+1) = data{i}([1 .1 3:end]);
            elseif data{i}{2}{1} == 'L'
                fullData.users{l}.session{end}.lose(end+1) = data{i}([1 .2 3:end]);
            elseif data{i}{2}{1} == 'B'
                fullData.users{l}.session{end}.pop(end+1) = data{i}([1 .3 3:end]);
            end
            %         for j = 1:9
            %             if iscell(data{i}{j})
            %                 if data{i}{j}{1} == 'V'%data{i}{j}{1}(1) == 'R'
            %                     dataOut(i-2,j) = -1;%dataOut(i-2,j) = -1*str2num(data{i}{j}{1}(2:end));
            %                 elseif data{i}{j}{1} == 'L'
            %                     dataOut(i-2,j) = 0;
            %                 elseif data{i}{j}{1} == 'B'
            %                     dataOut(i-2,j) = 1;
            %                 else
            %                     dataOut(i-2,j) = str2num(data{i}{j}{1});
            %                 end
            %             else
            %                 dataOut(i-2,j) = data{i}{j};
            %             end
            %         end
        end
    end
    fullData.users{l}.session{end}.position(:,6) = mod(fullData.users{l}.session{end}.position(:,6),360);
    fullData.users{l}.session{end}.keyPress(~ismember(fullData.users{l}.session{end}.keyPress(:,2),num)) = [];
dataOut(dataOut(:,9) == 0,:) = [];
%dataOut(dataOut(:,6) == 360,6) = 0;
%pos = dataOut(dataOut(:,2) == -1,:);
%dataOut(dataOut(:,2) == -1,:) = [];
%dataOut(~ismember(dataOut(:,2),num),:) = [];
%dataOut(dataOut(:,3) == 0,3) = 8;

% for l = 1:numel(names)
%     eq(l) = streq(data{1,1}{1},names{l});
% end
% figure(1);hold all;
% brk = find(abs(diff(pos(:,6))) > 100 | diff(pos(:,9) < 0));
% brk = [0 brk' size(pos,1)];
% posa = pos(:,9)  + (pos(:,3)-1)*40;%pos(:,1);
% for m = 1:numel(brk)-1
%     plot(posa((brk(m)+1):brk(m+1)),pos((brk(m)+1):brk(m+1),6),cols(find(eq)));
% end
% temp = dataOut(:,9)+(dataOut(:,3)-1)*40;%dataOut(:,1);
% drawnow;
% figure(2);hold all;plot(temp,dataOut(:,4)+dataOut(:,5),cols(find(eq)));drawnow;
% %figure(49);hold all;plot(dataOut(:,9));
% %figure(50);hold all;plot(dataOut(:,9)+40*(dataOut(:,3)-1));
end