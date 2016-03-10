function fullData = getKeyStrokesConstruct %[dataOut,pos]
%% reads all of the files in a directory and organizes the data into a struct
%% with a separate field for each player (ASSUMES EACH PLAYER ENTERS A DIFFERENT NAME
% 'Name' stores player name
% 'System' stores device type (eg. windows, linux, osx)
% 'session' stores info for each game played by player, with fields:
%   'Time' stores when game was played
%   'numVertices' stores number of vertices
%   'event' stores all events in a game, each with 9 entries:
%       Time
%       event type (key presses and other game events) *See below for how each event is represented*
%       wave number
%       current score
%       amount in bank
%       angle of spinner
%       size of polygon corresponding to score
%       size of polygon corresponding to bank
%       time within the current wave

% press (left, up, right, down); release (left up right down);land on vertex; lose life; lose bubble   
num = [-40,-39,-38,-37,37,38,39,40,.1,.2,.3];%;0,1,
num= [-18, -16, -15, -14, -1, 1, 14, 15, 16, 18,.1,.2,.3];
files = dir('*.txt');
fullData.users = {};
for k = 1:numel(files)
    data = loadjson(files(k).name);
    data = data.data;
    eq = [];
    for l = 1:numel(fullData.users)
        eq(l) = strcmp(data{1}{1},fullData.users{l}.name) && strcmp(data{1}{3},fullData.users{l}.system);
    end
    if ~sum(eq)
        if isempty(l)
            l = 1;
        else
            l = l+1;
        end
        fullData.users{l}.session = [];
    else
        l = find(eq);
    end
    fullData.users{l}.name = data{1}{1};
    fullData.users{l}.system = data{1}{3};
    fullData.users{l}.session{end+1}.time = data{1}{2};
    fullData.users{l}.session{end}.numVertices = data{1}{4};
    fullData.users{l}.session{end}.event = [];
    for i = 3:length(data)-1
        if ~iscell(data{i})
            fullData.users{l}.session{end}.event(end+1,:) = data{i};
        else
            if data{i}{2}{1} == 'V'
                fullData.users{l}.session{end}.event(end+1,:) = [cell2mat(data{i}(1)) .1 cell2mat(data{i}(3:end))];%vertex land
            elseif data{i}{2}{1} == 'L'
                fullData.users{l}.session{end}.event(end+1,:) = [cell2mat(data{i}(1)) .2  cell2mat(data{i}(3:end))];%lose game
            elseif data{i}{2}{1} == 'B'
                fullData.users{l}.session{end}.event(end+1,:) = [cell2mat(data{i}(1)) .3 cell2mat(data{i}(3:end))];%lose bubble
            end
        end
    end
    fullData.users{l}.session{end}.event(:,6) = mod(fullData.users{l}.session{end}.event(:,6),360);
    fullData.users{l}.session{end}.event(~ismember(fullData.users{l}.session{end}.event(:,2),num),:) = [];
    fullData.users{l}.session{end}.event(fullData.users{l}.session{end}.event(:,9) == 0,:) = [];
end