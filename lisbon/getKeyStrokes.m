function [dataOut,bonsDown] = getKeyStrokes(file1,file2)

isReversed = 1;

data = loadjson(file1);
data = data.data;
num = [-40,-39,-38,-37,37,38,39,40];
startFile1 = data{1}{2}/86400/1000 + datenum(1970,1,1);

for i = 3:length(data)-1
    if ~iscell(data{i})
        dataOut(i-2,:) = data{i};
    else
        for j = 1:6
            if iscell(data{i}{j})
                if data{i}{j}{1}(1) == 'R'
                    dataOut(i-2,j) = -1*str2num(data{i}{j}{1}(2:end));
                elseif data{i}{j}{1} == 'L'
                    dataOut(i-2,j) = 0;
                else
                    dataOut(i-2,j) = str2num(data{i}{j}{1});
                end
            else
                dataOut(i-2,j) = data{i}{j};
            end
        end
    end
end
dataOut(~ismember(dataOut(:,2),num),:) = [];
dataOut(dataOut(:,3) == 0,3) = 8;

pre = 'keydown_time';%up_time
fid = fopen([pre file2]);
temp = textscan(fid,'%s %f-%f-%fT%f:%f:%f+%f:00');
bonsUp = getTime(temp);
upKeys = temp{1};
fclose(fid);

pre = 'keyup_time';%down_time
fid = fopen([pre file2]);
temp = textscan(fid,'%s %f-%f-%fT%f:%f:%f+%f:00');
if isdst(datetime(temp{2}(1),temp{3}(1),temp{4}(1),'TimeZone','Europe/London'))
    startFile1 = startFile1 + 1/24;
end
bonsDown = getTime(temp);
downKeys = temp{1};
fclose(fid);

if isReversed
    tempt = bonsUp;
    tempk = upKeys;
    upKeys = downKeys;
    bonsUp = bonsDown;
    bonsDown = tempt;
    downKeys = tempk;
end

[bonsDown,downKeys] = fixPress(bonsUp,upKeys,bonsDown,downKeys);
for i = 1:numel(upKeys)
    upKeys{i} = ['R' upKeys{i}];
end
%[bonsUp,upKeys] = fixPress(bonsUp+.000001,upKeys,bonsUp,upKeys);

bonsTimes = [bonsUp; bonsDown'];
bonsKeys = [upKeys;downKeys'];
[bonsTimes,ind] = sort(bonsTimes);
bonsKeys = bonsKeys(ind);
%figure;plot(bonsDown(2:end)-dataOut(dataOut(:,2) < 0,1))
%dataOut(:,1) = dataOut(:,1) -bonsDown(1);
%% the following is to get game keylog w.r.t. initiation of data collectiondataOut(:,1) = dataOut(:,1) + getDiff(downKeys,bonsDown)*24*3600;
[bonsDown,keys] = getKey(bonsKeys,bonsTimes,dataOut(:,2),startFile1);

dataOut(ismember(dataOut(:,2),num),2) = keys(:,1);
dataOut(:,1) = dataOut(:,1)-bonsDown(1)*24*3600; %% THIS IS NOT VERIFIED TO BE THE RIGHT APPROACH
bonsDown = 24*3600*(bonsDown(2:end) - bonsDown(1));
bonsDown(:,2) = keys(:,2);
%bonsUp = 24*3600*getKey(downKeys,bonsDown,dataOut(:,2),startFile1,1);
%% the following is to find the accuracy of bonsai key-logging vs. game
%bonsDown = bonsDown(2:end) - bonsDown(1);

function [downTimeN,downKeyN] = fixPress(upTime,upKey,downTime,downKey)
keys = {'D1','Left','Right','Down','Up'};
downTimeN = [];
downKeyN = {};
for i = 1:numel(keys)
    f = find(streq(keys{i},upKey));
    f1 = find(streq(keys{i},downKey));
    downKeyTime = downTime(f1);
    %pressTimes = [upTime(f(2:end)); Inf];
    for j = 1:numel(f)%2:numel(f)
        downTimeN(end+1) = max(downKeyTime(downKeyTime < upTime(f(j))));%pressTimes(j)));
        downKeyN{end+1} = keys{i};
    end
end

function out = getTime(struc)
out = datenum([struc{2} struc{3} struc{4} struc{5} struc{6} struc{7}]);
%out = struc{5}*3600+struc{6}*60+struc{7};

function dif = getDiff(bk,bt)
for i = 1:length(bk)
    if streq('D1',bk(i))
        break
    end
end
dif = bt(i) - bt(1);

function [bt,k] = getKey(bk,bt,gk,gameTimeRef)
%isUp = 1-2*isUp;
gk(gk == 0) = [];
%gk = gk*-1;

let = {'RDown','RRight','RUp','RLeft','Left','Up','Right','Down'};
num = [-40,-39,-38,-37,37,38,39,40];
for i = 1:length(bk)
    if ismember(bk(i+1),let)%streq('RD1',bk(i))
        break
    end
end
gk1 = [];bk1 = [];
for j = 1:length(gk)
    if numel(find(gk(j) == num,1))
        gk1(end+1) = find(gk(j) == num,1);
    end
    if numel(find(streq(bk(i+j),let),1))
        bk1(end+1) = find(streq(bk(i+j),let),1);
    end
end
k = [gk1' bk1'];
%     if ~((gk(j) == 37 && streq(bk(i+j),'Left')) || (gk(j) == 39 && streq(bk(i+j),'Right')) ...
%             || (gk(j) == 38 && streq(bk(i+j),'Up')) || (gk(j) == 40 && streq(bk(i+j),'Down')) ...
%             || (gk(j) == -37 && streq(bk(i+j),'RLeft')) || (gk(j) == -39 && streq(bk(i+j),'RRight')) ...
%             || (gk(j) == -38 && streq(bk(i+j),'RUp')) || (gk(j) == -40 && streq(bk(i+j),'RDown')))
%         'Error!!! key identities from bonsai and game dont match!!!'
%         return
%     end
% end
bt = bt - gameTimeRef;%- bt(i);
bt(2:i) = [];
bt(j+2:end) = [];
