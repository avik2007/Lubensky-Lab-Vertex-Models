function [cloneMatx, NORMcloneMatx, BoxSize] = ReadCloneInfoAppend(fid, startPos, clonetotal, negclone, foundersize, totalsteps, n)
% Sister function of CloneStats written to read specific file to spit back
% into CloneStats

cloneMatx = zeros(clonetotal,totalsteps);
NORMcloneMatx = zeros(clonetotal,totalsteps);
founders = 1:foundersize;

%% Set starting position
fseek(fid, startPos, 'bof');

line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
temp = textscan(line,'%f');
BoxSize = temp{1};

for i = 1:totalsteps
    tempcloneMatx = []; %matrix which will store all info for the time step 
    line = fgetl(fid);
    line = fgetl(fid);
    line = fgetl(fid);
    temp = textscan(line,'%f');
    A0 = temp{1};    
    line = fgetl(fid);
    line = fgetl(fid);
    
    while ~strcmpi(line, 'End')             %Read until keyword 'End'
        line = fgetl(fid);
        tempdata = textscan(line,'%f');
        tempcloneMatx = [tempcloneMatx;tempdata{1}'];
    end
    for j = 1:numel(tempcloneMatx(:,1))
        if ismember(tempcloneMatx(j,1),negclone) == 0 & ismember(tempcloneMatx(j,2),founders) == 1 %if a cell does not belong to one of the cells to be neglected and is a part of the selected founders, add its norm area to the master array
            NORMcloneMatx(tempcloneMatx(j,1),i) = NORMcloneMatx(tempcloneMatx(j,1),i) + tempcloneMatx(j,3)/A0; %normalized area of the cell
            cloneMatx(tempcloneMatx(j,1),i) = cloneMatx(tempcloneMatx(j,1),i) + tempcloneMatx(j,3);
        end
    end
    
    %%%% This is to read only every n-th piece of output data:
    if n ~= 1
        for j = 1:n-1
            line = fgetl(fid);
            while ~strcmpi(line, 'End')
                line = fgetl(fid);
                if ~ischar(line) %return if end of file reached
                    return
                end
            end
        end
    end
end



