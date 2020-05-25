function [MaxLX, MaxLY] = FindMaxBoxSize(fid)
%% Set starting position
fseek(fid, 0, 'bof');

MaxLX = 0;
MaxLY = 0;
Lx = 0;
Ly = 0;

while true
    line = fgetl(fid);
    %return if end-of-file reached
    if ~ischar(line)
        return
    end

    while ~strcmpi(line, 'endGlobal')
        if strncmpi(line, 'LX', 2)
            temp = textscan(line, '%*s %f');
            Lx = temp{1};
        elseif strncmpi(line, 'LY', 2)
            temp = textscan(line, '%*s %f');
            Ly = temp{1};
        end
        
        %return if end-of-file reached
        line = fgetl(fid);
        if ~ischar(line)
            return
        end
        
        %compare values of Lx and Ly to those of MaxLX and MaxLY
        MaxLX = max(MaxLX, Lx);
        MaxLY = max(MaxLY, Ly);
    end
    
end