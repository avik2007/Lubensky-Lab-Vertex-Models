function GenerateOutputMovie(inputFileName, outputFileName, n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is just a script calling ReadSimulationOutput.m, MakePBC.m, and
% PlotSimulationOutput.m to generate movie. If 3rd argument 'n' is
% provided, only draw every n frame. (Default n = 1).
%
% filename = string of file name containing the output simulation data
% (including extension .txt).
%
% Jeremy Hadidjojo (hjeremy@umich.edu), 06/07/14
%
% Last Modified on: 062014  
% List of Changes:
%   062014 - added 3rd argument 'n'
%   021716 - meryl added tif output for vewing in ImageJ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%addpath(strcat(pwd,'\MATLAB_Codes'));   %add path to where the functions are located
mkdir(outputFileName)                   %make a folder to store tif images in 

%% Open file
fid = fopen(inputFileName);
if fid == -1
    error('Cannot open file!');
end


%% Get maximum box size
[MaxLX, MaxLY] = FindMaxBoxSize(fid);



%% if 3rd argument provided, draw only every n frame
nframe = 1;      %default value: draw every frame

if nargin == 3
    nframe = n;
end



%% Initialize VideoWriter object

vidObj = VideoWriter(outputFileName, 'MPEG-4');
vidObj.FrameRate = 24;
open(vidObj);


%% Plot and make video

figure('Position', [100 50 800 800*MaxLY/MaxLX])
count = 0;
frewind(fid);
currentPos = ftell(fid);

while (true)
    [V,E,C,BoxSize, currentPos, isCellIncluded] = ReadSimulationOutput(fid, currentPos);
    %check end of file
    if currentPos == -2
        break
    end

    %draw only every 'nframe' frame
    if mod(count, nframe) == 0
        [Vpbc,Epbc,Cpbc] = MakePBC(V,E,C,BoxSize,isCellIncluded);
        clf
        %Plot simulation output with periodic boundary conditions.
        PlotSimulationOutput(Vpbc,Epbc,Cpbc,BoxSize,false,isCellIncluded,MaxLX,MaxLY)
        
        %%Added by Meryl 2/17/2016
        framename=strcat(outputFileName,'/frame_',num2str(count,'%03d'));
        saveas(gcf,framename,'tif');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        frame = getframe; 
        writeVideo(vidObj, frame);
    end
    
    count = count + 1;
end

close(vidObj);
fclose(fid);
