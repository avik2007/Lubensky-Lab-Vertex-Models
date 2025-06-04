function GenerateOutputPlot(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is just a script calling ReadSimulationOutput.m, MakePBC.m, and
% PlotSimulationOutput.m to generate clickable plot.
% Reading of output data file by ReadSimulationOutput will stop at the
% first 'END' keyword encountered.
%
% filename = string of file name containing the output simulation data
% (including extension .txt).
%
% Jeremy Hadidjojo (hjeremy@umich.edu), 05/23/14
%
% Last Modified on:
% List of Changes:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(strcat(pwd,'\MATLAB_Codes'));   %add path to where the functions are located


%% Open file
fid = fopen(filename);
if fid == -1
    error('Cannot open file!');
end

currentPos = ftell(fid);

%%

[V,E,C,BoxSize, endPos,isCellIncluded] = ReadSimulationOutput(fid, currentPos);

if endPos == -2
    error('File is empty!');
end
if ~isCellIncluded
    str = ['Text file does not contain "_CWVertices". ', ...
        'Cells will not be clickable and some other information may not be displayed. ', ...
        'Use Tissue::WriteOutputData4Matlab() to ensure full functionality.'];
    warning(str)
end

[Vpbc,Epbc,Cpbc] = MakePBC(V,E,C,BoxSize,isCellIncluded);

hfig = figure;
PlotSimulationOutput(Vpbc, Epbc, Cpbc, BoxSize, true, isCellIncluded);

%% Make Plot Clickable
%set custom cursor callback function (the function DataCursor_Callback.m
%must be in MATLAB's search path for this to work)

datacursormode;
hdt = datacursormode(hfig);
set(hdt,'UpdateFcn',{@DataCursor_Callback,Vpbc,Epbc,Cpbc})   %set custom callback function
set(hdt,'SnapToDataVertex','off');

%% Close file
fclose(fid);
