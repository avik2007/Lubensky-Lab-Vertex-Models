function [ V, E, C, BoxSize ] = read_data_from_vertex_model_file( filename )


%% Open file
fid = fopen(filename);
if fid == -1
    error('Cannot open file!');
end

currentPos = ftell(fid);

%%

[V, E, C, BoxSize, endPos, ~] = ReadSimulationOutputClone(fid, currentPos, 1);

if endPos == -2
    error('File is empty!');
end

%% Close file
fclose(fid);

end

