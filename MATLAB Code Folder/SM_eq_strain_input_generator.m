%% generate parameter file for cluster (parameter1.txt)
% Get all text files in the current folder
files = dir('*EQ.txt');
parameter2filename = "parameters2.txt";
pfid = fopen(parameter2filename, 'wt');

% Loop through each file 
for id = 1:length(files)
    [~, f,ext] = fileparts(files(id).name);
    for deformnum = 0.001:0.001:0.01
        name = strcat(f,ext);
        newname = strcat(f, '_', num2str(deformnum,'%.3f'), ext);
        newjobname = strcat(f, '_', num2str(deformnum,'%.3f'));
        copyfile(name, strcat(f, '_', num2str(deformnum,'%.3f'), ext));
        fprintf(pfid, strcat(newname, "\t", newjobname, "\n"));
    end
    fid = fopen(strcat(f,'_dE', ext), 'wt');
    fclose(fid);
end
fclose(pfid);