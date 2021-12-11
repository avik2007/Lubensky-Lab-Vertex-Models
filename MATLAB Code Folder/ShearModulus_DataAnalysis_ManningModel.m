%% generating input files (round 1)
%{
rng = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
numCells = 100;
lim = 10;


%nomfichier= strcat('Tiss','_rng_', num2str(r,'%4d'), '_p0_', num2str(shape_index,'%.3g'), '.txt');
parameterFileName = "parameters1.txt";
fid=fopen(parameterFileName,'wt');
if fid == -1
    error('Cannot open file!');
end

for rngIndex = 50:10:100
    for p0Index = 3.7:0.01:3.80
        GenerationVoronoi_ManningPsqrModel(rngIndex, numCells, lim, "combined", p0Index); 
        %print name of file
        fprintf(fid, strcat('Tiss','_rng_', num2str(rngIndex,'%4d'), '_p0_', num2str(p0Index,'%.3g'), '.txt'));
        %on same line, print job name
        fprintf(fid, strcat('\t T', '_rng_',  num2str(rngIndex,'%4d'), '_p0_', num2str(p0Index,'%.3g'), '\n')); 
    end
end
fclose(fid);
%}
%% generate parameter file for cluster (parameter1.txt)
% Get all text files in the current folder
files = dir('*EQ.txt');
parameter2filename = "parameters2.txt";
pfid = fopen(parameter2filename, 'wt');
a = 0.01;
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