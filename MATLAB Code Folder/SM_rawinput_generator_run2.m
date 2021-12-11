%% generating input files (round 1)

rng = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
numCells = 100;
lim = 10;


%nomfichier= strcat('Tiss','_rng_', num2str(r,'%4d'), '_p0_', num2str(shape_index,'%.3g'), '.txt');
parameterFileName = "parameters1.txt";
fid=fopen(parameterFileName,'wt');
if fid == -1
    error('Cannot open file!');
end

for rngIndex = 110:10:200
    for p0Index = 3.7:0.01:3.80
        GenerationVoronoi_ManningPsqrModel(rngIndex, numCells, lim, "combined", p0Index); 
        %print name of file
        fprintf(fid, strcat('Tiss','_rng_', num2str(rngIndex,'%4d'), '_p0_', num2str(p0Index,'%.3g'), '.txt'));
        %on same line, print job name
        fprintf(fid, strcat('\t T', '_rng_',  num2str(rngIndex,'%4d'), '_p0_', num2str(p0Index,'%.3g'), '\n')); 
    end
end
fclose(fid);

