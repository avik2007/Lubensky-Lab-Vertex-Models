%% generating input files (round 1)

rng = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
numCells = 100;
lim = 10;


%nomfichier= strcat('Tiss','_rng_', num2str(r,'%4d'), '_p0_', num2str(shape_index,'%.3g'), '.txt');
parameterFileName = "parameters4.txt";
fid=fopen(parameterFileName,'wt');
if fid == -1
    error('Cannot open file!');
end

files = dir('*.txt');
for id = 1:length(files)
    [~, f,ext] = fileparts(files(id).name);
    namesplit = split(f, "_"); 
    if (length(namesplit) > 6) && (length(namesplit) < 8)
        if (str2double(namesplit{5}) < 3.96) && (str2double(namesplit{5}) > 3.8)
            fprintf(fid, strcat(f,ext, '\t', f, '\n')); 
        end
    end
end
fclose(fid);

