%% 
% Gather all the dE files
files = dir('*dE.txt');
%f = "Tiss_rng_30_p0_3.8_EQ_dE.txt";
Lx = 353.553391;
Ly = 353.553391; % box length to normalize energy by
column = 15;
row = 12;

p0num = 16;
rngnum = 20;
strainnum = 10;
p0byrng = zeros(p0num, rngnum); %matrix to stor fitted (non-averaged) moduli
r2s = zeros(p0num, rngnum);
dEstorage = zeros(10, 1);
map_to_file = ( ( column - 10 ) * p0num + ( row - 1 ) ) * strainnum
for id = 1:10  
    [~, f,ext] = fileparts(files(id+map_to_file).name);

    namesplit = split(f, "_");
    rng = str2double(namesplit{4});
    p0 = str2double(namesplit{6});
    ep = str2double(namesplit{10});
    rngindex = rng / 10;
    p0index = round((p0-3.75)*100 + 1);
    
    %open dE files to read out values
    fid = fopen(strcat(f, ext), 'r');
    %fid = fopen(f,'r');
    formatSpec = '%f\t%f';
    sizeA = [2 3];
    A = fscanf(fid, formatSpec, sizeA);
    A = A'; 
    Asize = size(A);
    if (Asize(1,2) < 2)
        dEstorage(id) = 0;
    else
        dEstorage(id) = A(2,1) ;
    end
    fclose(fid);
end
%%
epsilon = [1 2 3 4 5 6 7 8 9 10];
epsilon = epsilon' / 10000;
straindata = dEstorage;
nonzeros = find(straindata);
straindata = straindata(nonzeros);
epsilon = epsilon(nonzeros);
eps = epsilon(1:end-1);
dE = diff( straindata / ( Lx * Ly) );
if (length(straindata) > 5)
    [p,S] = polyfit(epsilon, straindata / ( Lx * Ly ) ,2  );
    [pp, pgof] = polyfit(eps, dE, 1);
end    
plot(epsilon, straindata / (Lx*Ly), 'o');
x1 = linspace(0, 0.01);
hold on;
y1 = polyval(p, epsilon, 'x');
plot(epsilon, y1)
hold off;
figure;
plot(eps, dE, 'o');
hold on;
y2 = polyval(pp, eps, 'x');
plot(eps, y2);
hold off;
goodness_of_fit = sum( ( (y1-straindata / (Lx*Ly) ).^2 )) %/ ( sum(y1-straindata / (Lx*Ly) ) ).^2 

