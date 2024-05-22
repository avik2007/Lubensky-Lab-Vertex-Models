%% 
% Gather all the dE files

function calculateShearModulusfromSMfiles(prefix, Lx, Ly, p0num, rngnum)
%strainlevel
%f = "Tiss_rng_30_p0_3.8_EQ_dE.txt";
if nargin < 1
    error('Give prefix for dE (does not have to match file name) folder from which you would calculate shear modulus plot');
end
 
if nargin < 3
    Lx = 353.553391; %default box length (for 100 cells)
    Ly = 353.553391; % box length to normalize energy by
end
if nargin < 5
    p0num = 16;
    rngnum = 20;
end

filename = strcat(prefix, '_SM');
files = dir(strcat( filename,'/*SM.txt'));
 
p0byrng = zeros(p0num, rngnum); %matrix to stor fitted (non-averaged) moduli
dEstorage = zeros(p0num, rngnum);
present = zeros(p0num, rngnum);

for id = 1:length(files)
    [~, f,ext] = fileparts(files(id).name);
    namesplit = split(f, "_");
    %there_are_nes = contains(namesplit, "NE");
    namesplit(ismember(namesplit, {'NE'})) = [];
    if (sum(ismember(namesplit,{'LP'})))%(size(namesplit,1) > 10)
        rng = str2double(namesplit{4});
        p0 = str2double(namesplit{6});
    else
        rng = str2double(namesplit{3});
        p0 = str2double(namesplit{5});
    end

    rngindex = rng / 10;
    p0index = round((p0-3.75)*100 + 1);

    present(p0index, rngindex) = 1;
    %open dE files to read out values
    fid = fopen(strcat(filename, '/',f, ext), 'r');
    %fid = fopen(f,'r');
    formatSpec = '%f\t%f';
    sizeA = [1];
    A = fscanf(fid, formatSpec, sizeA);
    A = A'; 
    Asize = size(A);
    if (Asize(1,1) < 1)
        dEstorage(p0index, rngindex) = 0;
    else
        %CHECK HOW OFTEN THIS EVENTUALITY IS NEEDED
        dEstorage(p0index, rngindex) = A(1,1) ;
        %if (strainlevel==1)
        %    dEstorage(p0index, rngindex, round(A(1,2)*1000)) = A(2,1) ;
        %else if (strainlevel==2):
        %    dEstorage(p0index, rngindex, round(A(1,2)*10000)) = A(2,1);
        %end
    end
    fclose(fid);
end

%%

%{
scatter(epsilon, A(:,2) / (2 * Lx * Ly ));
hold on;
plot(sort(epsilon), p(1)*sort(epsilon).^2 + p(2)*sort(epsilon) + p(3));
hold off;
%}
SM = zeros(size(p0byrng,1),1);
SMvar = zeros(size(p0byrng,1),1);
p0 = 3.75:0.01:3.90;
p0 = p0';
%now, we want to find averages for each row
for i = 1:size(p0,1)
    row = dEstorage(i,:);
    rowdata = row(find(row));
    if (length(rowdata) > 5)
        SM(i) = mean(rowdata);
        SMvar(i) = std(rowdata) / sqrt(length(rowdata)) ;
    else
        SM(i) = NaN;
        SMvar(i) = NaN;
    end
end
p0sm = horzcat(p0, SM);
p0smname = strcat('p0sm',prefix);
p0byrangename = strcat('p0byrng', prefix);
SMvarname = strcat('SMvar', prefix);
savename = strcat(filename, '/','smdata', prefix, '.mat');

feval(@()assignin('caller', p0smname, p0sm));
feval(@()assignin('caller', p0byrangename, p0byrng));
feval(@()assignin('caller', SMvarname, SMvar));

save(savename, p0smname, p0byrangename,SMvarname);

%%
%plot(p0sm(:,1), 0.02*p0sm(:,2),'o');
K = 1;
errorbar(p0sm(:,1), 0.02*p0sm(:,2) / K, 0.02*SMvar, 'o')
title('Shear Modulus v Shape Index');
ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$\mathbf{\frac{P_0}{\langle A \rangle }}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName','Times New Roman');

end 
