%% 
% Gather all the dE files

function calculateShearModulusfromDEDir(prefix, strainlevel, Lx, Ly, p0num, rngnum)
%strainlevel
%f = "Tiss_rng_30_p0_3.8_EQ_dE.txt";
if nargin < 1
    error('Give prefix for dE (does not have to match file name) folder from which you would calculate shear modulus plot');
end
if nargin < 2
    strainlevel = 1;
end 
if nargin < 3
    Lx = 353.553391; %default box length (for 100 cells)
    Ly = 353.553391; % box length to normalize energy by
end
if nargin < 5
    p0num = 16;
    rngnum = 20;
end

filename = strcat(prefix, '_dE');
files = dir(strcat( filename,'/*dE.txt'));
 
strainnum = 10;
p0byrng = zeros(p0num, rngnum); %matrix to stor fitted (non-averaged) moduli
dEstorage = zeros(p0num, rngnum, strainnum);
present = zeros(p0num, rngnum);
gof = zeros(p0num, rngnum);
for id = 1:length(files)
    [~, f,ext] = fileparts(files(id).name);
    namesplit = split(f, "_");
    %there_are_nes = contains(namesplit, "NE");
    namesplit(ismember(namesplit, {'NE'})) = [];
    if (size(namesplit,1) > 10)
        rng = str2double(namesplit{4});
        p0 = str2double(namesplit{6});
        ep = str2double(namesplit{10});
    else
        rng = str2double(namesplit{3});
        p0 = str2double(namesplit{5});
        ep = str2double(namesplit{9});
    end
    rngindex = rng / 10;
    p0index = round((p0-3.75)*100 + 1);
    %disp(f);
    %disp(p0index);
    %disp(rngindex);
    present(p0index, rngindex) = 1;
    %open dE files to read out values
    fid = fopen(strcat(filename, '/',f, ext), 'r');
    %fid = fopen(f,'r');
    formatSpec = '%f\t%f';
    sizeA = [2 3];
    A = fscanf(fid, formatSpec, sizeA);
    A = A'; 
    Asize = size(A);
    if (Asize(1,2) < 2)
        dEstorage(p0index, rngindex, round(ep*1000)) = 0;
    else
        %CHECK HOW OFTEN THIS EVENTUALITY IS NEEDED
        dEstorage(p0index, rngindex, round(A(1,2)*100*10^strainlevel)) = A(2,1) ;
        %if (strainlevel==1)
        %    dEstorage(p0index, rngindex, round(A(1,2)*1000)) = A(2,1) ;
        %else if (strainlevel==2):
        %    dEstorage(p0index, rngindex, round(A(1,2)*10000)) = A(2,1);
        %end
    end
    fclose(fid);
end
%%
for p0i = 1:1:p0num
    for rngi = 1:1:rngnum
        epsilon = [1 2 3 4 5 6 7 8 9 10];
        %if (strainlevel == 1)
        %    epsilon = epsilon' / 1000;
        %else
        %    epsilon = epsilon' / 10000;
        %end
        epsilon = epsilon' / 100*10^(strainlevel);
        straindata = squeeze(dEstorage(p0i, rngi, :));
        nonzeros = find(straindata);
        straindata = straindata(nonzeros);
        epsilon = epsilon(nonzeros); 
        if (length(straindata) > 5)
            [p,S] = polyfit(epsilon, straindata / ( Lx * Ly ) ,2  );
            %[ph, Sh] = polyfit(epsilon, normalize(straindata / (Lx*Ly) ),2);
            dE = diff(straindata / (Lx * Ly));
            eps = epsilon(1:end-1);
            p0byrng(p0i, rngi) = p(1) / 2; 
            y1 = polyval(p, epsilon, 'x');
            hold off;
            gof(p0i, rngi) = sum(  (y1-straindata / (Lx*Ly) ).^2 ) ;%/  sum( ( y1-straindata / (Lx*Ly) ) ).^2;
        end
    end
end
%%

goodfits = gof<1e-10;
p0byrng = p0byrng.*goodfits;
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
for i = 1:size(p0byrng,1)
    row = p0byrng(i,:);
    rowdata = p0byrng(i, find(p0byrng(i,:)));
    if (length(rowdata) > 0)
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
errorbar(p0sm(:,1), p0sm(:,2), SMvar, 'o')
title('Shear Modulus v Shape Index');
ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$\mathbf{\frac{P_0}{\langle A \rangle }}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName','Times New Roman');
%% easy way to see what data's getting plotted
problem_identifier = present(:,10:16) +(present(:,10:16).*goodfits(:,10:16));
figure;
x = [1:1:7];                  % Random data
y = [1:1:16];
C = problem_identifier;
xSplit = diff(x)/2;                 % Find edge points
ySplit = diff(y)/2;
xEdges = [x(1)-xSplit(1) x(2:end)-xSplit x(end)+xSplit(end)];
yEdges = [y(1)-ySplit(1) y(2:end)-ySplit y(end)+ySplit(end)];
[XGrid, YGrid] = meshgrid(xEdges,yEdges);
YGrid = flipud(YGrid);              % To match expected behavior
C = [[C zeros(size(C,1),1)] ; zeros(1,size(C,2)+1)];% Last row/col ignored
pcolor(XGrid,YGrid,C);
yticks([1:1:16])
xticklabels({'100', '110','120','130','140', '150', '160'});
yticklabels({'3.90', '3.89', '3.88', '3.87', '3.86', '3.85', '3.84', '3.83', '3.82', '3.81', '3.80', '3.79', '3.78', '3.77', '3.76', '3.75'})
hold on                             % Plot original data points
[X,Y] = meshgrid(x,y);
Y = flipud(Y);
map = [ 1 0 0
    0 0 1
    0 1 0];
colormap(map)
c = colorbar;
c.Ticks = [0.33, 1, 1.65];
c.TickLabels = {'Not present', 'T1 in fit', 'No issues'};
exportgraphics(gcf, strcat('Data_Check_Study/', prefix, '_datastudy.pdf'), 'ContentType', 'vector')
end 
