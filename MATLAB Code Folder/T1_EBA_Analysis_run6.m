
%% this is the round 6 (64 cells) analysis code (T1 energy barrier analysis) 
% Gather all the dE files
files = dir('*dE.txt');

Lx = sqrt(64*1250);
Ly = sqrt(64*1250); % box len gth to normalize energy by
p0num = 31;
rngnum = 20;

p0byrng = zeros(p0num, rngnum); %matrix to stor fitted (non-averaged) moduli


for id = 1:length(files)
    [~, f,ext] = fileparts(files(id).name);

    namesplit = split(f, "_");
    rng = str2double(namesplit{3});
    p0 = str2double(namesplit{5});
    rngindex = rng / 10;
    p0index = round((p0-3.7)*100 + 1);
    
    %open dE files to read out values
    fid = fopen(strcat(f, ext), 'r');
    %fid = fopen(f,'r');
    formatSpec = '%f';
    sizeA = [1 1];
    A = fscanf(fid, formatSpec, sizeA);
    A = A'; 
    if (isempty(find(A,1)))
        p0byrng(p0index, rngindex) = 0;
    else
        p0byrng(p0index, rngindex) = A(1,1) / (1.414*1250);
    end
    fclose(fid);
end

figure;

t1 = zeros(size(p0byrng,1),1);
t1var = zeros(size(p0byrng,1),1);
p0 = 3.7:0.01:4.0;
p0 = p0';
%now, we want to find averages for each row
for i = 1:size(p0byrng,1)
    row = p0byrng(i,:);
    rowdata = p0byrng(i, find(p0byrng(i,:)));
    if (length(rowdata) > 0)
        t1(i) = mean(rowdata);
        t1var(i) = std(rowdata) / sqrt(length(rowdata)) ;
    else
        t1(i) = 0;
        t1var(i) = 0;
    end
end
p0t1 = horzcat(p0, t1);
save('t1data6.mat', 'p0sm', 'p0byrng','t1var');
%save smdata p0byrng.mat
%%
errorbar(p0t1(:,1), p0t1(:,2),t1var(:)*2.086);
title('T1 Energy Barrier v Shape Index');
ylabel('T1 energy barrier', 'Interpreter', 'latex');
xlabel('$ \frac{P_0}{\sqrt{A_0}}$', 'Interpreter', 'latex');