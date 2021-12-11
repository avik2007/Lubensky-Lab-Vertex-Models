
%% this is the round 3 and 4 (500 cells) of fitting information that is collecting information for 
% Gather all the dE files
files = dir('*dE.txt');

Lx = 500;
Ly = 500; % box len gth to normalize energy by
p0num = 31;
rngnum = 20;
strainnum = 10;
p0byrng = zeros(p0num, rngnum); %matrix to stor fitted (non-averaged) moduli
dEstorage = zeros(p0num, rngnum, strainnum);

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
    formatSpec = '%f\t%f';
    sizeA = [2 1];
    A = fscanf(fid, formatSpec, sizeA);
    A = A'; 
    if (isempty(find(A,1)))
        dEstorage(p0index, rngindex, round(A(1,1)*1000)) = 0;
    else
        dEstorage(p0index, rngindex, round(A(1,1)*1000)) = A(1,2) / (1.414*1250);
    end
    fclose(fid);
end
%%
for p0i = 1:1:p0num
    for rngi = 1:1:rngnum
        epsilon = [1 2 3 4 5 6 7 8 9 10];
        epsilon = epsilon' / 1000;
        straindata = squeeze(dEstorage(p0i, rngi, :));
        nonzeros = find(straindata);
        straindata = straindata(nonzeros);
        epsilon = epsilon(nonzeros);
        if (length(straindata) > 5)
            [p,S] = polyfit(epsilon, straindata / ( Lx * Ly ) ,2  );
            p0byrng(p0i, rngi) = p(1) / 2; 
        end
    end
end
    
figure;
%{
scatter(epsilon, A(:,2) / (2 * Lx * Ly ));
hold on;
plot(sort(epsilon), p(1)*sort(epsilon).^2 + p(2)*sort(epsilon) + p(3));
hold off;
%}
SM = zeros(size(p0byrng,1),1);
SMvar = zeros(size(p0byrng,1),1);
p0 = 3.7:0.01:4.0;
p0 = p0';
%now, we want to find averages for each row
for i = 1:size(p0byrng,1)
    row = p0byrng(i,:);
    rowdata = p0byrng(i, find(p0byrng(i,:)));
    if (length(rowdata) > 0)
        SM(i) = mean(rowdata);
        SMvar(i) = std(rowdata) / sqrt(length(rowdata)) ;
    else
        SM(i) = 0;
        SMvar(i) = 0;
    end
end
p0sm = horzcat(p0, SM);
save('smdata5.mat', 'p0sm', 'p0byrng','SMvar');
%save smdata p0byrng.mat
%%
plot(p0sm(:,1), 0.02*p0sm(:,2),'o');
title('Shear Modulus v Shape Index');
ylabel('$r \langle G \rangle$', 'Interpreter', 'latex');
xlabel('$ \frac{P_0}{\sqrt{A_0}}$', 'Interpreter', 'latex');