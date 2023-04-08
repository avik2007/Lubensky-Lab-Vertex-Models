%% 
% Gather all the dE files
files = dir('*dE.txt');
%f = "Tiss_rng_30_p0_3.8_EQ_dE.txt";
Lx = 353.553391;
Ly = 353.553391; % box length to normalize energy by

p0num = 16;
rngnum = 20;
strainnum = 10;
p0byrng = zeros(p0num, rngnum); %matrix to stor fitted (non-averaged) moduli
dEstorage = zeros(p0num, rngnum, strainnum);

for id = 1:length(files)
    [~, f,ext] = fileparts(files(id).name);

    namesplit = split(f, "_");
    rng = str2double(namesplit{3});
    p0 = str2double(namesplit{5});
    ep = str2double(namesplit{9});
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
        dEstorage(p0index, rngindex, round(ep*1000)) = 0;
    else
        dEstorage(p0index, rngindex, round(A(1,2)*1000)) = A(2,1) ;
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
        SM(i) = 0;
        SMvar(i) = 0;
    end
end
p0sm = horzcat(p0, SM);
p0smM2Ex = p0sm;
p0byrngM2Ex = p0byrng;
SMvarM2Ex = SMvar;
save('smdataM2Ex.mat', 'p0smM2Ex', 'p0byrngM2Ex','SMvarM2Ex');
%save smdata p0byrng.mat
%%
%plot(p0sm(:,1), 0.02*p0sm(:,2),'o');
errorbar(p0smM2Ex(:,1), 0.02*p0smM2Ex(:,2), 0.02*SMvarM2Ex, 'o')
title('Shear Modulus v Shape Index (Exponential, High Polydispersity)');
ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\langle A \rangle }}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName','Times New Roman');