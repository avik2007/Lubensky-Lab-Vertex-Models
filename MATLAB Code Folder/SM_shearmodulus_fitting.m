
%% this is the first round (round 1) of fitting information that is collecting information for 
% Gather all the dE files
files = dir('*EQ_dE.txt');
%f = "Tiss_rng_30_p0_3.8_EQ_dE.txt";
Lx = 790.569415;
Ly = 790.569415; % box len gth to normalize energy by
p0byrng = zeros(11,20); %matrix to stor fitted (non-averaged) moduli

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
    sizeA = [2 10];
    A = fscanf(fid, formatSpec, sizeA);
    A = A';
    sort(A,1);
    %fxx = ones(size(A,1), 1) + A(:,1);
    %strain = 1/4*( ( fxx.^4 + fxx.^(-4) ) - 2*( fxx.^2 + fxx.^(-2) ) + 2*ones(size(A,1),1) );
    epsilon = A(:,1);
    
    if (size(A,1) > 5)
        [p,S] = polyfit(epsilon, A(:,2) / (2 * Lx * Ly ) ,2);
        p0byrng(p0index, rngindex ) = p(1) / 2; 
    end
    fclose(fid);
end

figure;
scatter(epsilon, A(:,2) / (2 * Lx * Ly ));
hold on;
plot(sort(epsilon), p(1)*sort(epsilon).^2 + p(2)*sort(epsilon) + p(3));
hold off;
SM = zeros(11,1);
p0 = 3.7:0.01:3.8;
p0 = p0';
%now, we want to find averages for each row
for i = 1:size(p0byrng,1)
    row = p0byrng(i,:);
    if (length(find(row)) > 0)
        SM(i) = sum(row) / length(find(row));
    else
        SM(i) = 0;
    end
end
p0sm = horzcat(p0, SM);
save('smdata2.mat', 'p0sm');
%save smdata p0byrng.mat
%%
plot(p0, 0.05*SM,'o')