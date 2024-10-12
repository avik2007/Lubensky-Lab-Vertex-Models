filename = 'shearmodulus_check0.6_EQenergytracking.txt'; %'shearmodulus_check0.95_2by2_EQenergytracking.txt';%'shearmodulus_check0.95_EQenergytracking.txt';%shearmodulus_check_0.8ps_EQenergytracking.txt'%'SWAB_rng_100_p0_3.76_ELS_3_EQ_EQenergytracking.txt';%
fid = fopen(filename, 'r');
%fid = fopen(f,'r');
formatSpec = '%f\t%f';
sizeA = [3 10];
A = fscanf(fid, formatSpec, sizeA);
A = A'; 
Asize = size(A);

gamma = A(:,1)';
nondenergies = A(:, 2)';

[p, S] = polyfit(gamma, nondenergies, 2);

disp(strcat(filename, " ", string(p(1) / 2)));