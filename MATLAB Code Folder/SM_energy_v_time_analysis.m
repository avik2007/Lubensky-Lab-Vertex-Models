%% reading the energy dynamics during the equilibration step in calculating shear moduli
f1 = 'Tiss_rng_200_p0_4_error';
ftemp = 'Tiss_rng_200_p0_4_error';
f2 = 'Tiss_rng_10_p0_3.6_error';
f3 = 'Tiss_rng_10_p0_4_EQ_error';
f4 = 'Tiss_rng_10_p0_4_EQ_EQ_0.001_error';
f5 = 'Tiss_rng_10_p0_4_EQ_EQ_0.001_neg_0.001_error';
f6 = 'Tiss_rng_10_p0_3.6_EQ_0.001_error';
f7 = 'Tiss_rng_10_p0_3.6_EQ_0.001_neg_0.001_error';
ext = '.txt';

fidtemp = fopen(strcat(ftemp, ext), 'r');
formatSpectemp = '%f\t%f\t%f\t%f';
sizeTemp = [5 Inf];

tdata = fscanf(fidtemp, formatSpectemp, sizeTemp);
fclose(fidtemp);
%%
fid1 = fopen(strcat(f1, ext), 'r');
    %fid = fopen(f,'r');fr
formatSpec = '%f\t%f\t%f';

sizeA = [3 Inf];


A = fscanf(fid1, formatSpec, sizeA);
fclose(fid1);

fid2 = fopen(strcat(f2,ext), 'r');
B = fscanf(fid2, formatSpec, sizeA);
fclose(fid2);

fid3 = fopen(strcat(f3, ext),'r');
C = fscanf(fid3, formatSpec, sizeA);
fclose(fid3);
fid4 = fopen(strcat(f4,ext), 'r');
D = fscanf(fid4, formatSpec, sizeA);
fclose(fid4);
fid5 = fopen(strcat(f5,ext), 'r');
E = fscanf(fid5, formatSpec, sizeA);
fclose(fid5);




Cdelay = A(1, end);
for index = 1:length(C)
   C(1,index) = C(1,index) + Cdelay; 
end
Ddelay = C(1,end);
for index = 1:length(D)
   D(1,index) = D(1,index) + Ddelay; 
end
Edelay = D(1,end);
for index = 1:length(E)
    E(1,index) = E(1,index) + Edelay;
end


A = cat(2, A, C);
A = cat(2, A, D);
A = cat(2, A, E);

%%
figure;
hold on
plot(A(1,:), A(2,:), 'o');
plot(B(1,:), B(2,:), 'o');
set(gca, 'YScale', 'log');
title('Energy vs Time');
legend('p0 = 4', 'p0 = 3.6');
xlabel('Time');
ylabel('Total Energy');
hold off;
%%
figure; 
plot(tdata(1,:), tdata(4,:), 'o');
set(gca, 'YScale', 'log');
title('Energy vs Time');
xlabel('Time');
ylabel('Total Tissue Energy');

