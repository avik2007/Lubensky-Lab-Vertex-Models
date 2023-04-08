
%% Read file

% solid tissue, first run
file = "N_compresses_rng_100_p0_3.75_ELS_3_EQ_iter_test_error";

ext = ".txt";
fid = fopen(strcat(file, ext), 'r');
formatSpec = '%f\t%f';
sizeA = [4 Inf];
A = fscanf(fid, formatSpec, sizeA) ;

A = A';

%% Plotting energies    
plot(A(:,1), A(:, 2) );

%% Now, a fluid sample to compare

% fluid tissue, first run
file = "N_compresses_rng_100_p0_3.75_ELS_3_EQ_compress_iter_EQ_error";

ext = ".txt";
fid = fopen(strcat(file, ext), 'r');
formatSpec = '%f\t%f';
sizeB = [4 Inf];
B = fscanf(fid, formatSpec, sizeB) ;

B = B';
%% Plotting
%hold on;
figure
plot(B(:,1), B(:, 2));

%% Just testing the relaxation post squeeze


file = "N_compresses_rng_100_p0_3.75_ELS_3_EQ_compress_iter_EQ_error";

ext = ".txt";
fid = fopen(strcat(file, ext), 'r');
formatSpec = '%f\t%f';
sizeC = [4 Inf];
C = fscanf(fid, formatSpec, sizeC) ;

C = C';

%% 
hold on;
plot(C(:,1), C(:, 2));

%% Just testing the relaxation post squeeze


file = "N_compresses_rng_100_p0_3.75_ELS_3_EQ_compress_iter_EQ_3.94_error";

ext = ".txt";
fid = fopen(strcat(file, ext), 'r');
formatSpec = '%f\t%f';
sizeC = [4 Inf];
C = fscanf(fid, formatSpec, sizeC) ;

C = C';

%% 
%%hold on;
plot(C(:,1), C(:, 2));

%% 

file = "Overlap_Central_rng_100_p0_3.75_ELS_3_EQ_SOLID_EQ_error";

ext = ".txt";
fid = fopen(strcat(file, ext), 'r');
formatSpec = '%f\t%f';
sizeD = [4 Inf];
D = fscanf(fid, formatSpec, sizeD) ;

D = D';

%% 
hold on;
plot(D(:,1), D(:, 2));

%% 

file = "Overlap_Central_rng_100_p0_3.75_ELS_3_EQ_4.03_error";

ext = ".txt";
fid = fopen(strcat(file, ext), 'r');
formatSpec = '%f\t%f';
sizeE = [4 Inf];
E = fscanf(fid, formatSpec, sizeE) ;

E = E';

%% 
%hold on;
plot(E(:,1), E(:, 2));

%%
file = "CID_400cell_rng_100_p0_3.75_ELS_3_EQ_SOLID_error";

ext = ".txt";
fid = fopen(strcat(file, ext), 'r');
formatSpec = '%f\t%f';
sizeF = [4 Inf];
F = fscanf(fid, formatSpec, sizeF) ;

F = F';

%% 
%hold on;
plot(F(:,1), F(:, 2))
%%
file = "CID_400cell_rng_100_p0_3.75_ELS_3_EQ_SOLID_500_error";

ext = ".txt";
fid = fopen(strcat(file, ext), 'r');
formatSpec = '%f\t%f';
sizeG = [4 Inf];
G = fscanf(fid, formatSpec, sizeG) ;

G = G';

%%
plot(G(:,1), G(:,2))
%%
file = "CID_400cell_rng_100_p0_3.75_ELS_3_EQ_FLUID_4.07_error";

ext = ".txt";
fid = fopen(strcat(file, ext), 'r');
formatSpec = '%f\t%f';
sizeH = [4 Inf];
H = fscanf(fid, formatSpec, sizeH) ;

H = H';

%%
plot(H(:,1), H(:,2))

%% CID_400cell_rng_100_p0_3.75_ELS_3_EQ_FLUID_4.07_100strain_error
%%
file = "CID_400cell_rng_100_p0_3.75_ELS_3_EQ_FLUID_4.07_100strain_error";

ext = ".txt";
fid = fopen(strcat(file, ext), 'r');
formatSpec = '%f\t%f';
sizeI = [4 Inf];
I = fscanf(fid, formatSpec, sizeI) ;

I = I';

%%
plot(I(:,1), I(:,2))