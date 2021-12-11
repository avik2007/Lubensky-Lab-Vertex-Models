%% this is a script to calculate the shear modulus for two extreme points - p0 =3.6 and p0 = 4.0
%% 
% Gather all the dE files
files = dir('*dE.txt');
%f = "Tiss_rng_30_p0_3.8_EQ_dE.txt";
Lx = 353.553391;
Ly = 353.553391; % box len gth to normalize energy by
strainnum = 10;
strain4 = zeros(10,1);
strain36 = zeros(10,1);

for id = 1:length(files)
    [~, f,ext] = fileparts(files(id).name);

    namesplit = split(f, "_");
    rng = str2double(namesplit{3});
    p0 = str2double(namesplit{5});
    
    
    %open dE files to read out values
    fid = fopen(strcat(f, ext), 'r');
    %fid = fopen(f,'r');
    formatSpec = '%f\t%f';
    sizeA = [2 1];
    A = fscanf(fid, formatSpec, sizeA);
    A = A'; 
    if ( p0 == 4)
       strain4(round(A(1,1)*1000)) = A(1,2) / (1.414*1250);
    else
       strain36(round(A(1,1)*1000)) = A(1,2) / (1.414*1250);
    end
    fclose(fid);
end
%%
epsilon = [1 2 3 4 5 6 7 8 9 10];
epsilon = epsilon' / 1000;

[p4,S4] = polyfit(epsilon, strain4/ ( Lx * Ly ), 2  );
[p36,S36] = polyfit(epsilon, strain36 / ( Lx * Ly ) , 2);    

y4 = polyval(p4, epsilon);
y36 = polyval(p36, epsilon);
figure;
plot(epsilon, strain4 / (Lx*Ly), 'o');
hold on 
plot(epsilon,strain36 / (Lx * Ly), 'o');
plot(epsilon, y4);
plot(epsilon, y36);

title('Energy v Strain');
xlabel('$\epsilon$', 'Interpreter', 'latex');
ylabel('Normalised Energy Density');
legend('p0 = 4.0', 'p0 = 3.6', strcat('fit 4.0; SM = ',  num2str(p4(1)/ 2)),strcat('fit 3.6; SM = ', num2str(p36(1)/2) ));
hold off;
%%
SM4 = p4(1)/ 2;
SM36 = p36(1) / 2;