function [  ] = GenerateParameterList( fid, parameter_code, A0, compression )

%Find maximum possible int16
v=intmax('int16');
%Generate random int16
SEED = randi(v);
%Calculate characteristic length.
characteristic_length = sqrt(A0 * 2 / 3 / sqrt(3));
%Define inverse viscosities for mechanical and protein dynamics.
alphap = 100;
alpham = 100;
%Define inverse viscosity for boxsize dynamics.
TISSUE_ALPHA = 10000;
%Define cutoff lengths for T1s
MINLENGTH = 0.11 * characteristic_length;
%Calculate new edge length.
NEWEDGELENGTH = MINLENGTH * 2;
%Store tension of precursor edges.
GAMMA = 17.4859;
%Parameter for redistribution of protein after T1
KT1 = 0.1;
%Specify total time of the simulation, and the current time.
CURRENTTIME = 0.0;
%Define external tissue stress.
if compression == 0
    TISSUE_EXT_STRESS_X= GAMMA / sqrt(A0);
    TISSUE_EXT_STRESS_Y= GAMMA / sqrt(A0);
else
    TISSUE_EXT_STRESS_X= 2/3 * GAMMA / sqrt(A0);
    TISSUE_EXT_STRESS_Y= 3/2 * GAMMA / sqrt(A0) * cos(pi/6);    
end
%Define absolute tolerance for mechanics and protein dynamics.
error_tolerance_absolute_protein=0.01;

%If the you want the parameter from 2012 PLOS paper
if parameter_code == 2012
    epsab = 4 * alphap * GAMMA / characteristic_length / alpham;
    kappa = 0.93 * GAMMA / characteristic_length ^ 3;
    te = epsab * 0.033 / characteristic_length;
    jone = epsab * 0.06 * characteristic_length;
    jtwo = epsab * 0.076 * characteristic_length;
    epsnaught = epsab * 0.02;
    epsone = epsab * 0.41;
    epstwo = epsab * 0.41;
    gamma0 = 2 * GAMMA;
    gammac = GAMMA * characteristic_length;
%If the you want the parameter from 2014 PLOS paper
elseif parameter_code == 2014
    epsab = 4.86 * alphap * GAMMA / characteristic_length / alpham;
    kappa = 1.97 * GAMMA / characteristic_length ^ 3;
    te = epsab * 0.068 / characteristic_length;
    jone = epsab * 0.025 * characteristic_length;
    jtwo = epsab * 0.03 * characteristic_length;
    epsnaught = epsab * 0.02;
    epsone = epsab * 0.41;
    epstwo = epsab * 0.41;
    gamma0 = 2 * GAMMA;
    gammac = 0.92 * GAMMA * characteristic_length;
elseif parameter_code == 2017
    epsab = 1.00;
    kappa = 1.00;
    te = 0.068;
    jone = 0.025;
    jtwo = 0.03;
    epsnaught = 0.02;
    epstwo = 0.41;
    epsone = 100 * 0.41;
    gamma0 = 1.00;
    gammac = 1.00;
end
%Define protein time scale.
protein_time_scale = alphap / epsab;

%Define increment for propagation.
INCREMENT_FOR_PROPAGATION = 15.0 * protein_time_scale;

%Define max_step_size (PROBABLY ALTER THIS NOW).
max_step_size = 3.0e-04 * protein_time_scale;

%Define INITIALSTEPSIZE (PROBABLY ALTER THIS NOW).
INITIALSTEPSIZE = 5.5e-10 * protein_time_scale;

error_tolerance_absolute_mechanical = 0.01;

%Perhaps alter total time based on how many differentiation cycles.
number_of_diff_cycles = 8;
TOTALTIME = number_of_diff_cycles * INCREMENT_FOR_PROPAGATION;

%Define the time interval.
TIMEINTERVAL = TOTALTIME / 100.00;

% %Define increment for renormalization as total time (don't renormalize within a given run.
% INCREMENT_FOR_RENORMALIZATION = TOTALTIME;

%Write all of these parameters to file.
fprintf(fid, 'SEED\t%d\n', SEED);
fprintf(fid, 'TOTALTIME\t%.15g\n', TOTALTIME);
fprintf(fid, 'CURRENTTIME\t%.15g\n', CURRENTTIME);
fprintf(fid, 'TIMEINTERVAL\t%.15g\n', TIMEINTERVAL);
fprintf(fid, 'INITIALSTEPSIZE\t%.15g\n', INITIALSTEPSIZE);
fprintf(fid, 'MAX_STEP_SIZE\t%.15g\n', max_step_size);
fprintf(fid, 'error_tolerance_absolute_mechanical\t%.15g\n', error_tolerance_absolute_mechanical);
fprintf(fid, 'error_tolerance_absolute_protein\t%.15g\n', error_tolerance_absolute_protein);

% fprintf(fid, 'INCREMENT_FOR_PROPAGATION\t%f\n', INCREMENT_FOR_PROPAGATION);
% fprintf(fid, 'INCREMENT_FOR_RENORMALIZATION\t%f\n', INCREMENT_FOR_RENORMALIZATION);

fprintf(fid, 'ALPHAP\t%.15g\n', alpham);
fprintf(fid, 'TISSUE_ALPHA\t%.15g\n', TISSUE_ALPHA);
fprintf(fid, 'TISSUE_EXT_STRESS_X\t%.15g\n', TISSUE_EXT_STRESS_X);
fprintf(fid, 'TISSUE_EXT_STRESS_Y\t%.15g\n', TISSUE_EXT_STRESS_Y);
fprintf(fid, 'MINLENGTH\t%.15g\n', MINLENGTH);
fprintf(fid, 'NEWEDGELENGTH\t%.15g\n', NEWEDGELENGTH);
fprintf(fid, 'K\t%.15g\n', kappa);
fprintf(fid, 'GAMMA\t%.15g\n', GAMMA);
fprintf(fid, 'GAMMA0\t%.15g\n', gamma0);
fprintf(fid, 'GAMMAC\t%.15g\n', gammac);
fprintf(fid, 'EPSAB\t%.15g\n', epsab);
fprintf(fid, 'EPSNAUGHT\t%.15g\n', epsnaught);
fprintf(fid, 'EPSONE\t%.15g\n', epsone);
fprintf(fid, 'EPSTWO\t%.15g\n', epstwo);
fprintf(fid, 'TE\t%.15g\n', te);
fprintf(fid, 'JONE\t%.15g\n', jone);
fprintf(fid, 'JTWO\t%.15g\n', jtwo);
fprintf(fid, 'KT1\t%.15g\n', KT1);
fprintf(fid, 'endGlobal\n'); 

end

