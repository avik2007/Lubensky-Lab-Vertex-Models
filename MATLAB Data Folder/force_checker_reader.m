force_data = readmatrix('LNbC_NbPts_80_rng_100_pot_tcj_lim_3_error.txt');
%%
meanForce = force_data(100:length(force_data),1);
varianceForce = force_data(100:length(force_data),2);
%%
time = (100:length(force_data))';
%%
figure();
plot(time, meanForce);
hold on;
plot(time, varianceForce);
legend('Mean', 'Variance');
%% RECHECKING DIVISION PROBABILITY
probdata = readmatrix('probdiv.txt');
%%
histogram(probdata);