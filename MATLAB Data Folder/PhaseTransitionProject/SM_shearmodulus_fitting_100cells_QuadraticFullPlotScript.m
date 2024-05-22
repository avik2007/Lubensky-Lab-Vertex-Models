

% I want to non-dimensionalize the shear modulus with a pressure units
K = 50 / sqrt(1250);
%% This creates the plot that overlays the low polydispersity quadratic model shear moduli
load('L2Ts_dE/smdataL2Ts.mat')
errorbar(p0smL2Ts(:,1), 1*p0smL2Ts(:,2) / K, 1*SMvarL2Ts / K, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')

ylabel('$\mathbf{r \langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;
%grid on;
load('M2Ts_dE/smdataM2Ts.mat')
errorbar(p0smM2Ts(:,1), 1*p0smM2Ts(:,2) / K, 1*SMvarM2Ts / K, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('T2ss_dE/smdata7.mat')
errorbar(p0sm(:,1), 1*p0sm(:,2) / K, 1*SMvar / K, 'diamond', 'MarkerSize', 8, 'MarkerEdgeColor', 'red', 'Color', 'red')
%title('Shear Modulus v Shape Index (Quadratic)');

legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', 18)

savefig('L2Ts_M2Ts_T2ss_overlay.fig')
exportgraphics(gcf, 'L2Ts_M2Ts_T2ss_overlay.pdf', 'ContentType', 'vector')
hold off;
%% This creates the plot that overlays the lower polydispersity logarithmic model shear moduli
load('L2Lo_dE/smdataL2Lo.mat')

errorbar(p0smL2Lo(:,1), 1*p0smL2Lo(:,2) / K, 1*SMvarL2Lo / K, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{r \langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;
%grid on;

load('M2Lo_dE/smdataM2Lo.mat')
errorbar(p0smM2Lo(:,1), 1*p0smM2Lo(:,2) / K, 1*SMvarM2Lo / K, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('dE_LoTs/smdataLoTs.mat')
errorbar(p0smLoTs(:,1) , 1*p0smLoTs(:,2)/ K, 1*SMvarLoTs/ K, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', 18)

savefig('L2Lo_M2Lo_LoTs_overlay.fig')
exportgraphics(gcf, 'L2Lo_M2Lo_LoTs_overlay.pdf', 'ContentType', 'vector')
hold off;
%% This creates the plot that overlays the low polydispersity exponential model shear moduli

load('L2Ex_dE/smdataL2Ex.mat')

errorbar(p0smL2Ex(:,1), 1*p0smL2Ex(:,2)/ K, 1*SMvarL2Ex/ K, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{r \langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('M2Ex_dE/smdataM2Ex.mat')
errorbar(p0smM2Ex(:,1), 1*p0smM2Ex(:,2)/ K, 1*SMvarM2Ex/ K, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('ExTs_dE/smdataExTs.mat')
errorbar(p0smExTs(:,1), 1*p0smExTs(:,2)/ K, 1*SMvarExTs/ K, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', 18)

savefig('L2Ex_M2Ex_ExTs_overlay.fig')
exportgraphics(gcf, 'L2Ex_M2Ex_ExTs_overlay.pdf', 'ContentType', 'vector')
hold off;

%% This creates the plot that overlays the low polydispersity shear moduli plots where all three plots have A0 = <A>

load('T2ss_dE/smdata7.mat')
errorbar(p0sm(:,1), 1*p0sm(:,2)/ K, 1*SMvar/ K, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')

ylabel('$\mathbf{r \langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('dE_LoTs/smdataLoTs.mat')
errorbar(p0smLoTs(:,1), 1*p0smLoTs(:,2)/ K, 1*SMvarLoTs/ K, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('ExTs_dE/smdataExTs.mat')
errorbar(p0smExTs(:,1), 1*p0smExTs(:,2)/ K, 1*SMvarExTs/ K, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('Quadratic $E^A$','Logarithmic $E^A$ ', 'Exponential $E^A$' , 'Interpreter', 'latex' ,'FontSize', 18)
savefig('T2ss_LoTs_ExTs_overlay.fig')
exportgraphics(gcf, 'T2ss_LoTs_ExTs_overlay.pdf', 'ContentType', 'vector')

%% This creates the plot that overlays the high polydiserpsity quadratic shear moduli plot

load('L2Ts_LP_dE_2.0/smdataL2TsHP.mat') %-- I never properly made this up
errorbar(p0smL2TsHP(:,1), 1*p0smL2TsHP(:,2)/ K, 1*SMvarL2TsHP/ K, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{r \langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('M2Ts_LP_dE/smdataM2TsHP.mat')
errorbar(p0smM2TsHP(:,1), 1*p0smM2TsHP(:,2)/ K, 1*SMvarM2TsHP/ K, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('T2ss_LP_dE_2.0/smdataTshp.mat')
errorbar(p0smTsHP(:,1), 1*p0smTsHP(:,2)/ K, 1*SMvarTsHP/ K, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', 18)

% RERUN L2TS_LP_DE AND GET THAT FIXED UP
savefig('L2Ts_M2Ts_T2ss_HighPoly_overlay.fig')
exportgraphics(gcf, 'L2Ts_M2Ts_T2ss_HighPoly_overlay.pdf', 'ContentType', 'vector')

%%  This creates the plot that overlays the high polydispersity logarithmic shear moduli plot
K = 50 / sqrt(1250);
load('L2Lo_LP_dE/smdataL2LoHP.mat')
errorbar(p0smL2LoHP(:,1), 1*p0smL2LoHP(:,2)/ K, 1*SMvarL2LoHP/ K, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{r \langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('M2Lo_LP_dE/smdataM2LoHP.mat')
errorbar(p0smM2LoHP(:,1), 1*p0smM2LoHP(:,2)/ K, 1*SMvarM2LoHP/ K, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('LoTs_LP_dE/smdataLoTsHP.mat')
errorbar(p0smLoTsHP(:,1), 1*p0smLoTsHP(:,2)/ K, 1*SMvarLoTsHP/ K, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', 18)

savefig('L2Lo_M2Lo_LoTs_HighPoly_overlay.fig')
exportgraphics(gcf, 'L2Lo_M2Lo_LoTs_HighPoly_overlay.pdf', 'ContentType', 'vector')


%% This creates the plot that overlays the high polydispsersity exponential shear moduli plot

load('L2Ex_LP_dE/smdataL2ExHP.mat')
errorbar(p0smL2ExHP(:,1), 1*p0smL2ExHP(:,2)/ K, 1*SMvarL2ExHP/ K, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{r \langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName','Times New Roman');
hold on;

load('M2Ex_LP_dE/smdataM2ExHP.mat')
errorbar(p0smM2ExHP(:,1), 1*p0smM2ExHP(:,2)/ K, 1*SMvarM2ExHP/ K, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('ExTs_LP_dE/smdataExTsHP.mat')
errorbar(p0smExTsHP(:,1), 1*p0smExTsHP(:,2)/ K, 1*SMvarExTsHP/ K, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', 18)


savefig('L2Ex_M2Ex_ExTs_HighPoly_overlay.fig')
exportgraphics(gcf, 'L2Ex_M2Ex_ExTs_HighPoly_overlay.pdf', 'ContentType', 'vector')

%% This creates the plot that overlays the high polydispersity shear moduli plots where all three plots have A0 = <A>

load('T2ss_LP_dE_2.0/smdataTshp.mat')
errorbar(p0smTsHP(:,1), 1*p0smTsHP(:,2)/ K, 1*SMvarTsHP/ K, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')

ylabel('$\mathbf{r \langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$\mathbf{\frac{P_0}{A_{avg} }}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;
Tmore
load('LoTs_LP_dE/smdataLoTsHP.mat')
errorbar(p0smLoTsHP(:,1), 1*p0smLoTsHP(:,2)/ K, 1*SMvarLoTsHP/ K, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')
%ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
%xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');


load('ExTs_LP_dE/smdataExTsHP.mat')
errorbar(p0smExTsHP(:,1), 1*p0smExTsHP(:,2)/ K, 1*SMvarExTsHP/ K, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('Quadratic $E^A$','Logarithmic $E^A$ ', 'Exponential $E^A$' , 'Interpreter', 'latex' ,'FontSize', 18)
savefig('T2ss_LoTs_ExTs_HighPoly_overlay.fig')
exportgraphics(gcf, 'T2ss_LoTs_ExTs_overlay_HighPoly_overlay.pdf', 'ContentType', 'vector')

%% This combines the plots that show A0=<A> for high and low polydispersity cases
K = 50 / sqrt(1250);
subgplotlabelindex_x = 3.73;
subgplotlabelindex_y = 3.7;
plot_width = 1200;
plot_height = 360;
markersize = 6;
fontsize = 18;
titlefontsize = 16;
subplot(1,2,1); %low polydispersity plot first 

load('T2ss_dE/smdata7.mat')
errorbar(p0sm(:,1), 1*p0sm(:,2)/ K, 1*SMvar/ K, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'black', 'Color', 'black')
%yticks('manual')

hF = gcf;
hF.Position(1:4) = [0 0 plot_width plot_height];

ylabel('$\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName','Times New Roman');
hold on;
grid on;
load('LoTs_dE/smdataLoTs.mat')
errorbar(p0smLoTs(:,1), 1*p0smLoTs(:,2)/ K, 1*SMvarLoTs/ K, 'square', 'MarkerSize', markersize, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('ExTs_dE/smdataExTs.mat')
errorbar(p0smExTs(:,1), 1*p0smExTs(:,2)/ K, 1*SMvarExTs/ K, 'diamond', 'MarkerSize', markersize, 'MarkerEdgeColor', 'red', 'Color', 'red')
ylim([-0.5 4.0])
yticks(-0.5:0.5:4.0)
yticklabels({'-0.5','0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0','3.5', '4.0'})
xticks(3.75:0.03:3.9)
legend('Quadratic $E^A$','Logarithmic $E^A$ ', 'Exponential $E^A$' , 'Interpreter', 'latex' ,'FontSize', fontsize)
title('$\textbf{Low Polydispersity}$, $\mathbf{A_0 = A_{avg}}$', 'Interpreter', 'latex' ,'FontSize', titlefontsize, 'FontName','Times New Roman');
text(subgplotlabelindex_x ,subgplotlabelindex_y , '$\textbf{a}$', 'FontSize', 20,'Interpreter', 'latex');
%------------------------------------------------------------------
subplot(1,2,2); %then the high polydispersity case

load('T2ss_LP_dE_2.0/smdataTshp.mat')
errorbar(p0smTsHP(:,1), 1*p0smTsHP(:,2)/ K, 1*SMvarTsHP/ K, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName','Times New Roman');
hold on;
grid on;
load('LoTs_LP_dE/smdataLoTsHP.mat')

errorbar(p0smLoTsHP(:,1), 1*p0smLoTsHP(:,2)/ K, 1*SMvarLoTsHP/ K, 'square', 'MarkerSize', markersize, 'MarkerEdgeColor', 'blue', 'Color', 'blue')
%yticks(-0.01:0.01:0.12)

load('ExTs_LP_dE/smdataExTsHP.mat')
errorbar(p0smExTsHP(:,1), 1*p0smExTsHP(:,2)/ K, 1*SMvarExTsHP/ K, 'diamond', 'MarkerSize', markersize, 'MarkerEdgeColor', 'red', 'Color', 'red')
ylim([-0.5 4.0])
yticks(-0.5:0.5:4.0)
yticklabels({'-0.5','0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0','3.5', '4.0'})
xticks(3.75:0.03:3.9)
legend('Quadratic $E^A$','Logarithmic $E^A$ ', 'Exponential $E^A$' , 'Interpreter', 'latex' ,'FontSize', fontsize)
title('$\textbf{High Polydispersity}$, $\mathbf{A_0 = A_{avg}}$', 'Interpreter', 'latex' ,'FontSize', 16, 'FontName','Times New Roman');
text(subgplotlabelindex_x, subgplotlabelindex_y, '$\textbf{b}$', 'FontSize', 20,'Interpreter', 'latex');

savefig('T2ss_LoTs_ExTs_overlay_HP_AND_LP.fig')
exportgraphics(gcf, 'T2ss_LoTs_ExTs_overlay_HP_AND_LP.pdf', 'ContentType', 'vector')

%% This will combine the 6 plots that show the three different types of potentials for the three different values of A0 for high and low polydispersity 
%------------------------------------------------------------------
K = 50 / sqrt(1250);

subgplotlabelindex_x = 3.73;
subgplotlabelindex_y = 3.7;
subgplotlabelfontsize = 20;
plot_width = 1200;
plot_height = 1200;
horizontal_title_x = 3.715;
horizontal_title_y = 0.5;
quadratic_title_y_offset = 0.1;
log_title_y_offset = -0.1;
exp_title_y_offset = -0.1;
title_font = 18;
legend_font = 18;
markersize = 6;
xlabelfont = 16;
ylabelfont = 14;
fontsize = 18;
subplot(3,2,1) % first low polydispersity Tiss case

hF = gcf;
hF.Position(1:4) = [0 0 plot_width plot_height];
 
load('L2Ts_dE/smdataL2Ts.mat')
errorbar(p0smL2Ts(:,1), 1*p0smL2Ts(:,2)/ K, 1*SMvarL2Ts/ K, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'black', 'Color', 'black')

ylabel('$\mathbf{ \langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', ylabelfont, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', xlabelfont, 'FontName','Times New Roman');

hold on;
grid on;
load('M2Ts_dE/smdataM2Ts.mat')
errorbar(p0smM2Ts(:,1), 1*p0smM2Ts(:,2)/ K, 1*SMvarM2Ts/ K, 'square', 'MarkerSize', markersize, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('T2ss_dE/smdata7.mat')
errorbar(p0sm(:,1), 1*p0sm(:,2)/ K, 1*SMvar/ K, 'diamond', 'MarkerSize', markersize, 'MarkerEdgeColor', 'red', 'Color', 'red')
%title('Shear Modulus v Shape Index (Quadratic)');

legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', legend_font)
ylim([-0.5 4.0])
yticks(-0.5:0.5:4.0)
yticklabels({'-0.5','0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0','3.5', '4.0'})
xticks(3.75:0.03:3.9)

title('$\textbf{Low Polydispersity}$', 'Interpreter', 'latex' ,'FontSize', title_font, 'FontName','Times New Roman');
%subtitle('$\textbf{Quadratic } \mathbf{E^A}$', 'Interpreter', 'latex' ,'FontSize', 14, 'FontName','Times New Roman');
text(horizontal_title_x, horizontal_title_y+quadratic_title_y_offset, '$\textbf{Quadratic}$', 'FontSize', title_font,'Interpreter', 'latex', 'Rotation', 90);
text(subgplotlabelindex_x, subgplotlabelindex_y, '$\textbf{a}$', 'FontSize', subgplotlabelfontsize,'Interpreter', 'latex');

%------------------------------------------------------------------
subplot(3,2,3)  % second, low polydisperisty LoTs case

load('L2Lo_dE/smdataL2Lo.mat')

errorbar(p0smL2Lo(:,1), 1*p0smL2Lo(:,2)/ K, 1*SMvarL2Lo/ K, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'black', 'Color', 'black')
%ylabel('$\mathbf{r \langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', xlabelfont, 'FontName','Times New Roman');
ylabel(' $\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', ylabelfont, 'FontName', 'Times New Roman');
%subtitle('$\textbf{Logarithmic } \mathbf{E^A}$', 'Interpreter', 'latex' ,'FontSize', 14, 'FontName','Times New Roman');
hold on;
grid on;

load('M2Lo_dE/smdataM2Lo.mat')
errorbar(p0smM2Lo(:,1), 1*p0smM2Lo(:,2)/ K, 1*SMvarM2Lo/ K, 'square', 'MarkerSize', markersize, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('LoTs_dE/smdataLoTs.mat')
errorbar(p0smLoTs(:,1), 1*p0smLoTs(:,2)/ K, 1*SMvarLoTs/ K, 'diamond', 'MarkerSize', markersize, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', legend_font)
ylim([-0.5 4.0])
yticks(-0.5:0.5:4.0)
yticklabels({'-0.5','0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0','3.5', '4.0'})
xticks(3.75:0.03:3.9)
text(horizontal_title_x, horizontal_title_y+log_title_y_offset, '$\textbf{Logarithmic}$', 'FontSize', title_font,'Interpreter', 'latex', 'Rotation', 90);
text(subgplotlabelindex_x, subgplotlabelindex_y, '$\textbf{c}$', 'FontSize', subgplotlabelfontsize,'Interpreter', 'latex');
%------------------------------------------------------------------
subplot(3,2,5)  % third, low polydisperisty ExTs case


load('L2Ex_dE/smdataL2Ex.mat')

errorbar(p0smL2Ex(:,1), 1*p0smL2Ex(:,2)/ K, 1*SMvarL2Ex/ K, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', ylabelfont, 'FontName', 'Times New Roman');
%ylabel({'$\textbf{Exponential}$ ',' $\mathbf{r \langle G \rangle}$'}, 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', xlabelfont, 'FontName','Times New Roman');
%subtitle('$\textbf{Exponential } \mathbf{E^A}$', 'Interpreter', 'latex' ,'FontSize', 14, 'FontName','Times New Roman');
hold on;
grid on;
load('M2Ex_dE/smdataM2Ex.mat')
errorbar(p0smM2Ex(:,1), 1*p0smM2Ex(:,2)/ K, 1*SMvarM2Ex/ K, 'square', 'MarkerSize', markersize, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('ExTs_dE/smdataExTs.mat')
errorbar(p0smExTs(:,1), 1*p0smExTs(:,2)/ K, 1*SMvarExTs/ K, 'diamond', 'MarkerSize', markersize, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', legend_font)

ylim([-0.5 4.0])
yticks(-0.5:0.5:4.0)
yticklabels({'-0.5','0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0','3.5', '4.0'})
xticks(3.75:0.03:3.9)
text(horizontal_title_x, horizontal_title_y+exp_title_y_offset, '$\textbf{Exponential}$', 'FontSize', title_font,'Interpreter', 'latex', 'Rotation', 90);
text(subgplotlabelindex_x, subgplotlabelindex_y, '$\textbf{e}$', 'FontSize', subgplotlabelfontsize,'Interpreter', 'latex');
%------------------------------------------------------------------
subplot(3,2,2); % now going to the second column with the high polydispersity T2ss plot

load('L2Ts_LP_dE_2.0/smdataL2TsHP.mat') %-- I never properly made this up
errorbar(p0smL2TsHP(:,1), 1*p0smL2TsHP(:,2)/ K, 1*SMvarL2TsHP/ K, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('M2Ts_LP_dE/smdataM2TsHP.mat')
errorbar(p0smM2TsHP(:,1), 1*p0smM2TsHP(:,2)/ K, 1*SMvarM2TsHP/ K, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('T2ss_LP_dE_2.0/smdataTshp.mat')
errorbar(p0smTsHP(:,1), 1*p0smTsHP(:,2)/ K, 1*SMvarTsHP/ K, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')


ylabel('$\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
title('$\textbf{High Polydispersity}$', 'Interpreter', 'latex' ,'FontSize', 18, 'FontName','Times New Roman');
%subtitle('$\textbf{Quadratic } \mathbf{E^A}$', 'Interpreter', 'latex' ,'FontSize', 14, 'FontName','Times New Roman');
legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', 18)
ylim([-0.5 4.0])
yticks(-0.5:0.5:4.0)
yticklabels({'-0.5','0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0','3.5', '4.0'})
xticks(3.75:0.03:3.9)
text(subgplotlabelindex_x, subgplotlabelindex_y, '$\textbf{b}$', 'FontSize', 20,'Interpreter', 'latex');
grid on;

%-------------------------------------------------------------------
subplot(3,2,4); % this is the high polydispiersity log plot

load('L2Lo_LP_dE/smdataL2LoHP.mat')
errorbar(p0smL2LoHP(:,1), 1*p0smL2LoHP(:,2)/ K, 1*SMvarL2LoHP/ K, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', ylabelfont, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', xlabelfont, 'FontName','Times New Roman');
%subtitle('$\textbf{Logarithmic } \mathbf{E^A}$', 'Interpreter', 'latex' ,'FontSize', 14, 'FontName','Times New Roman');
hold on;
grid on;
load('M2Lo_LP_dE/smdataM2LoHP.mat')
errorbar(p0smM2LoHP(:,1), 1*p0smM2LoHP(:,2)/ K, 1*SMvarM2LoHP/ K, 'square', 'MarkerSize', markersize, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('LoTs_LP_dE/smdataLoTsHP.mat')
errorbar(p0smLoTsHP(:,1), 1*p0smLoTsHP(:,2)/ K, 1*SMvarLoTsHP/ K, 'diamond', 'MarkerSize', markersize, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', legend_font)
ylim([-0.5 4.0])
yticks(-0.5:0.5:4.0)
yticklabels({'-0.5','0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0','3.5', '4.0'})
xticks(3.75:0.03:3.9)
text(subgplotlabelindex_x, subgplotlabelindex_y, '$\textbf{d}$', 'FontSize', subgplotlabelfontsize,'Interpreter', 'latex');
%-------------------------------------------------------------------
subplot(3,2,6); %high polydispiersity expoential plot

load('L2Ex_LP_dE/smdataL2ExHP.mat')
errorbar(p0smL2ExHP(:,1), 1*p0smL2ExHP(:,2)/ K, 1*SMvarL2ExHP/ K, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel(' $\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', ylabelfont, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', xlabelfont, 'FontName','Times New Roman');
%subtitle('$\textbf{Exponential } \mathbf{E^A}$', 'Interpreter', 'latex' ,'FontSize', 14, 'FontName','Times New Roman');
hold on;
grid on;
load('M2Ex_LP_dE/smdataM2ExHP.mat')
errorbar(p0smM2ExHP(:,1), 1*p0smM2ExHP(:,2)/ K, 1*SMvarM2ExHP/ K, 'square', 'MarkerSize', markersize, 'MarkerEdgeColor', 'blue', 'Color', 'blue')
load('ExTs_LP_dE/smdataExTsHP.mat')
errorbar(p0smExTsHP(:,1), 1*p0smExTsHP(:,2)/ K, 1*SMvarExTsHP/ K, 'diamond', 'MarkerSize', markersize, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', legend_font)
ylim([-0.5 4.0])
yticks(-0.5:0.5:4.0)
yticklabels({'-0.5','0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0','3.5', '4.0'})
xticks(3.75:0.03:3.9)
text(subgplotlabelindex_x, subgplotlabelindex_y, '$\textbf{f}$', 'FontSize', subgplotlabelfontsize,'Interpreter', 'latex');

savefig('quad_log_exp_cum_overlay.fig')
exportgraphics(gcf, 'quad_log_exp_cum_overlay.pdf', 'ContentType', 'vector')

%% Polydispersity plot

la05 = [0.215695; 0.184582; 0.198962; 0.190442; 0.190386; 0.213417; 0.172848];
la10 = [0.13892;0.116656; 0.123921; 0.117274; 0.11225; 0.143902; 0.110445];
la = [la05 la10];
boxplot(la, 'Labels', {'high polydispsersity', 'low polydispersity'})
ylabel('\boldmath${\sigma_a}$', 'Interpreter', 'latex', 'FontSize', 28)
set(gca,'XTickLabel',{'high polydispersity', 'low polydispersity'},'FontSize',xlabelfont-2, 'FontName', 'Times New Roman')
set(gca,'YTickLabel', get(gca, 'YTickLabel'), 'FontSize', ylabelfont, 'FontName', 'Times New Roman')
%title('$\textbf{Polydispersity Distribution}$', 'Interpreter', 'latex' ,'FontSize', 18, 'FontName','Times New Roman');

savefig('polydispersity_distribution.fig')
exportgraphics(gcf, 'polydispersity_distribution.pdf', 'ContentType', 'vector')%% Polydispersity plot

%% Polydispersity plot with no lloyd iteration for reference
la05 = [0.215695; 0.184582; 0.198962; 0.190442; 0.190386; 0.213417; 0.172848];
la10 = [0.13892;0.116656; 0.123921; 0.117274; 0.11225; 0.143902; 0.110445];
la01 = [0.496707; 0.511931; 0.528366; 0.443767; 0.511317; 0.519427; 0.52317];
la = [la01 la05 la10];
boxplot(la, 'Labels', {'no Lloyd iters', 'high polydispsersity', 'low polydispersity'})
ylabel('\boldmath${\sigma_a}$', 'Interpreter', 'latex', 'FontSize', 28)
set(gca,'XTickLabel',{'no Lloyd iters', 'high polydispersity', 'low polydispersity'},'FontSize',xlabelfont-2, 'FontName', 'Times New Roman')
set(gca,'YTickLabel', get(gca, 'YTickLabel'), 'FontSize', ylabelfont, 'FontName', 'Times New Roman')
%title('$\textbf{Polydispersity Distribution}$', 'Interpreter', 'latex' ,'FontSize', 18, 'FontName','Times New Roman');

savefig('polydispersity_distribution_ref.fig')
exportgraphics(gcf, 'polydispersity_distribution_ref.pdf', 'ContentType', 'vector')

%% Polydispersity plot more than just la05 and la10
subgplotlabelindex_x = 3.73;
subgplotlabelindex_y = 3.7;
subgplotlabelfontsize = 20;
plot_width = 1200;
plot_height = 1200;
horizontal_title_x = 3.715;
horizontal_title_y = 0.5;
quadratic_title_y_offset = 0.1;
log_title_y_offset = -0.1;
exp_title_y_offset = -0.1;
title_font = 18;
legend_font = 18;
markersize = 6;
xlabelfont = 16;
ylabelfont = 14;
fontsize = 18;

la05 = [0.215695; 0.184582; 0.198962; 0.190442; 0.190386; 0.213417; 0.172848];
la10 = [0.13892; 0.116656; 0.123921; 0.117274; 0.11225; 0.143902; 0.110445];
la04 = [0.245162; 0.21775; 0.229866; 0.219575; 0.223599; 0.241915; 0.205338]; 
la03 = [0.287166; 0.268531; 0.275529; 0.258159; 0.271914; 0.281516; 0.261287 ];
la02 = [0.354672; 0.349919; 0.353379; 0.316223; 0.349813; .346665; 0.358769;];
la01 = [0.496707; 0.511931; 0.528366; 0.443767; 0.511317; 0.519427; 0.52317];
la = [la01 la02 la03 la04 la05 la10];
boxplot(la, 'Labels', {'1', '2', '3', '4', '5', '10'})
ylabel('\boldmath${\sigma_a}$', 'Interpreter', 'latex', 'FontSize', 28)
xlabel('Number of Lloyd Iterations', 'FontSize', 28)
set(gca,'XTickLabel',{'1','2', '3', '4', '5', '10'},'FontSize',xlabelfont-2, 'FontName', 'Times New Roman')
set(gca,'YTickLabel', get(gca, 'YTickLabel'), 'FontSize', ylabelfont, 'FontName', 'Times New Roman')
%title('$\textbf{Polydispersity Distribution}$', 'Interpreter', 'latex' ,'FontSize', 18, 'FontName','Times New Roman');

savefig('polydispersity_distribution_higher.fig')
exportgraphics(gcf, 'polydispersity_distribution_higher.pdf', 'ContentType', 'vector')
%% Sweeping Abox plot
subgplotlabelindex_x = 3.73;
subgplotlabelindex_y = 3.7;
subgplotlabelfontsize = 20;
plot_width = 600;
plot_height = 360;
horizontal_title_x = 3.715;
horizontal_title_y = 0.5;
quadratic_title_y_offset = 0.1;
log_title_y_offset = -0.1;
exp_title_y_offset = -0.1;
title_font = 18;
legend_font = 18;
markersize = 6;
xlabelfont = 16;
ylabelfont = 14;
fontsize = 18;
K = 50 / sqrt(1250);
figure
t = tiledlayout(1,1);
load('SWAB_dE/smdataSWAB.mat');

ax1 = axes(t);
errorbar(ax1, p0smSWAB(:,1), 1*p0smSWAB(:,2)/ K, 1*SMvarSWAB/ K, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor',[0 0.4 0], 'Color', [0 0.4 0]);

ylim(ax1, [-0.5 2.6])
yticks(ax1, -0.5:0.5:2.5)
%yticklabels(ax1, {'-0.5','0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0'})
xticks(ax1, 3.75:0.03:3.9)
grid on;
ax2 = axes(t);

ax2.XAxisLocation = 'top';
set(ax2, 'color', 'none');
set(ax2, 'ytick', []);
ax2.XDir = 'reverse';
S0 = 3.75:0.03:3.90;

P0 = 135;
Abox = P0^2 ./ S0.^2;
xlim(ax2, [1198 1296])
xticks(ax2, 1198:(1296-1198)/5:1296);
%xticklabels(ax2, {'1296', '1276', '1256', '1236', '1217', '1198'})

ylabel(ax1, '$\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', 'Times New Roman');
xlabel(ax1, '$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName','Times New Roman');
xlabel(ax2, '$ \mathbf{A_{box}}$ (pixels)', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', 'Times New Roman');

savefig('SWAB.fig')
exportgraphics(gcf, 'SWAB.pdf', 'ContentType', 'vector')
%% Sweeping Abox plot with comparison for reference
subgplotlabelindex_x = 3.73;
subgplotlabelindex_y = 3.7;
subgplotlabelfontsize = 20;
plot_width = 600;
plot_height = 360;
horizontal_title_x = 3.715;
horizontal_title_y = 0.5;
quadratic_title_y_offset = 0.1;
log_title_y_offset = -0.1;
exp_title_y_offset = -0.1;
title_font = 18;
legend_font = 18;
markersize = 6;
xlabelfont = 16;
ylabelfont = 14;
fontsize = 18;
K = 50 / sqrt(1250);
figure
t = tiledlayout(1,1);
load('SWAB_dE/smdataSWAB.mat');
load('T2ss_dE/smdata7.mat')

ax1 = axes(t); 

%errorbar(ax1, p0smSWAB(:,1), 1*p0smSWAB(:,2)/ K, 1*SMvarSWAB/ K, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor',[0 0.4 0], 'Color', [0 0.4 0]);
errorbar(ax1, p0sm(:,1), 1*p0sm(:,2) / K, 1*SMvar / K, 'diamond', 'MarkerSize', markersize, 'MarkerEdgeColor', 'r', 'Color', 'r')
hold on;
errorbar(ax1, p0smSWAB(:,1), 1*p0smSWAB(:,2)/ K, 1*SMvarSWAB/ K, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor',[0.4660 0.6740 0.2880], 'Color', [0.4660 0.6740 0.2880]);
legend('sweeping $P_0$ ','sweeping $A_{box}$', 'Interpreter', 'latex' ,'FontSize', legend_font)

ylim(ax1, [-0.5 2.6])
yticks(ax1, -0.5:0.5:2.5)
%yticklabels(ax1, {'-0.5','0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0'})
xticks(ax1, 3.75:0.03:3.9)
grid on;
ax2 = axes(t);

ax2.XAxisLocation = 'top';
set(ax2, 'color', 'none');
set(ax2, 'ytick', []);
ax2.XDir = 'reverse';
S0 = 3.75:0.03:3.90;

P0 = 135;
Abox = P0^2 ./ S0.^2;
xlim(ax2, [1198 1296])
xticks(ax2, 1198:(1296-1198)/5:1296);
%xticklabels(ax2, {'1296', '1276', '1256', '1236', '1217', '1198'})

ylabel(ax1, '$\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', 'Times New Roman');
xlabel(ax1, '$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName','Times New Roman');
xlabel(ax2, '$ \mathbf{A_{box}}$ (pixels)', 'Interpreter', 'latex', 'FontSize', fontsize, 'FontName', 'Times New Roman');

savefig('SWAB_with_ref.fig')
exportgraphics(gcf, 'SWAB_with_ref.pdf', 'ContentType', 'vector')

%% T2ss, LBTS, BBTs
legend_font = 18;

K1 =  50 / sqrt(1250);
LAMBDA = 70.710678;
load('T2ss_dE/smdata7.mat')
errorbar(p0sm(:,1), 1*p0sm(:,2)/ K1, 1*SMvar/ K1, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')

ylabel('$\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

K2 = K1; %50 / sqrt(625);

load('LBTs_dE/smdataLBTs.mat')
errorbar(p0smLBTs(:,1), 1*p0smLBTs(:,2)/ K2, 1*SMvarLBTs/ K2, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

K3 =  K1; %50 / sqrt(2500);

load('BBTs_dE/smdataBBTs.mat')
errorbar(p0smBBTs(:,1), 1*p0smBBTs(:,2)/ K3, 1*SMvarBBTs/ K3, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')


ylim([-0.5 3.0])
yticks(-0.5:0.5:2.5)
yticklabels({'-0.5','0', '0.5', '1.0', '1.5', '2.0', '2.5'})
xticks(3.75:0.03:3.9)

legend('$A_{avg} = A_0$','$A_{avg} = 0.5 A_0 $ ', '$A_{avg} = 2 A_0 $ ' , 'Interpreter', 'latex' ,'FontSize', legend_font)
grid on;

savefig('T2ss_LBTs_BBTs_overlay.fig')
exportgraphics(gcf, 'T2ss_LBTs_BBTs_overlay.pdf', 'ContentType', 'vector')
%% ExTs, LBEx, BBEx with inverse perimeter modulus 
legend_font = 18;

K = 50 / sqrt(1250)^3; 
A0 = 1250;
shape_index = 0.02;
LAMBDA = K *A0 / shape_index;

r1 = 1;%K*A0 / LAMBDA;
load('ExTs_dE/smdataExTs.mat')
errorbar(p0smExTs(:,1), r1*p0smExTs(:,2), r1*SMvarExTs, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')

ylabel('$r \mathbf{\langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

K2 = 50 / sqrt(625);
%r2 = K2 / LAMBDA * sqrt(1250 / 625);
r2 = 1;%r1*(A0 / 625);
load('LBEx_dE/smdataLBEx.mat')
errorbar(p0smLBEx(:,1), r2*p0smLBEx(:,2), r2*SMvarLBEx, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

K3 = 50 / sqrt(2500);
%r3 = K3 / LAMBDA * sqrt(1250 / 2500);
r3 = 1;%r1*(A0 / 2500);
load('BBEx_dE/smdataBBEx.mat')
errorbar(p0smBBEx(:,1), r3*p0smBBEx(:,2), r3*SMvarBBEx, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_{avg} = A_0$','$A_{avg} = 0.5 A_0 $ ', '$A_{avg} = 2 A_0 $ ' , 'Interpreter', 'latex' ,'FontSize', legend_font)

savefig('ExTs_LBEx_BBEx_overlay_scaled_r.fig')
exportgraphics(gcf, 'ExTs_LBEx_BBEx_overlay_scaled_r.fig.pdf', 'ContentType', 'vector')
%% linear fit of LBTs, BBTs, T2ss
load('T2ss_dE/smdata7.mat')
load('LBTs_dE/smdataLBTs.mat')
load('BBTs_dE/smdataBBTs.mat')
l = 6;
K1 =  50 / sqrt(1250);
LAMBDA = 70.710678;
r = 0.02;
pBB = polyfit(3.813 - p0smBBTs(1:l,1), p0smBBTs(1:l,2)*r, 1);
pLB = polyfit(3.813 - p0smLBTs(1:l,1), p0smLBTs(1:l,2)*r, 1);
pT2 = polyfit(3.813 - p0sm(1:l,1), p0sm(1:l,2)*r, 1);

pBB_slope = pBB(1);
pLB_slope = pLB(1);
pT2_slope = pT2(1);


%% linear fit of LBEx, BBEx, ExTs
load('BBEx_dE/smdataBBEx.mat')
load('ExTs_dE/smdataExTs.mat')
load('LBEx_dE/smdataLBEx.mat')
K1 =  50 / sqrt(1250);
LAMBDA = 70.710678;
r = 0.02;
potential_slope = LAMBDA / K1;
l = 6;
pBB = polyfit(p0smBBEx(1:l,1), p0smBBEx(1:l,2) *r, 1);
pLB = polyfit(p0smLBEx(1:l,1), p0smLBEx(1:l,2) *r, 1);
pT2 = polyfit(p0smExTs(1:l,1), p0smExTs(1:l,2) *r, 1);

%% linear fit of LP functions

load('L2Ex_LP_dE/smdataL2ExHP.mat')
load('M2Ex_LP_dE/smdataM2ExHP.mat')
load('ExTs_LP_dE/smdataExTsHP.mat')
K1 = 50 / sqrt(1250);
LAMBDA = 70.710678;
r = 0.02;
l= 9;
pL2ExLP = polyfit(p0smL2ExHP(1:l, 1), p0smL2ExHP(1:l, 2) *r, 1);
pM2ExLP = polyfit(p0smM2ExHP(1:l, 1), p0smM2ExHP(1:l, 2) *r, 1);
pExTsLP = polyfit(p0smExTsHP(1:l, 1), p0smExTsHP(1:l, 2) *r, 1);
pL2ExLP_slope = pL2ExLP(1);
pM2ExLP_slope = pM2ExLP(1);
pExTsLP_slope = pExTsLP(1);
%% LoTs, LBLo, BBLo with inverse perimeter modulus 
legend_font = 18;

K = 50 / sqrt(1250)^3; 
A0 = 1250;
shape_index = 0.02;
LAMBDA = K *A0 / shape_index;

r1 = 1;%K*A0 / LAMBDA;
load('LoTs_dE/smdataLoTs.mat')
errorbar(p0smLoTs(:,1), r1*p0smLoTs(:,2), r1*SMvarLoTs, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')

ylabel('$r \mathbf{\langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

K2 = 50 / sqrt(625);
%r2 = K2 / LAMBDA * sqrt(1250 / 625);
r2 = 1;%r1*(A0 / 625);
load('LBLo_dE/smdataLBLo.mat')
errorbar(p0smLBLo(:,1), r2*p0smLBLo(:,2), r2*SMvarLBLo, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

K3 = 50 / sqrt(2500);
%r3 = K3 / LAMBDA * sqrt(1250 / 2500);
r3 = 1;%r1*(A0 / 2500);
load('BBLo_dE/smdataBBLo.mat')
errorbar(p0smBBLo(:,1), r3*p0smBBLo(:,2), r3*SMvarBBLo, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_{avg} = A_0$','$A_{avg} = 0.5 A_0 $ ', '$A_{avg} = 2 A_0 $ ' , 'Interpreter', 'latex' ,'FontSize', legend_font)

savefig('LoTs_LBLo_BBLo_overlay_scaled_r.fig')
exportgraphics(gcf, 'LoTs_LBLo_BBLo_overlay_scaled_r.fig.pdf', 'ContentType', 'vector')
%% LoTs_LL, M2Lo_LL, L2Lo_LL
K = 50 / sqrt(1250);

load('L2Lo_LL_dE/smdataL2LoLL.mat')
errorbar(p0smL2LoLL(:,1), 1*p0smL2LoLL(:,2) / K, 1*SMvarL2LoLL / K, 'diamond', 'MarkerSize', 8, 'MarkerEdgeColor', 'black', 'Color', 'black')

ylabel('$\mathbf{r \langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;
%grid on;
load('M2Lo_LL_dE/smdataM2LoLL.mat')
errorbar(p0smM2LoLL(:,1), 1*p0smM2LoLL(:,2) / K, 1*SMvarM2LoLL / K, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('LoTs_LL_dE/smdataLoTsLL.mat')
errorbar(p0smLoTsLL(:,1), 1*p0smLoTsLL(:,2) / K, 1*SMvarLoTsLL / K, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

%title('Shear Modulus v Shape Index (Quadratic)');

legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', 18)

savefig('LoTs_M2Lo_L2Lo_LL_overlay.fig')
exportgraphics(gcf, 'LoTs_M2Lo_L2Lo_LL_overlay.pdf', 'ContentType', 'vector')
hold off;

%% ExTs_LL, M2Ex_LL, L2Ex_LL
K = 50 / sqrt(1250);

load('L2Ex_LL_dE/smdataL2ExLL.mat')
errorbar(p0smL2ExLL(:,1), 1*p0smL2ExLL(:,2) / K, 1*SMvarL2ExLL / K, 'diamond', 'MarkerSize', 8, 'MarkerEdgeColor', 'black', 'Color', 'black')

ylabel('$\mathbf{r \langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;
%grid on;
load('M2Ex_LL_dE/smdataM2ExLL.mat')
errorbar(p0smM2ExLL(:,1), 1*p0smM2ExLL(:,2) / K, 1*SMvarM2ExLL / K, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('ExTs_LL_dE/smdataExll.mat')
errorbar(p0smExll(:,1), 1*p0smExll(:,2) / K, 1*SMvarExll / K, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

%title('Shear Modulus v Shape Index (Quadratic)');

legend('$A_0 = 0.5A_{avg}$','$A_0 = 2A_{avg}$ ', '$A_0 = A_{avg}$' , 'Interpreter', 'latex' ,'FontSize', 18)

savefig('ExTs_M2Ex_L2Ex_LL_overlay.fig')
exportgraphics(gcf, 'ExTs_M2Ex_L2Ex_LL_overlay.pdf', 'ContentType', 'vector')
hold off;

%% MTQ3 and QuTs
K = 50 / sqrt(1250);

load('QuTs_dE/smdataQuTs.mat')
errorbar(p0smQuTs(:,1), p0smQuTs(:,2) / K, SMvarQuTs / K, 'diamond', 'MarkerSize', 8, 'MarkerEdgeColor', 'black', 'Color', 'black')

ylabel('$\mathbf{r \langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on; 
load('MTQ3_dE/smdataMTQ3.mat');
errorbar(p0smMTQ3(:,1), p0smMTQ3(:,2) / K, SMvarMTQ3 / K, 'square', 'MarkerSize', 8, 'MarkerEdgeColor', 'blue', 'Color', 'blue')
%% QuTw, Qu2w, QuOw
K = 50 / sqrt(1250);

load('QuTw_dE/smdataQuTw.mat');
errorbar(p0smQuTw(:,1), p0smQuTw(:,2) / K, SMvarQuTw / K, 'diamond', 'markerSize', 8, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylim([-0.1, 3])
ylabel('$\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('Qu2w_dE/smdataQu2w.mat');
errorbar(p0smQu2w(:,1), p0smQu2w(:,2) / K, SMvarQu2w / K, 'square', 'MarkerSize', 8, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('QuOw_dE/smdataQuOw.mat');
errorbar(p0smQuOw(:,1), p0smQuOw(:,2) / K, SMvarQuOw / K, 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'red', 'Color', 'red')

load('T2ss_dE/smdata7.mat')
errorbar(p0sm(:,1), 1*p0sm(:,2)/ K, 1*SMvar/ K, 'pentagram', 'MarkerSize', 8, 'MarkerEdgeColor', [0.07 0.81 0.16], 'Color', [0.07 0.81 0.16])

yline(0)
title('Quartic Sweep, $C_4 = 1,2, 3$', 'Interpreter', 'latex', 'FontSize', 22)
legend('$C_4 = 3$', '$C_4 = 2$', '$C_4 = 1$', '$C_4 = 0$', 'Interpreter', 'latex', 'FontSize', 18)
savefig('T2ss_QuTw_Qu2w_QuOw_overlay_C4_0_1_2_3.fig')
exportgraphics(gcf, 'T2ss_QuTw_Qu2w_QuOw_overlay_C4_0_1_2_3.pdf', 'ContentType', 'vector')
%% Qu2w, QuOw
K = 50 / sqrt(1250);

load('Qu2w_dE/smdataQu2w.mat');
errorbar(p0smQu2w(:,1), p0smQu2w(:,2) / K, SMvarQu2w / K, 'diamond', 'markerSize', 8, 'MarkerEdgeColor', 'blue', 'Color', 'blue')
ylim([-0.1, 3])
ylabel('$\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('QuOw_dE/smdataQuOw.mat');
errorbar(p0smQuOw(:,1), p0smQuOw(:,2) / K, SMvarQuOw / K, 'square', 'MarkerSize', 8, 'MarkerEdgeColor', 'red', 'Color', 'red')

load('T2ss_dE/smdata7.mat')
errorbar(p0sm(:,1), 1*p0sm(:,2)/ K, 1*SMvar/ K, 'pentagram', 'MarkerSize', 6, 'MarkerEdgeColor', [0.07 0.81 0.16], 'Color', [0.07 0.81 0.16])

yline(0)
title('Quartic Sweep, $C_4 = 0,1,2$', 'Interpreter', 'latex', 'FontSize', 22)
legend('$C_4 = 2$', '$C_4 = 1$', '$C_4 = 0$', 'Interpreter', 'latex', 'FontSize', 18)
savefig('Qu2w_QuOw_overlay_C4_0_1_2.fig')
exportgraphics(gcf, 'Qu2w_QuOw_overlay_C4_0_1_2.pdf', 'ContentType', 'vector')

%% QuTw and T2ss

K = 50 / sqrt(1250);

load('QuTw_dE/smdataQuTw.mat');
errorbar(p0smQuTw(:,1), p0smQuTw(:,2) / K, SMvarQuTw / K, 'diamond', 'markerSize', 8, 'MarkerEdgeColor', 'black', 'Color', 'black')
hold on; 

load('T2ss_dE/smdata7.mat')
errorbar(p0sm(:,1), 1*p0sm(:,2)/ K, 1*SMvar/ K, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', [0.07 0.81 0.16], 'Color', [0.07 0.81 0.16])

yline(0)
title('Quartic Sweep, $C_4 = 0, 3$', 'Interpreter', 'latex', 'FontSize', 22)
legend('$C_4 = 3$', '$C_4 = 0$', 'Interpreter', 'latex', 'FontSize', 18)
savefig('QuTw_T2ss_overlay_C4_0_3.fig')
exportgraphics(gcf, 'QuTw_T2ss_overlay_C4_0_3.pdf', 'ContentType', 'vector')

%% QuOw, AqOw, LqOw
K = 50 / sqrt(1250);

load('QuOw_dE/smdataQuOw.mat');
errorbar(p0smQuOw(:,1), p0smQuOw(:,2) / K, SMvarQuOw / K, 'diamond', 'markerSize', 8, 'MarkerEdgeColor', 'blue', 'Color', 'blue')
ylim([-0.1, 3])
ylabel('$\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('AqOw_dE/smdataAqOw.mat');
errorbar(p0smAqOw(:,1), p0smAqOw(:,2) / K, SMvarAqOw / K, 'square', 'MarkerSize', 8, 'MarkerEdgeColor', 'red', 'Color', 'red')

load('LqOw_dE/smdataLqOw.mat')
errorbar(p0smLqOw(:,1), 1*p0smLqOw(:,2)/ K, 1*SMvarLqOw/ K, 'pentagram', 'MarkerSize', 6, 'MarkerEdgeColor', [0.07 0.81 0.16], 'Color', [0.07 0.81 0.16])

yline(0)
title('Quartic Sweep, $C_4 = 1$, A0 sweep', 'Interpreter', 'latex', 'FontSize', 22)
legend('$A_0 = 2 A_{avg}$', '$A_0 = A_{avg}$', '$A_0 = 0.5 A_{avg}$', 'Interpreter', 'latex', 'FontSize', 18)
savefig('QuOw_AqOw_LqOw_overlay_A0_1_1_0.5_Ab.fig')
exportgraphics(gcf, 'QuOw_AqOw_LqOw_overlay_A0_2_1_0.5_Ab.pdf', 'ContentType', 'vector')

%% Qu2w, Aq2w, Lq2w
K = 50 / sqrt(1250);

load('Qu2w_dE/smdataQu2w.mat');
errorbar(p0smQu2w(:,1), p0smQu2w(:,2) / K, SMvarQu2w / K, 'diamond', 'markerSize', 8, 'MarkerEdgeColor', 'blue', 'Color', 'blue')
ylim([-0.1, 3])
ylabel('$\mathbf{\langle g \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{A_{avg} }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('Aq2w_dE/smdataAq2w.mat');
errorbar(p0smAq2w(:,1), p0smAq2w(:,2) / K, SMvarAq2w / K, 'square', 'MarkerSize', 8, 'MarkerEdgeColor', 'red', 'Color', 'red')

load('Lq2w_dE/smdataLq2w.mat')
errorbar(p0smLq2w(:,1), 1*p0smLq2w(:,2)/ K, 1*SMvarLq2w/ K, 'pentagram', 'MarkerSize', 6, 'MarkerEdgeColor', [0.07 0.81 0.16], 'Color', [0.07 0.81 0.16])

yline(0)
title('Quartic Sweep, $C_4 = 2$, A0 sweep', 'Interpreter', 'latex', 'FontSize', 22)
legend('$A_0 = 2 A_{avg}$', '$A_0 = A_{avg}$', '$A_0 = 0.5 A_{avg}$', 'Interpreter', 'latex', 'FontSize', 18)
savefig('Qu2w_Aq2w_Lq2w_overlay_A0_2_1_0.5_Ab.fig')
exportgraphics(gcf, 'Qu2w_Aq2w_Lq2w_overlay_A0_2_1_0.5_Ab.pdf', 'ContentType', 'vector')