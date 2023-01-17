
%% This creates the plot that overlays the low polydispersity quadratic model shear moduli
load('L2Ts_dE/smdataL2Ts.mat')
errorbar(p0smL2Ts(:,1), 0.02*p0smL2Ts(:,2), 0.02*SMvarL2Ts, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')

ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;
%grid on;
load('M2Ts_dE/smdataM2Ts.mat')
errorbar(p0smM2Ts(:,1), 0.02*p0smM2Ts(:,2), 0.02*SMvarM2Ts, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('T2ss_dE/smdata7.mat')
errorbar(p0sm(:,1), 0.02*p0sm(:,2), 0.02*SMvar, 'diamond', 'MarkerSize', 8, 'MarkerEdgeColor', 'red', 'Color', 'red')
%title('Shear Modulus v Shape Index (Quadratic)');

legend('$A_0 = 0.5\langle A\rangle$','$A_0 = 2\langle A\rangle$ ', '$A_0 = \langle A\rangle$' , 'Interpreter', 'latex' ,'FontSize', 18)

savefig('L2Ts_M2Ts_T2ss_overlay.fig')
exportgraphics(gcf, 'L2Ts_M2Ts_T2ss_overlay.pdf', 'ContentType', 'vector')
hold off;
%% This creates the plot that overlays the lower polydispersity logarithmic model shear moduli
load('L2Lo_dE/smdataL2Lo.mat')

errorbar(p0smL2Lo(:,1), 0.02*p0smL2Lo(:,2), 0.02*SMvarL2Lo, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;
%grid on;

load('M2Lo_dE/smdataM2Lo.mat')
errorbar(p0smM2Lo(:,1), 0.02*p0smM2Lo(:,2), 0.02*SMvarM2Lo, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('dE_LoTs/smdataLoTs.mat')
errorbar(p0smLoTs(:,1), 0.02*p0smLoTs(:,2), 0.02*SMvarLoTs, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5\langle A\rangle$','$A_0 = 2\langle A\rangle$ ', '$A_0 = \langle A\rangle$' , 'Interpreter', 'latex' ,'FontSize', 18)

savefig('L2Lo_M2Lo_LoTs_overlay.fig')
exportgraphics(gcf, 'L2Lo_M2Lo_LoTs_overlay.pdf', 'ContentType', 'vector')
hold off;
%% This creates the plot that overlays the low polydispersity exponential model shear moduli

load('L2Ex_dE/smdataL2Ex.mat')

errorbar(p0smL2Ex(:,1), 0.02*p0smL2Ex(:,2), 0.02*SMvarL2Ex, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('M2Ex_dE/smdataM2Ex.mat')
errorbar(p0smM2Ex(:,1), 0.02*p0smM2Ex(:,2), 0.02*SMvarM2Ex, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('ExTs_dE/smdataExTs.mat')
errorbar(p0smExTs(:,1), 0.02*p0smExTs(:,2), 0.02*SMvarExTs, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5\langle A\rangle$','$A_0 = 2\langle A\rangle$ ', '$A_0 = \langle A\rangle$' , 'Interpreter', 'latex' ,'FontSize', 18)

savefig('L2Ex_M2Ex_ExTs_overlay.fig')
exportgraphics(gcf, 'L2Ex_M2Ex_ExTs_overlay.pdf', 'ContentType', 'vector')
hold off;

%% This creates the plot that overlays the low polydispersity shear moduli plots where all three plots have A0 = <A>

load('T2ss_dE/smdata7.mat')
errorbar(p0sm(:,1), 0.02*p0sm(:,2), 0.02*SMvar, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')

ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('dE_LoTs/smdataLoTs.mat')
errorbar(p0smLoTs(:,1), 0.02*p0smLoTs(:,2), 0.02*SMvarLoTs, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('ExTs_dE/smdataExTs.mat')
errorbar(p0smExTs(:,1), 0.02*p0smExTs(:,2), 0.02*SMvarExTs, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('Quadratic $E^A$','Logarithmic $E^A$ ', 'Exponential $E^A$' , 'Interpreter', 'latex' ,'FontSize', 18)
savefig('T2ss_LoTs_ExTs_overlay.fig')
exportgraphics(gcf, 'T2ss_LoTs_ExTs_overlay.pdf', 'ContentType', 'vector')

%% This creates the plot that overlays the high polydiserpsity quadratic shear moduli plot

% load('L2Ts_LP_dE/smdataL2Ts.mat') -- I never properly made this up

load('M2Ts_LP_dE/smdataM2TsHP.mat')
errorbar(p0smM2TsHP(:,1), 0.02*p0smM2TsHP(:,2), 0.02*SMvarM2TsHP, 'o')
ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('T2ss_LP_dE/smdataTshp.mat')
errorbar(p0smTsHP(:,1), 0.02*p0smTsHP(:,2), 0.02*SMvarTsHP, 'o')

% RERUN L2TS_LP_DE AND GET THAT FIXED UP

%%  This creates the plot that overlays the high polydispersity logarithmic shear moduli plot

load('L2Lo_LP_dE/smdataL2LoHP.mat')
errorbar(p0smL2LoHP(:,1), 0.02*p0smL2LoHP(:,2), 0.02*SMvarL2LoHP, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('M2Lo_LP_dE/smdataM2LoHP.mat')
errorbar(p0smM2LoHP(:,1), 0.02*p0smM2LoHP(:,2), 0.02*SMvarM2LoHP, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('LoTs_LP_dE/smdataLoTsHP.mat')
errorbar(p0smLoTsHP(:,1), 0.02*p0smLoTsHP(:,2), 0.02*SMvarLoTsHP, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5\langle A\rangle$','$A_0 = 2\langle A\rangle$ ', '$A_0 = \langle A\rangle$' , 'Interpreter', 'latex' ,'FontSize', 18)

savefig('L2Lo_M2Lo_LoTs_HighPoly_overlay.fig')
exportgraphics(gcf, 'L2Lo_M2Lo_LoTs_HighPoly_overlay.pdf', 'ContentType', 'vector')


%% This creates the plot that overlays the high polydispsersity exponential shear moduli plot

load('L2Ex_LP_dE/smdataL2ExHP.mat')
errorbar(p0smL2ExHP(:,1), 0.02*p0smL2ExHP(:,2), 0.02*SMvarL2ExHP, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName','Times New Roman');
hold on;

load('M2Ex_LP_dE/smdataM2ExHP.mat')
errorbar(p0smM2ExHP(:,1), 0.02*p0smM2ExHP(:,2), 0.02*SMvarM2ExHP, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('ExTs_LP_dE/smdataExTsHP.mat')
errorbar(p0smExTsHP(:,1), 0.02*p0smExTsHP(:,2), 0.02*SMvarExTsHP, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5\langle A\rangle$','$A_0 = 2\langle A\rangle$ ', '$A_0 = \langle A\rangle$' , 'Interpreter', 'latex' ,'FontSize', 18)


savefig('L2Ex_M2Ex_ExTs_HighPoly_overlay.fig')
exportgraphics(gcf, 'L2Ex_M2Ex_ExTs_HighPoly_overlay.pdf', 'ContentType', 'vector')

%% This creates the plot that overlays the high polydispersity shear moduli plots where all three plots have A0 = <A>

%load('T2ss_LP_dE/smdataTshp.mat')
%errorbar(p0smTsHP(:,1), 0.02*p0smTsHP(:,2), 0.02*SMvarTsHP, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')

%ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
%xlabel('$\mathbf{\frac{P_0}{\langle A \rangle }}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
%hold on;

load('LoTs_LP_dE/smdataLoTsHP.mat')
errorbar(p0smLoTsHP(:,1), 0.02*p0smLoTsHP(:,2), 0.02*SMvarLoTsHP, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')
ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
hold on;

load('ExTs_LP_dE/smdataExTsHP.mat')
errorbar(p0smExTsHP(:,1), 0.02*p0smExTsHP(:,2), 0.02*SMvarExTsHP, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('Quadratic $E^A$','Logarithmic $E^A$ ', 'Exponential $E^A$' , 'Interpreter', 'latex' ,'FontSize', 18)
savefig('T2ss_LoTs_ExTs_HighPoly_overlay.fig')
exportgraphics(gcf, 'T2ss_LoTs_ExTs_overlay_HighPoly_overlay.pdf', 'ContentType', 'vector')

%% This combines the plots that show A0=<A> for high and low polydispersity cases

subgplotlabelindex_x = 3.73;
subgplotlabelindex_y = .105;
plot_width = 1200;
plot_height = 360;

subplot(1,2,1); %low polydispersity plot first 

load('T2ss_dE/smdata7.mat')
errorbar(p0sm(:,1), 0.02*p0sm(:,2), 0.02*SMvar, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'black', 'Color', 'black')
%yticks('manual')

hF = gcf;
hF.Position(1:4) = [0 0 plot_width plot_height];

ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', 18, 'FontName','Times New Roman');
hold on;
grid on;
load('dE_LoTs/smdataLoTs.mat')
errorbar(p0smLoTs(:,1), 0.02*p0smLoTs(:,2), 0.02*SMvarLoTs, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('ExTs_dE/smdataExTs.mat')
errorbar(p0smExTs(:,1), 0.02*p0smExTs(:,2), 0.02*SMvarExTs, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')
ylim([-0.01 0.11])
yticks(-0.01:0.01:0.11)
yticklabels({'-0.01', '0', '0.01', '0.02', '0.03', '0.04', '0.05','0.06', '0.07', '0.08', '0.09', '0.10', '0.11'})
xticks(3.75:0.03:3.9)
legend('Quadratic $E^A$','Logarithmic $E^A$ ', 'Exponential $E^A$' , 'Interpreter', 'latex' ,'FontSize', 18)
title('$\textbf{Low Polydispersity Shear Moduli}$, $\mathbf{A_0 = \langle A\rangle}$', 'Interpreter', 'latex' ,'FontSize', 16, 'FontName','Times New Roman');
text(subgplotlabelindex_x ,subgplotlabelindex_y , '$\textbf{a}$', 'FontSize', 20,'Interpreter', 'latex');
%------------------------------------------------------------------
subplot(1,2,2); %then the high polydispersity case

load('LoTs_LP_dE/smdataLoTsHP.mat')

errorbar(p0smLoTsHP(:,1), 0.02*p0smLoTsHP(:,2), 0.02*SMvarLoTsHP, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')
%yticks(-0.01:0.01:0.12)
ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', 18, 'FontName','Times New Roman');
hold on;
grid on;
load('ExTs_LP_dE/smdataExTsHP.mat')
errorbar(p0smExTsHP(:,1), 0.02*p0smExTsHP(:,2), 0.02*SMvarExTsHP, 'diamond', 'MarkerSize', 6, 'MarkerEdgeColor', 'red', 'Color', 'red')
ylim([-0.01 0.11])
yticks(-0.01:0.01:0.11)
yticklabels({'-0.01', '0', '0.01', '0.02', '0.03', '0.04', '0.05','0.06', '0.07', '0.08', '0.09', '0.10', '0.11'})
xticks(3.75:0.03:3.9)
legend('Logarithmic $E^A$ ', 'Exponential $E^A$' , 'Interpreter', 'latex' ,'FontSize', 18)
title('$\textbf{High Polydispersity Shear Moduli}$, $\mathbf{A_0 = \langle A\rangle}$', 'Interpreter', 'latex' ,'FontSize', 16, 'FontName','Times New Roman');
text(subgplotlabelindex_x, subgplotlabelindex_y, '$\textbf{b}$', 'FontSize', 20,'Interpreter', 'latex');

savefig('T2ss_LoTs_ExTs_overlay_HP_AND_LP.fig')
exportgraphics(gcf, 'T2ss_LoTs_ExTs_overlay_HP_AND_LP.pdf', 'ContentType', 'vector')



%% This will combine the 6 plots that show the three different types of potentials for the three different values of A0 for high and low polydispersity 
%------------------------------------------------------------------
subgplotlabelindex_x = 3.73;
subgplotlabelindex_y = .105;
subgplotlabelfontsize = 20;
plot_width = 1200;
plot_height = 1800;
horizontal_title_x = 3.715;
horizontal_title_y = 0.015;
title_font = 18;
legend_font = 18;
markersize = 6;
xlabelfont = 16;
ylabelfont = 14;

subplot(3,2,1) % first low polydispersity Tiss case

hF = gcf;
hF.Position(1:4) = [0 0 plot_width plot_height];

load('L2Ts_dE/smdataL2Ts.mat')
errorbar(p0smL2Ts(:,1), 0.02*p0smL2Ts(:,2), 0.02*SMvarL2Ts, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'black', 'Color', 'black')

ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', ylabelfont, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', xlabelfont, 'FontName','Times New Roman');

hold on;
grid on;
load('M2Ts_dE/smdataM2Ts.mat')
errorbar(p0smM2Ts(:,1), 0.02*p0smM2Ts(:,2), 0.02*SMvarM2Ts, 'square', 'MarkerSize', markersize, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('T2ss_dE/smdata7.mat')
errorbar(p0sm(:,1), 0.02*p0sm(:,2), 0.02*SMvar, 'diamond', 'MarkerSize', markersize, 'MarkerEdgeColor', 'red', 'Color', 'red')
%title('Shear Modulus v Shape Index (Quadratic)');

legend('$A_0 = 0.5\langle A\rangle$','$A_0 = 2\langle A\rangle$ ', '$A_0 = \langle A\rangle$' , 'Interpreter', 'latex' ,'FontSize', legend_font)
ylim([-0.01 0.11])
yticks(-0.01:0.01:0.11)
yticklabels({'-0.01', '0', '0.01', '0.02', '0.03', '0.04', '0.05','0.06', '0.07', '0.08', '0.09', '0.10', '0.11'})
xticks(3.75:0.03:3.9)

title('$\textbf{Low Polydispersity Shear Moduli}$', 'Interpreter', 'latex' ,'FontSize', title_font, 'FontName','Times New Roman');
%subtitle('$\textbf{Quadratic } \mathbf{E^A}$', 'Interpreter', 'latex' ,'FontSize', 14, 'FontName','Times New Roman');
text(horizontal_title_x, horizontal_title_y+0.01, '$\textbf{Quadratic}$', 'FontSize', title_font,'Interpreter', 'latex', 'Rotation', 90);
text(subgplotlabelindex_x, subgplotlabelindex_y, '$\textbf{a}$', 'FontSize', subgplotlabelfontsize,'Interpreter', 'latex');

%------------------------------------------------------------------
subplot(3,2,3)  % second, low polydisperisty LoTs case

load('L2Lo_dE/smdataL2Lo.mat')

errorbar(p0smL2Lo(:,1), 0.02*p0smL2Lo(:,2), 0.02*SMvarL2Lo, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'black', 'Color', 'black')
%ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', xlabelfont, 'FontName','Times New Roman');
ylabel(' $\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', ylabelfont, 'FontName', 'Times New Roman');
%subtitle('$\textbf{Logarithmic } \mathbf{E^A}$', 'Interpreter', 'latex' ,'FontSize', 14, 'FontName','Times New Roman');
hold on;
grid on;

load('M2Lo_dE/smdataM2Lo.mat')
errorbar(p0smM2Lo(:,1), 0.02*p0smM2Lo(:,2), 0.02*SMvarM2Lo, 'square', 'MarkerSize', markersize, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('dE_LoTs/smdataLoTs.mat')
errorbar(p0smLoTs(:,1), 0.02*p0smLoTs(:,2), 0.02*SMvarLoTs, 'diamond', 'MarkerSize', markersize, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5\langle A\rangle$','$A_0 = 2\langle A\rangle$ ', '$A_0 = \langle A\rangle$' , 'Interpreter', 'latex' ,'FontSize', legend_font)
ylim([-0.01 0.11])
yticks(-0.01:0.01:0.11)
yticklabels({'-0.01', '0', '0.01', '0.02', '0.03', '0.04', '0.05','0.06', '0.07', '0.08', '0.09', '0.10', '0.11'})
xticks(3.75:0.03:3.9)
text(horizontal_title_x, horizontal_title_y, '$\textbf{Logarithmic}$', 'FontSize', title_font,'Interpreter', 'latex', 'Rotation', 90);
text(subgplotlabelindex_x, subgplotlabelindex_y, '$\textbf{c}$', 'FontSize', subgplotlabelfontsize,'Interpreter', 'latex');
%------------------------------------------------------------------
subplot(3,2,5)  % third, low polydisperisty ExTs case


load('L2Ex_dE/smdataL2Ex.mat')

errorbar(p0smL2Ex(:,1), 0.02*p0smL2Ex(:,2), 0.02*SMvarL2Ex, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', ylabelfont, 'FontName', 'Times New Roman');
%ylabel({'$\textbf{Exponential}$ ',' $\mathbf{r \langle G \rangle}$'}, 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', xlabelfont, 'FontName','Times New Roman');
%subtitle('$\textbf{Exponential } \mathbf{E^A}$', 'Interpreter', 'latex' ,'FontSize', 14, 'FontName','Times New Roman');
hold on;
grid on;
load('M2Ex_dE/smdataM2Ex.mat')
errorbar(p0smM2Ex(:,1), 0.02*p0smM2Ex(:,2), 0.02*SMvarM2Ex, 'square', 'MarkerSize', markersize, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('ExTs_dE/smdataExTs.mat')
errorbar(p0smExTs(:,1), 0.02*p0smExTs(:,2), 0.02*SMvarExTs, 'diamond', 'MarkerSize', markersize, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5\langle A\rangle$','$A_0 = 2\langle A\rangle$ ', '$A_0 = \langle A\rangle$' , 'Interpreter', 'latex' ,'FontSize', legend_font)

ylim([-0.01 0.11])
yticks(-0.01:0.01:0.11)
yticklabels({'-0.01', '0', '0.01', '0.02', '0.03', '0.04', '0.05','0.06', '0.07', '0.08', '0.09', '0.10', '0.11'})
xticks(3.75:0.03:3.9)
text(horizontal_title_x, horizontal_title_y, '$\textbf{Exponential}$', 'FontSize', title_font,'Interpreter', 'latex', 'Rotation', 90);
text(subgplotlabelindex_x, subgplotlabelindex_y, '$\textbf{e}$', 'FontSize', subgplotlabelfontsize,'Interpreter', 'latex');
%------------------------------------------------------------------
subplot(3,2,2); % now going to the second column with the high polydispersity T2ss plot


load('M2Ts_LP_dE/smdataM2TsHP.mat')
errorbar(p0smM2TsHP(:,1), 0.02*p0smM2TsHP(:,2), 0.02*SMvarM2TsHP, 'square', 'MarkerSize', 6, 'MarkerEdgeColor', 'blue', 'Color', 'blue')
ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', 14, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontName','Times New Roman');
title('$\textbf{High Polydispersity Shear Moduli}$', 'Interpreter', 'latex' ,'FontSize', 18, 'FontName','Times New Roman');
%subtitle('$\textbf{Quadratic } \mathbf{E^A}$', 'Interpreter', 'latex' ,'FontSize', 14, 'FontName','Times New Roman');
legend('$A_0 = 2\langle A\rangle$ ' , 'Interpreter', 'latex' ,'FontSize', 18)
ylim([-0.01 0.11])
yticks(-0.01:0.01:0.11)
yticklabels({'-0.01', '0', '0.01', '0.02', '0.03', '0.04', '0.05','0.06', '0.07', '0.08', '0.09', '0.10', '0.11'})
xticks(3.75:0.03:3.9);
text(subgplotlabelindex_x, subgplotlabelindex_y, '$\textbf{b}$', 'FontSize', 20,'Interpreter', 'latex');
grid on;

%-------------------------------------------------------------------
subplot(3,2,4); % this is the high polydispiersity log plot

load('L2Lo_LP_dE/smdataL2LoHP.mat')
errorbar(p0smL2LoHP(:,1), 0.02*p0smL2LoHP(:,2), 0.02*SMvarL2LoHP, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel('$\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', ylabelfont, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', xlabelfont, 'FontName','Times New Roman');
%subtitle('$\textbf{Logarithmic } \mathbf{E^A}$', 'Interpreter', 'latex' ,'FontSize', 14, 'FontName','Times New Roman');
hold on;
grid on;
load('M2Lo_LP_dE/smdataM2LoHP.mat')
errorbar(p0smM2LoHP(:,1), 0.02*p0smM2LoHP(:,2), 0.02*SMvarM2LoHP, 'square', 'MarkerSize', markersize, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('LoTs_LP_dE/smdataLoTsHP.mat')
errorbar(p0smLoTsHP(:,1), 0.02*p0smLoTsHP(:,2), 0.02*SMvarLoTsHP, 'diamond', 'MarkerSize', markersize, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5\langle A\rangle$','$A_0 = 2\langle A\rangle$ ', '$A_0 = \langle A\rangle$' , 'Interpreter', 'latex' ,'FontSize', legend_font)
ylim([-0.01 0.11])
yticks(-0.01:0.01:0.11)
yticklabels({'-0.01', '0', '0.01', '0.02', '0.03', '0.04', '0.05','0.06', '0.07', '0.08', '0.09', '0.10', '0.11'})
xticks(3.75:0.03:3.9);
text(subgplotlabelindex_x, subgplotlabelindex_y, '$\textbf{d}$', 'FontSize', subgplotlabelfontsize,'Interpreter', 'latex');
%-------------------------------------------------------------------
subplot(3,2,6); %high polydispiersity expoential plot

load('L2Ex_LP_dE/smdataL2ExHP.mat')
errorbar(p0smL2ExHP(:,1), 0.02*p0smL2ExHP(:,2), 0.02*SMvarL2ExHP, 'o', 'MarkerSize', markersize, 'MarkerEdgeColor', 'black', 'Color', 'black')
ylabel(' $\mathbf{r \langle G \rangle}$', 'Interpreter', 'latex', 'FontSize', ylabelfont, 'FontName', 'Times New Roman');
xlabel('$ \mathbf{\frac{P_0}{\sqrt{\langle A \rangle }}}$', 'Interpreter', 'latex', 'FontSize', xlabelfont, 'FontName','Times New Roman');
%subtitle('$\textbf{Exponential } \mathbf{E^A}$', 'Interpreter', 'latex' ,'FontSize', 14, 'FontName','Times New Roman');
hold on;
grid on;
load('M2Ex_LP_dE/smdataM2ExHP.mat')
errorbar(p0smM2ExHP(:,1), 0.02*p0smM2ExHP(:,2), 0.02*SMvarM2ExHP, 'square', 'MarkerSize', markersize, 'MarkerEdgeColor', 'blue', 'Color', 'blue')

load('ExTs_LP_dE/smdataExTsHP.mat')
errorbar(p0smExTsHP(:,1), 0.02*p0smExTsHP(:,2), 0.02*SMvarExTsHP, 'diamond', 'MarkerSize', markersize, 'MarkerEdgeColor', 'red', 'Color', 'red')

legend('$A_0 = 0.5\langle A\rangle$','$A_0 = 2\langle A\rangle$ ', '$A_0 = \langle A\rangle$' , 'Interpreter', 'latex' ,'FontSize', legend_font)
ylim([-0.01 0.11])
yticks(-0.01:0.01:0.11)
yticklabels({'-0.01', '0', '0.01', '0.02', '0.03', '0.04', '0.05','0.06', '0.07', '0.08', '0.09', '0.10', '0.11'})
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
exportgraphics(gcf, 'polydispersity_distribution.pdf', 'ContentType', 'vector')