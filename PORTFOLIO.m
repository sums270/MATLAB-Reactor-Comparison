%% COMPARE MODELLING OF BATCH, CSTR AND PFR REACTOR
% Hydrogen peroxide decomposition
% Reaction: 2H2O2 -> 2H2O + O2
% Kinetic model: -rA = k*CA
% Assumptions: 
% Isothermal operation
% Ireversable first-order reaction
% Ideal reactor behaviour
% Constant liquid voloume/density
% Gas evolutionn and mass transfer effects neglected
clear;
clc;
close all;
%% Define portfolio colors
burgandy = [104 21 64]/255;
cream = [247 243 240]/255;
charcoal = [36 33 36]/255;
pink = [235 145 179]/255;
gold = [235 191 111]/255;
%% Define reaction parameters
k = 0.0648; %(first order rate constant)
CA0 = 0.3;  %(initial concentration)
T = 20;  %(temperature)
% common caracteristic time range
time = linspace(0,80,500);
% for batch, time represents reaction time, 
% for CSTR/PFR, time represents space time
tau = time;
%% Batch reactor model
% Equation = dCA/dt = -k*CA
% intergrated equation = CA = CA0*exp(-k*t)
% conversion = X = 1 - CA/CA0
CA_batch = CA0*exp(-k*time);
X_batch = 1 - CA_batch/CA0;
%% CSTR model (STEADY STATE)
% MASS BALANCE- F*CA0 - F*CA - V*K*CA = 0 
% SPACE TIME - tau = V/F
X_cstr = (k.*tau)./(1+k.*tau);
CA_cstr = CA0.*(1-X_cstr);
%% PFR MODEL (STEADY STATE)
% Equation = dCA/dtau = -k*CA
% Intergrated equastion = CA= CA0*exp(-k*tau)
% Convertion X = 1-exp(-k*tau)
CA_pfr = CA0.*exp(-k.*tau);
X_pfr = 1-CA_pfr/CA0;
%% CONVERSION COMPARISON GRAPH
fig1 = figure('color', charcoal,...
    'Position',[100 100 1150 680]);
ax = axes(fig1);
hold(ax,'on');
%plot batch
plot(time,X_batch*100,'-',...
    'Color',pink,'LineWidth',3);
%plot CSTR
plot(tau,X_cstr*100,'-',...
    'Color',gold,'Linewidth',3);
%plot PFR (added dashed line due to overlap)
plot(tau,X_pfr*100,'--',...
    'Color',cream,'LineWidth',2);
% FORMAT GRAPH
ax.Color = charcoal;
ax.XColor = cream;
ax.YColor = cream;
ax.GridColor = cream;
ax.GridAlpha = 0.12;
ax.FontSize = 12;
ax.LineWidth = 1;
grid on;
box on;
xlabel('Reaction Time / Space Time (h)',...
    'Color',cream,'FontSize',13);
ylabel('Conversion (%)',...
    'Color',cream,'FontSize',13);
title('BATCH vs CSTR vs PFR',...
    'Color',cream,...
    'FontSize',20,...
    'FontWeight','bold');
subtitle('Hydrogen Peroxide Decomposition | 20 celcius',...
    'Color',pink);
ylim([0 105]);
xlim([0 80]);
lgd = legend('Batch','CSTR','PFR',...
    'Location','southeast');
lgd.TextColor = cream;
lgd.Color = charcoal;
lgd.EdgeColor = pink;
%% ADD REACTION AND EQUASION TO GRAPH
reactionText = {
    '$2H_2O_2 \ rightarrow 2H_2O + O_2$'
    '$-rA = kC_A$'
    '$k = 0.0648\ \mathrm{h}^{-1}$'
    };
annotation(fig1,'textbox',[0.18 0.68 0.35 0.16],...
    'String',reactionText,...
    'Interpreter','latex',...
    'Color',cream,...
    'BackgroundColor', burgandy,...
    'EdgeColor',pink,...
    'LineWidth',1.2,...
    'Margin',10,...
    'FontSize',14);
%export graph
exportgraphics(fig1,'reactor_conversion_comparison.png','Resolution',300);
%% CALCULATE CONVERSION AT 10 HOURS
t_check = 10;
X_batch_10 = 1-exp(-k*t_check);
X_cstr_10 = (k*t_check)/(1 + k*t_check);
X_pfr_10 = 1 - exp(-k*t_check);
%% CALCULATE TIME FOR 80% CONVERSION
X_target = 0.80;
%batch reaction time 
t_batch_80 = -log(1-X_target)/k;
%CSTR space time 
tau_cstr_80 = X_target/(k*(1-X_target));
%PFR space time
tau_pfr_80 = -log(1-X_target)/k;
%% DISPLAY RESULTS
fprintf('\n==========\n');
fprintf(' REACTOR COMPARISON RESULTS\n ');
fprintf('==========\n');
fprintf('Reaction: Hydrogen Peroxide decomposition\n');
fprintf('Temperature: %.1f C/n',k);
fprintf('Rate constant: %.4f h^-1\n', k);
fprintf('Initial concentration: %.2f mol/L\n',CA0);
fprintf('\nConversion at 10 hours:\n');
fprintf('Batch: %.2f%% | CSTR: %.2f%% | PFR: %.2f%%\n', ...
    100*[X_batch_10, X_cstr_10, X_pfr_10]);
fprintf('\nTime/space time for 80%% conversion:\n');
fprintf('Batch: %.2f h | CSTR: %.2f h | PFR: %.2f h\n', ...
    t_batch_80, tau_cstr_80, tau_pfr_80);
%% TARGET CONVERSION BAR CHART 
fig2 = figure('Color',charcoal,...
    'Position',[150 150 950 600]);
values = [t_batch_80 tau_cstr_80 tau_pfr_80];
b = bar(values,0.55,'FaceColor','flat');
b.CData = [pink;gold;cream;];
ax2 = gca;
ax2.Color = charcoal;
ax2.XColor = cream;
ax2.YColor = cream;
ax2.FontSize = 13;
ax2.GridColor = cream;
ax2.GridAlpha = 0.12;
xticks(1:3);
xticklabels({'Batch','CSTR','PFR'});
ylabel('Time / Residence Time (h)', ...
    'Color',cream);
title('TIME REQUIRED FOR 80% CONVERSION', ...
    'Color',cream,'FontSize',18);
grid on;
for i = 1:3
    text(i,values(i)+2, ...
        sprintf('%.2f h',values(i)), ...
        'HorizontalAlignment','center', ...
        'Color',cream,'FontSize',13);
end
ylim([0 max(values)*1.2]);
exportgraphics(fig2, ...
    'Reactor_80pct_Comparison.png', ...
    'Resolution',300);
%% REACTOR EQUATIONS FIGURE
fig3 = figure('Color',burgandy, ...
    'Position',[100 100 1100 700]);
axis off;
text(0.07,0.91,'REACTOR DESIGN EQUATIONS', ...
    'Color',cream, ...
    'FontSize',23, ...
    'FontWeight','bold');
text(0.07,0.80, ...
    '2H_2O_2  ->  2H_2O + O_2', ...
    'Interpreter','tex', ...
    'Color',pink, ...
    'FontSize',18);
text(0.07,0.68,'BATCH REACTOR', ...
    'Color',pink, ...
    'FontSize',16, ...
    'FontWeight','bold');
text(0.07,0.61,'X = 1 - e^{-kt}', ...
    'Interpreter','tex', ...
    'Color',cream, ...
    'FontSize',20);
text(0.07,0.48,'CSTR', ...
    'Color',pink, ...
    'FontSize',16, ...
    'FontWeight','bold');
text(0.07,0.41,'X = k\tau/(1+k\tau)', ...
    'Interpreter','tex', ...
    'Color',cream, ...
    'FontSize',20);
text(0.07,0.28,'PFR', ...
    'Color',pink, ...
    'FontSize',16, ...
    'FontWeight','bold');
text(0.07,0.21,'X = 1 - e^{-k\tau}', ...
    'Interpreter','tex', ...
    'Color',cream, ...
    'FontSize',20);
text(0.07,0.08, ...
    'Ideal first-order, isothermal reactor models', ...
    'Color',cream, ...
    'FontSize',12);
exportgraphics(fig3, ...
    'Reactor_Design_Equations.png', ...
    'Resolution',300);
%% VERIFY MODEL CONSISTENCY
assert(max(abs(X_batch-X_pfr)) < 1e-10, ...
    'Batch and PFR comparison failed.');
% CSTR conversion should not exceed PFR conversion.
assert(all(X_cstr <= X_pfr+1e-10), ...
    'CSTR/PFR comparison failed.');
fprintf('\nModel consistency checks passed.\n');
fprintf('Three portfolio figures exported.\n');