clear, clc, close all

%% add necessary folders to the path
addpath(genpath(fullfile(cd ,'..\..')))

%% parameters
% ----------  constant parameters
R  = 8314.0;	% [mJmol-1K-1]	gas constant
T = 293;   % [K] absolute temperature
F  = 96485.0;	% [Cmol-1] Faraday's constant
RT_F = R*T/F;   % RT/F
z_K  = 1;    % K ion valence
Cm = 8;     %[pF] membrane capacitance
Ki = 150;   %[mM] intracellular potassium

%----------  Kir channel characteristic
delta_V_kir  = 25; % [mV]	voltage diff at half-max. eI_K_i
G_kirbar = 0.18;   % [nS/mM^0.5] inward rectifier constant
n_kir = 0.5;	% inward rectifier constant
k_kir = 7;	% [mV]	inward rectifier slope factor

%----------- Background current
E_bg = -30;  % [mV] resting membrane potential

%% plot kir and bg at 3, 5, and 8

Ko_vals = [3, 5, 8]; %[mM]
ratio = 0.7;    % Gbg to Gkirbar ratio [mM0.5]
G_bg_desired = G_kirbar*ratio;
Vmplot = linspace(-120, 0, 300);

for i = 1:numel(Ko_vals)
    Ko = Ko_vals(i);
    E_K = RT_F/z_K*log(Ko/Ki);
    
    I_kirMean = G_kirbar*(Ko).^n_kir .*...
        ((Vmplot - E_K)./(1 + exp((Vmplot - E_K - delta_V_kir)./k_kir)));   %[pA] whole cell kir current
    I_bgMean = G_bg_desired*(Vmplot - E_bg);
    
    fig = figure;
    plot(Vmplot, [I_kirMean; -I_bgMean],'linewidth',3)
    axis([-110,0,-5,8])
    xlabel('V_m(mV)'), ylabel('I(pA)'),
    ax = gca;
    ax.FontSize = 24;
    ax.FontName = 'arial';
    ax.LineWidth = 3;
    grid off; box off
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    ax.XTick = [-90:30:0];
    ax.XLabel.Position = [-43.8596   -3.8881   -1.0000];
    fig.Color = 'w';
    ax.YTick = [-4:4:8];
    ax.YLim = [-4,8];

end

