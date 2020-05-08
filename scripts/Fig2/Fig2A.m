clear, clc, close all

%% add necessary folders to the path
addpath(genpath(fullfile(cd ,'..\..')))

%% 

ratios = [0.3, 0.7];       % [mM0.5] Gbg to Gkirbar ratio 
lineStyles = {'--','-'};

fig = figure('units','inches');
fig.Renderer = 'painters';

hold all

for i = 1:numel(ratios)
    ratio = ratios(i);
    lineStyle = lineStyles{i};
    [t, Vm] = SingleCellModel(ratio);
    
    indx1 = find(t>=60.3,1,'first');
    indx2 = find(t>=90,1,'first');
    
    tplot1 = t(1:indx2);
    Vm_post1 = Vm(1:indx2);
    tplot2 = t(indx1:end);
    Vm_post2 = Vm(indx1:end);
    
    %% plot hysteresis
    plot(tplot1, Vm_post1,'b','linestyle',lineStyle, 'linewidth',2), hold on
    plot(tplot2-t(indx1), flipud(Vm_post2),'r','linestyle',lineStyle,'linewidth',2)
    
end

xlabel('Time(s)'),
ylabel('V_m(mV)')

ax = gca; box off;
ax.LineWidth = 3;
ax.FontSize = 22;
ax.FontName = 'arial';
ax.YLim = [-100, -20];
ax.YTick = [-100:20:-20];
ax.XTick = [0:20:100];

%% Equations

function [tplot, Vm_post] = SingleCellModel(ratio)

%% Simulation time
Tmax = 150;  % [s] Total time of the simulation
tspan = [0:0.01:Tmax]*1e3;      % converted to [ms]

%% Stimulation protocol
current_stim = false;     % current stimulus

%% Potassium stimulation
K_i = 150;                       % [mM] intracellular potassium concentration
K_out = @(t) 3*heaviside(t+1e3) + 2*heaviside(t-30e3) + 3*heaviside(t-60e3)...
    - 3*heaviside(t - 90e3) -2* heaviside(t - 120e3);

%% Current stimulation
I_stim_onset = 5;           % [s] current stimulation onset
I_stim_end = 10;            % [s] current stimulation end
current = -1;        % [pA]   Injected current

%% model parameters 

%% constant parameters
R  = 8314.0;	% [mJmol-1K-1]	gas constant
T = 293; % [K] absolute temperature
F  = 96485.0;	% [Cmol-1] Faraday's constant
RT_F = R*T/F;   % RT/F
z_K  = 1;       % K ion valence

%% Kir channel characteristic
delta_V_kir  = 25; % [mV]	voltage diff at half-max.
G_kirbar = 0.18;% [nS/mM^0.5] inward rectifier constant
n_kir = 0.5;	% inward rectifier constant
k_kir = 7;	% [mV]	inward rectifier slope factor

%% Background current
E_bg = -30;  % [mV] resting membrane potential
G_bg = ratio*G_kirbar;  % [nS] background conductance

%% Cell capacitance
Cm = 8;   % [pF] capillary membrane capacitance

%% Solution 
% initial condition
Vm0 = -30;      %[mV] initial condition for membrane potentials
init = Vm0;
opt = optimset;
opt.Display = 'off';
X0 = fsolve(@(X) Eqns(0,X), init, opt);

% Solve
tic
[t,X] = ode15s(@Eqns,tspan,X0);
toc
tplot = t/1000;     % convert to seconds

Vm_post	= X;

%% main differential equations
    function dXdt = Eqns(t,X)
        
        Vm = X;
        
        %% Stimulation protocol
        
        K_o = K_out(t);
        
        % Current stimulation
        I_stim = 0;
        if current_stim
            if t >= I_stim_onset*1e3 && t <= I_stim_end*1e3
                I_stim = current;
            else
                I_stim = 0;
            end
        end
        
        %% Potassium reversal potential
        E_K = RT_F/z_K*log(K_o./K_i);     %[mV]
        
        %% Membrane currents
        I_bg = G_bg.*(Vm - E_bg);        %[pA] lumped background current
        I_kir = G_kirbar.*(K_o).^n_kir .*((Vm - E_K)./(1 + exp((Vm - E_K - delta_V_kir)./k_kir)));   %[pA] whole cell kir current
        
        %% Total transmembrane current
        I_tot = I_bg + I_kir;
        
        %% Differential equations
        dVmdt = -1./Cm .* (I_tot - I_stim);
        dXdt = dVmdt;
    end
end

