clear, clc, close all

%% add necessary folders to the path
addpath(genpath(fullfile(cd ,'..\..')))

%%
G_vals = [0.2, 0.4];    %[nS] Gkirbar values
stim_conds = [1, 1];    % extracelluar potassium stimulation 1-on, 0-off

for i = 1:numel(G_vals)
    stim = stim_conds(i);                % extracellular potassium stimulation
    G = G_vals(i);
    CapillarySegment(G,stim)
end

%% functions
function CapillarySegment(G,stim)

%% Simulation time
Tmax = 8;                % [s] Total time of the simulation
tspan = [0 Tmax*1e3];      % converted to [ms]
tot_cell_number = 20;

%% Stimulation protocol
K_stim = stim;                % extracellular potassium stimulation
wash_out = false;          % true if you want to remove the potassium stimulus to a new value (K_o_washout) after stim_time
current_stim = false;     % current stimulus

%% Potassium stimulation
K_stim_cell_indx = 1:2;
K_o_Stim = 10;              % [mM] extracellular potassium concentration for stimulated cells
K_o_rest = 3;                 % [mM] extracellular potassium concentration for unstimulated cells
K_stim_onset = 2;         % [s] K Stimulation onset
K_stim_end = 6;     % [s] K Stimulation end
K_i = 150;               % [mM] intracellular potassium concentration

%% potassium washout
K_o_washout = K_o_rest;       % [mM] potassium level after the stimulus is removed (after stim_time)
K_washout_onset = 25;
K_washout_end = Tmax;

%% Current stimulation
I_stim_cell_indx = [1:10];
I_stim_onset = 5;           % [s] current stimulation onset
I_stim_end = 10;            % [s] current stimulation end
current = -1;        % [pA]   Injected current

%% model parameters

%% constant parameters
R  = 8314.0;	% [mJmol-1K-1]	gas constant
T = 293.0;      % [K] absolute temperature
F  = 96487.0;	% [Cmol-1] Faraday's constant
RT_F = R*T/F;   % RT/F
z_K  = 1;       % K ion valence

%% Kir channel characteristic
delta_V_kir  = 25;% [mV]	voltage diff at half-max. eI_K_i
G_kirbar = G/sqrt(K_o_rest);  % [nS/mM^0.5] inward rectifier constant
n_kir = 0.5;	% inward rectifier constant (0.5)
k_kir = 7;	% [mV]	inward rectifier slope factor

%% Background current
E_bg = -30;        % [mV] resting membrane potential
G_bg = 0.06*ones(1,tot_cell_number); % [nS]  lumped background conductance

%% Gap junctional resistance
R_gj = 1e-2; % [Gohm]  gap junctional resistance

%% Cell capacitance
Cm = 8;   % [pF] capillary membrane capacitance

%% find adjacency matrix
parent_CaP = [0:tot_cell_number-1];
n = length(parent_CaP);
m = find(parent_CaP);
Index_CaP = sparse (parent_CaP(m), m, 1, n, n);     % creates an upper triangular matrix with 1 if cells i and j are connected and zero if not (sparse matrix)
Index_CaP = full(Index_CaP);
Index_CaP = Index_CaP+triu(Index_CaP,1).';              % converts the upper triangular matrix to symmetric

%% find neighbors
G = graph(Index_CaP);
nnod = numnodes(G);
neighbors = zeros(nnod, max(G.degree));

for i = 1:tot_cell_number
    neighborNodes = G.neighbors(i);
    neighborNodes = neighborNodes';
    neighbors(i,1:numel(neighborNodes)) = neighborNodes;
end

%% pattern of the jacobian matrix
S = Index_CaP;
S = S + eye(size(S));
S = sparse(S);

%% Solution
% initial condition
Vm0 = -20*ones(1,tot_cell_number);      %[mV] initial condition for membrane potentials
opt = optimset;
opt.Display = 'off';
Vm0 = fsolve(@(X)Eqns(0,X), Vm0, opt);

% Solve using odesolver
tic
opts=odeset('JPattern',S);
[t,Vm] = ode15s(@Eqns,tspan,Vm0,opts);
toc

%% plot Vm vs t
fig = figure(1);
fig.Name = 'Vmvst';
fig.Renderer = 'painters';
fig.Units = 'inches';
fig.Position = [1 1 4.5 4.5];

tplot = t/1000;     % convert to seconds
if tot_cell_number > 1
    unstim_cells = 1:tot_cell_number;
    ylabel('V_m(mV)'),
    hold on
    plot(tplot, Vm(:,unstim_cells(1:end-1)),'k','linewidth',2), hold on
    plot(tplot, Vm(:,unstim_cells(end)),'r','linewidth',1), hold on
    plot(tplot, Vm(:,K_stim_cell_indx),'g','linewidth',1),
    xlabel('Time(s)'),
else
    plot(tplot, Vm,'k','linewidth',2), grid on, xlabel('t(s)'),ylabel('V_m(mV)'), hold on
end
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 22;
ax.YLim = [-100,-20];
ax.YTick = [-100:20:-20];
ax.XLim = [0,8];
ax.XTick = [0:2:8];
ax.FontName = 'arial';
box off; grid off

%% plot steady state
fig = figure(2);
fig.Name = 'ssVm';
fig.Renderer = 'painters';
fig.Units = 'inches';
fig.Position = [1 1 4.5 4.5];%4.5];
timePoint = 5;
indx = find(tplot>=timePoint,1,'first');

if tot_cell_number > 1
    unstim_cells = 1:tot_cell_number;
    unstim_cells(K_stim_cell_indx) = [];
    plot(1:tot_cell_number, Vm(indx,:),'k-o','linewidth',2),hold on
    plot(K_stim_cell_indx, Vm(indx,K_stim_cell_indx),'g-o','linewidth',2)
    plot(unstim_cells, Vm(indx,unstim_cells),'k-o','linewidth',2)
    plot(tot_cell_number, Vm(indx,end),'r-o','linewidth',2), hold on
    plot(1:tot_cell_number,-30*ones(1,tot_cell_number),'k-.','linewidth',2)
else
    plot(tplot, Vm,'k','linewidth',2), grid on,
    hold on
end
xlabel('Cell#'),
ylabel('V_m(mV)'),
ax = gca;
ax.XTick = [0:10:20];
ax.YLim = [-100,-20];
ax.YTick = [-100:20:-20];
ax.LineWidth = 3;
ax.FontSize = 22;
ax.FontName = 'arial';
box off; grid off


    function dVmdt = Eqns(t,X)
        
        Vm(1:tot_cell_number) = X;
        
        %% --------------------------------------------------- Stimulation protocol
        % Potassium stimulus
        K_o = K_o_rest*ones(1,tot_cell_number);
        
        if K_stim
            if t > K_stim_onset*1e3 && t <= K_stim_end*1e3
                K_o(K_stim_cell_indx) = K_o_Stim;
            else
                K_o(K_stim_cell_indx) = K_o_rest;
            end
        end
        
        % Potassium washout
        if wash_out
            if t > K_washout_onset*1e3   && t <= K_washout_end*1e3
                K_o(K_stim_cell_indx) = K_o_washout;
            end
        end
        
        % Current stimulation
        I_stim = zeros(1,tot_cell_number);
        
        if current_stim
            if t >= I_stim_onset*1e3 && t <= I_stim_end*1e3
                I_stim(I_stim_cell_indx) = current;
            else
                I_stim = 0;
            end
        end
        
        %% Potassium equilibrium potential
        E_K= RT_F/z_K*log(K_o./K_i);     %[mV]
        
        %% Membrane currents [pA]
        I_bg = G_bg.*(Vm - E_bg);        %[pA] lumped background current
        I_kir = G_kirbar*(K_o).^n_kir .*((Vm - E_K)./(1 + exp((Vm - E_K - delta_V_kir)/k_kir)));   %[pA] whole cell kir current
        
        %% Gap junctional current [pA]
        
        I_gj = zeros(1,nnod);
        for i = 1:nnod
            for j = 1:nnz(neighbors(i,:))
                I_gj(i) = I_gj(i) + 1/R_gj*(Vm(i) - Vm(neighbors(i,j)));
            end
        end
        
        %% Total membrane current
        I_tot = I_bg + I_kir;
        
        %% Differential equation
        dVmdt =  -1/Cm * (I_tot + I_gj - I_stim);
        dVmdt = dVmdt';
        
    end
end