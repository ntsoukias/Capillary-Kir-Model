clear,clc, close all

%% add necessary folders to the path
addpath(genpath(fullfile(cd ,'..\..')))

%% Simulation time
Tmax = 100;                % [s] Total time of the simulation
tspan = [0, Tmax*1e3];      % converted to [ms]
tot_cell_number = 20;
probeCell = 20;     % cell index for Vm reading

%% %%%%%%%%%%%%%%%% Stimulation protocol %%%%%%%%%%%%%%%%%%%%%%%%%
Options.K_stim = true;                % extracellular potassium stimulation

%% Potassium stimulation
K_o_Stim = 10;              % [mM] extracellular potassium concentration for stimulated cells
K_o_rest = 3;                 % [mM] extracellular potassium concentration for unstimulated cells
K_stim_onset = 20;         % [s] K Stimulation onset
K_stim_end = 80;     % [s] K Stimulation end
K_i = 150;                  % [mM] intracellular potassium concentration

%% %%%%%%%%%%%%%%%%%% model parameters %%%%%%%%%%%%%%%%%%%%%%%%

%% constant parameters
R  = 8314.0;	% [mJmol-1K-1]	gas constant
T = 293.0;      % [K] absolute temperature
F  = 96487.0;	% [Cmol-1] Faraday's constant
RT_F = R*T/F;   % RT/F
z_K  = 1;       % K ion valence

%% Kir channel characteristic
delta_V_kir  = 25;   % [mV]	voltage diff at half-max. eI_K_i 
n_kir = 0.5;	% inward rectifier constant (0.5)
k_kir = 7;	% [mV]	inward rectifier slope factor

%% Background current
E_bg = -30;        % [mV] resting membrane potential
G_bg = 0.06*ones(1,tot_cell_number); % [nS]  lumped background conductance

%% Cell capacitance
Cm = 8;   % [pF] capillary membrane capacitance

%% find adjacency matrix
parent_CaP = [0:tot_cell_number-1];
n = length(parent_CaP);
m = find(parent_CaP);
Adj = sparse (parent_CaP(m), m, 1, n, n);     % creates an upper triangular matrix with 1 if cells i and j are connected and zero if not (sparse matrix)
Adj = full(Adj);
Adj = Adj+triu(Adj,1).';              % converts the upper triangular matrix to symmetric

%% find neighbors
G = graph(Adj);
nnod = numnodes(G);
edgelist = G.Edges{:,1};
neighbors = zeros(nnod, max(G.degree));

for i = 1:nnod
    neighborNodes = G.neighbors(i);
    neighborNodes = neighborNodes';
    neighbors(i,1:numel(neighborNodes)) = neighborNodes;
end

%% pattern of the jacobian matrix
S = Adj;
S = S + eye(size(S));
S = sparse(S);

%%
Vm0 = -30*ones(1,tot_cell_number);      %[mV] initial condition for membrane potentials
R_gj_vals = logspace(-3,1.4,200);
G_kirbar_vals = [0:0.1:8.1]*G_bg(1);   % [nS/mM^0.5] inward rectifier constant

stim_number = 2;
all_deltaVms = zeros(numel(G_kirbar_vals),tot_cell_number,numel(R_gj_vals));

%% store all parameters
ParNames = who();
for i = 1:numel(ParNames)
    eval(['Param.' ParNames{i} '=' ParNames{i} ';' ]);
end

%% solve

opts=odeset('JPattern',S);

f = waitbar(0, 'Solving ...');

for i = 1:numel(G_kirbar_vals)
    for j = 1:numel(R_gj_vals)
        number = stim_number;
        Param.K_stim_cell_indx = 1:number;
        Param.R_gj = R_gj_vals(j);      % [Gohm]  gap junctional resistance
        Param.G_kirbar = G_kirbar_vals(i);   % [nS/mM^0.5] inward rectifier 

        
        [t,Vm] = ode15s(@(t,X)Eqns(t,X,Param,Options),tspan,Vm0,opts);
        tplot = t/1000;
        t_indx1 = find(tplot<Param.K_stim_onset,1,'last');
        t_indx2 = find(tplot<Param.K_stim_end,1,'last');
        all_deltaVms(i,:,j) = Vm(t_indx2,:) - Vm(t_indx1,:);    
    end
    waitbar(i/numel(G_kirbar_vals), f);
end
close(f)

%% plots
x = R_gj_vals;
y = G_kirbar_vals*sqrt(K_o_rest)/G_bg(1);
[X,Y] = meshgrid(x,y);

fig = figure;
fig.Name = 'Vmvst';
fig.Renderer = 'painters';
fig.Units = 'inches';
fig.Position = [1 1 4.5 4.5];

contourf(X,Y,squeeze(all_deltaVms(:,probeCell,:)),40);

axis tight
ax = gca;
ax.XAxis.Scale = 'log';
ax.LineWidth = 3;
ax.FontSize = 16;
ax.FontName = 'arial';
ax.XTick = [1e-3,1e-1,1e1];
ax.YTick = 4:4:20;
xlabel('R_{gj}(G\Omega)')
ylabel('G_{K_{IR}}/G_{bg}')

box off
colors = flipud(jet(50));
colormap(colors)
c = colorbar;
c.FontSize = 12;
ylabel(c,'\DeltaV_m(mV)');

%% main equations
function dVmdt = Eqns(t,X,Param,Options)

X = X(:);
Vm = X';

%% options
K_stim = Options.K_stim;               % extracellular potassium stimulation

%% parameters
tot_cell_number = Param.tot_cell_number;
K_stim_cell_indx = Param.K_stim_cell_indx; 
K_o_Stim = Param.K_o_Stim;
K_o_rest = Param.K_o_rest; 
K_stim_onset = Param.K_stim_onset; 
K_stim_end = Param.K_stim_end;
K_i = Param.K_i; 
RT_F = Param.RT_F;
z_K = Param.z_K;
E_bg = Param.E_bg;
G_bg = Param.G_bg;
delta_V_kir = Param.delta_V_kir; 
G_kirbar = Param.G_kirbar;
n_kir = Param.n_kir;
k_kir = Param.k_kir;
R_gj = Param.R_gj; 
Cm = Param.Cm;
neighbors = Param.neighbors;

%% --------------------------------------------------- Stimulation protocol
% Potassium stimulus
K_o = K_o_rest*ones(1,tot_cell_number);

if K_stim
    if t > K_stim_onset*1e3 && t <= K_stim_end*1e3
        K_o(K_stim_cell_indx) = K_o_Stim;
    end
end

% Current stimulation
I_stim = zeros(1,tot_cell_number);

%% Potassium equilibrium potential
E_K= RT_F/z_K*log(K_o./K_i);     %[mV]

%% Membrane currents [pA]
I_bg = G_bg.*(Vm - E_bg);        %[pA] lumped background current
I_kir = G_kirbar*(K_o).^n_kir .*((Vm - E_K)./(1 + exp((Vm - E_K - delta_V_kir)/k_kir)));   %[pA] whole cell kir current

%% Gap junctional current [pA]

nnod = tot_cell_number;

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



