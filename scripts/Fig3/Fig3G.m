clear,clc, close all

%% add necessary folders to the path
addpath(genpath(fullfile(cd ,'..\..')))

%%
stim_cell_number = [2, 4];  %number of stimulated cells
ratio = 6;      %Gkir/Gbg
for i = 1:numel(stim_cell_number)
    number = stim_cell_number(i);
    CapillarySegment(ratio, number)
end

%% Equations 
function CapillarySegment(ratio, number)
%% Simulation time
Tmax = 8;                % [s] Total time of the simulation
tspan = [0, Tmax*1e3];      % converted to [ms]
tot_cell_number = 20;

%% %%%%%%%%%%%%%%%% Stimulation protocol %%%%%%%%%%%%%%%%%%%%%%%%%
Options.K_stim = true;                % extracellular potassium stimulation

%% Potassium stimulation
K_stim_cell_indx = 1:number;
K_o_Stim = 10;              % [mM] extracellular potassium concentration for stimulated cells
K_o_rest = 3;                 % [mM] extracellular potassium concentration for unstimulated cells
K_stim_onset = 1;         % [s] K Stimulation onset
K_stim_end = 80;     % [s] K Stimulation end
K_i = 150;               % [mM] intracellular potassium concentration

%% %%%%%%%%%%%%%%%%%% model parameters %%%%%%%%%%%%%%%%%%%%%%%%

%% constant parameters
R  = 8314.0;	% [mJmol-1K-1]	gas constant
T = 293.0;      % [K] absolute temperature
F  = 96487.0;	% [Cmol-1] Faraday's constant
RT_F = R*T/F;   % RT/F
z_K  = 1;       % K ion valence

%% Background curren0t
E_bg = -30;        % [mV] resting membrane potential
G_bg = 0.06;%*ones(1,tot_cell_number); % [nS]  lumped background conductance

%% Kir channel characteristic
delta_V_kir  = 25;% [mV]	voltage diff at half-max. eI_K_i
G_kirbar = ratio*G_bg/sqrt(K_o_rest);%0.18;   % [nS/mM^0.5] inward rectifier constant
n_kir = 0.5;	% inward rectifier constant (0.5)
k_kir = 7;	% [mV]	inward rectifier slope factor

%% Gap junctional resistance
R_gj = 1e-2;  % [Gohm]  gap junctional resistance

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

%% store all parameters
ParNames = who();
for i = 1:numel(ParNames)
    eval(['Param.' ParNames{i} '=' ParNames{i} ';' ]);
end

%% pattern of the jacobian matrix
S = Adj;
S = S + eye(size(S));
S = sparse(S);

%% Solve
tic
Vm0 = -30*ones(1,tot_cell_number);      %[mV] initial condition for membrane potentials

tic
opts=odeset('JPattern',S);
[t,Vm] = ode15s(@(t,X)Eqns(t,X,Param,Options),tspan,Vm0,opts);
toc

%% ------------ plot Vm vs t
Fig = figure(1);
hold all
Fig.Renderer= 'painters';
Fig.Units = 'inches';
Fig.Position = [1 1 4.5 4.5];


tplot = t/1000;     % convert to seconds
unstim_cells = 1:tot_cell_number;
unstim_cells(Param.K_stim_cell_indx) = [];
plot(tplot, Vm(:,Param.K_stim_cell_indx),'g','linewidth',1), grid on, xlabel('t(s)'),ylabel('V_m (mV)'), hold on
plot(tplot, Vm(:,unstim_cells(1:end)),'k','linewidth',1), grid on, xlabel('t(s)'),ylabel('V_m (mV)'), hold on
plot(tplot, Vm(:,end),'r','linewidth',1)

ax = gca;
ax.FontSize = 16;
ax.LineWidth = 2;
ax.FontName = 'arial';
grid off, box off, ax.Color = 'none'; ax.YLim = [-100,-20];

end

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

