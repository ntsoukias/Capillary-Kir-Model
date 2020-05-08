clear, clc, close all

%% add necessary folders to the path
addpath(genpath(fullfile(cd ,'..\..')))

%% %%%%%%%%%%%%%%%% Stimulation time %%%%%%%%%%%%%%%%%%%%%%%%%
Tmax = 100;% [s] Total time of the simulation
dt = 0.0001; % [s] time increment
tspan = [0:dt*1e3:Tmax*1e3];      % [ms] time span

%% %%%%%%%%%%%%%%%% Stimulation protocol %%%%%%%%%%%%%%%%%%%%%%%%%
K_stim = true;           % bath potassium stimulation
current_stim = false;    % current stimulus
Options.PIP2 = 'on';    % dynamic regulation of Kir and TRPV4

%% geometry

axial_ArtEC = 10; % number of SMCs in the axial direction
circ_ArtEC = 3; % number of SMCs in the circumferential direction

number = axial_ArtEC*circ_ArtEC;  % total number of cells

tot_ArtECs = number;
ArtECs = 1:tot_ArtECs;  ArtECs = reshape(ArtECs,circ_ArtEC,axial_ArtEC);    %arrangement of SMCs (column = circular)
AdjArtECs = zeros(tot_ArtECs);
a = ones(1,length(diag(AdjArtECs,1)));
a(circ_ArtEC:circ_ArtEC:numel(a)) = 0;
AdjArtECs = diag(a,1);
a = ones(1,length(diag(AdjArtECs,circ_ArtEC)));
AdjArtECs = AdjArtECs + diag(a,circ_ArtEC);
a = zeros(1,length(diag(AdjArtECs,circ_ArtEC - 1)));
a(1:circ_ArtEC:end) = 1;
AdjArtECs = AdjArtECs + diag(a,circ_ArtEC - 1);
AdjArtECs = AdjArtECs+triu(AdjArtECs,1).';
AdjArtECs = AdjArtECs - diag(diag(AdjArtECs));

%% capillary branch
s = load('Geometry_Fig4_only_cap.mat');
AdjCapECs = s.A;

%% combined
AdjCombined = zeros(length(AdjArtECs) + length(AdjCapECs));
AdjCombined(1:length(AdjCapECs),1:length(AdjCapECs)) = AdjCapECs;
AdjCombined(length(AdjCapECs)+1:end,length(AdjCapECs)+1:end) = AdjArtECs;
ArtECsInds = length(AdjCapECs) + ArtECs;

AdjCombined(1,ArtECsInds(1,5)) = 1;
AdjCombined(ArtECsInds(1,5),1) = 1;
tot_cell_number = length(AdjCombined);
Index_CaP = AdjCombined;

g = graph(AdjCombined);

%% find neighbors
nnod = g.numnodes;
Param.edgelist = g.Edges{:,1};
neighbors = zeros(nnod, max(g.degree));

for i = 1:nnod
    neighborNodes = g.neighbors(i);
    neighborNodes = neighborNodes';
    neighbors(i,1:numel(neighborNodes)) = neighborNodes;
end

Param.neighbors = neighbors;
Param.nnod = nnod;

%% index of parenchymal arterioles
PAECs = ArtECsInds;
CapECs = 1:tot_cell_number;
CapECs(PAECs(:)) = [];

%% Bath Potassium stimulation

K_stim_cell_indx = [34:36 42:45];

K_bath_stim = 10;        % [mM] bath potassium concentration for stimulated cells
K_bath_rest = 3;        % [mM] normal bath potassium concentration
K_stim_onset = 20;	% [s] K Stimulation onset
K_stim_end = 50;	  % [s] K Stimulation end
Ki = 150;

PIP2_cells = 1:tot_cell_number;%randi(tot_cell_number, 1, round(tot_cell_number/4));
tau_kir = 3e3*60;    % [ms] inhibition by PIP2 depletion time constant
tau_bg = 6e3*60;     % [ms] activation by PIP2 depletion time constant
lag = 5;  %[s] lag in activation of PIP2

%% Current stimulation

I_stim_cell_indx = [34:36 42:45];
I_stim_onset = 20;           % [s] current stimulation onset
I_stim_end = 50;            % [s] current stimulation end
I_ext = -4;      % [pA] Injected current

%% %%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%

%% constant parameters
R = 8314;  %[mJ/mol/K]  universal gas constant
T = 293; %[K] temperature
F = 96485;    %[C/mol]    Faraday's constant
RT_F = R*T/F;   %[mV] RT/F
z_Na = 1;    %valence of sodium
z_Ca = 2;   %valence of chloride
z_K = 1; %valence of potassium

%% KIR CHANNEL CHARACTERISTICS

delta_V_kir  = 25;      % [mV]	voltage diff at half-max. eI_K_i
G_kirbarMean = 0.8/sqrt(3)*ones(1,tot_cell_number);   % [nS/mM^0.5] inward rectifier constant
G_kirbarMean(ArtECsInds(:)) = 0.65/sqrt(3);

G_kirbarStd = 0.1;
G_kirbar = G_kirbarMean + G_kirbarStd * randn(1,tot_cell_number);
G_kirbar(ArtECsInds(:)) = G_kirbarMean(ArtECsInds(:));
G_kirbar(G_kirbar<=0) = 0;

n_kir = 0.5;	% inward rectifier constant
k_kir = 7;	% [mV]	inward rectifier slope factor

%% BACKGROUND CURRENT
E_bg = -30; % [mV] resting membrane potential
G_bg = 0.06*ones(1,tot_cell_number); % [nS]  lumped background conductance
G_bg(ArtECsInds(:)) = 0.2;

%% TRPV4 currents

TRPV4_type = '3 states 1 channel';

Cm = 8;   % [pF] capillary membrane capacitance
A_m = Cm./1e6;	%[cm2]	membrane area
Cao = 2.0;  %[mM] [Ca2+]o extracellular calcium
Cai = 70e-6;    %[mM] [Ca2+]i intracellular calcium
Nai = 10;
Nao = 150;
p = 1e-7;  %[cm./sec]
PTRPV4Ca= p; %[cm./sec]
PTRPV4Na=p./6.9; %%% PCa:PK:PNa= 6.9:1.35:1
PTRPV4K=1.35.*p./6.9;

%% generate new traces
disp('Generating stochastic TRPV4 traces....')
meanshutlifetime= 240e3;%msec
meanopenlifetime = 90;  %msec
meanburstlength = 20e3; %msec
meanblockedlifetime = 33; %msec
TRPV4_pars = [meanshutlifetime, meanopenlifetime, ...
    meanburstlength, meanblockedlifetime];
TRPV_P0 = zeros(tot_cell_number, length(tspan));
for cells = CapECs
    TRPV_P0(cells,:) = TRPV4_P0_threestate1channel(tspan,TRPV4_pars);
end

%% load existing traces
% s = load('Fig4BCDE_TRPV_P0.mat');
% TRPV_P0 = s.TRPV_P0;

%% Gap junctional resistance
R_gj = 1e-2;  % [Gohm]  gap junctional resistance

%% membrane capacitance
Cm = 8*ones(1,tot_cell_number); %pF
Cm(ArtECsInds(:)) = 16;

%% %%%%%%%%%%%%%%%% Initial conditions %%%%%%%%%%%%%%%%%%%%%%%%%
Vm_in = -30*ones(1,tot_cell_number);  %[mV]
X0 = Vm_in;

%% store pars
ParNames = who();
for i = 1:numel(ParNames)
    eval(['Param.' ParNames{i} '=' ParNames{i} ';' ]);
end

%% %%%%%%%%%%%%%%%%%% Solution %%%%%%%%%%%%%%%%%%%%%%%%%
disp('solving...')

% pattern of the jacobian matrix
S = AdjCombined;
S = S + eye(size(S));
S = sparse(S);

tic
opts=odeset('OutputFcn',@odeprog,'JPattern',S);
[t,X] = ode15s(@(t,X)EqnsNew(t,X,Param,Options),[0 Tmax*1e3],X0, opts);
toc

%% visualization
ii = 1;
Vm_post = X(:,ii:tot_cell_number);     % [mV] membrane voltage

%% plot Vm vs t

fig = figure;
FigName = 'VmvsTime';
fig.Name = FigName;
fig.Renderer = 'painters';
fig.Units = 'inches';
fig.Position = [1 1 5 4.5];
hold all


tplot = t/1000;     % convert to seconds
cells = 1:tot_cell_number;
cells(PAECs) = [];
cECinds = cells;
plot(tplot, Vm_post(:,cells),'color',[0.7 0.7 0.7],'linewidth',2), grid on, 
xlabel('t(s)'),ylabel('V_m(mV)'), hold on
plot(tplot, Vm_post(:,PAECs),'k','linewidth',3), grid on, 

ylim([-90,-20])
xlim([0, 100])
ax = gca;
ax.YTick = [-80:20:-20];
ax.XTick = [0:20:100];
ax.FontSize = 22;
ax.FontName = 'arial';
ax.LineWidth = 3;
box off
grid off

fig = gcf;
fig.Color = 'w';

%% %%%%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%
function dXdt = EqnsNew(t,X,Param,Options)

%% parameters
tspan = Param.tspan;
tot_cell_number = Param.tot_cell_number;
K_stim =  Param.K_stim;
current_stim= Param.current_stim;
K_stim_cell_indx = Param.K_stim_cell_indx;
K_bath_stim= Param.K_bath_stim;
K_bath_rest= Param.K_bath_rest;
K_stim_onset= Param.K_stim_onset;
K_stim_end= Param.K_stim_end;
I_stim_cell_indx = Param.I_stim_cell_indx;
I_stim_onset= Param.I_stim_onset;
I_stim_end= Param.I_stim_end;
I_ext= Param.I_ext;
F= Param.F;
R = Param.R;
T = Param.T;
RT_F = Param.RT_F;
Ki = Param.Ki;
z_Na= Param.z_Na;
z_K= Param.z_K;
delta_V_kir= Param.delta_V_kir;
G_kirbar= Param.G_kirbar;
n_kir= Param.n_kir;
k_kir= Param.k_kir;
E_bg= Param.E_bg;
G_bg= Param.G_bg;
Cm= Param.Cm;
A_m= Param.A_m;
Cao= Param.Cao;
Cai= Param.Cai;
Nai = Param.Nai;
Nao = Param.Nao;
PTRPV4Ca= Param.PTRPV4Ca;
PTRPV4Na= Param.PTRPV4Na;
PTRPV4K= Param.PTRPV4K;
TRPV_P0= Param.TRPV_P0;
R_gj= Param.R_gj;
PIP2_cells = Param.PIP2_cells;
z_Ca = Param.z_Ca;
neighbors = Param.neighbors;
nnod = Param.nnod;
tau_kir = Param.tau_kir;
tau_bg = Param.tau_bg;
lag = Param.lag;

%% state variabels
X = X(:);
ii = 1;
Vm = X(ii:tot_cell_number)';     % [mV] membrane voltage
dXdt = zeros(numel(X), 1);  % stores the vector of first order differential EqnsNew

%% stimulation protocol
% ------------------ bath potassium stimulus
K_bath = K_bath_rest*ones(1,tot_cell_number); % current stimulus (uA/cm^2)

fact_kir_basal = ones(1,tot_cell_number);
fact_bg_basal = ones(1,tot_cell_number);

fact_kir = fact_kir_basal;
fact_bg = fact_bg_basal;

if K_stim
    if t>= K_stim_onset*1e3 && t<=K_stim_end*1e3
        K_bath(K_stim_cell_indx) = K_bath_stim;
    end
    if t > (K_stim_onset + lag)*1e3 && strcmp(Options.PIP2,'on')
        fact_kir(PIP2_cells) = fact_kir_basal(PIP2_cells) - (1 - exp(-(t - K_stim_onset*1e3)/tau_kir));
        fact_bg(PIP2_cells)= fact_bg_basal(PIP2_cells) + 2*(1 - exp(-(t - K_stim_onset*1e3)/tau_bg));
    end
end

Ko = K_bath;

% ------------------- current stimulus
I_stim = zeros(1,tot_cell_number); % current stimulus (uA/cm^2)
if current_stim
    if t>= I_stim_onset*1e3 && t<=I_stim_end*1e3
        I_stim(I_stim_cell_indx) = I_ext;
    end
end

%% Reversal potentials
EK = RT_F/z_K*log(Ko./Ki);    %[mV] reversal potential for sodium ions

%% Membrane currents

I_bg = fact_bg.*G_bg.*(Vm - E_bg);        %[pA] nonspecific background current
I_Kir = fact_kir.*G_kirbar.*(Ko).^n_kir .*((Vm - EK)./(1 + exp((Vm - EK - delta_V_kir)./k_kir)));   %[pA] whole cell kir current

%% TRPV4 current
if Vm ~= 0
    I_TRPV4Na = 1e6.*A_m.*PTRPV4Na.*z_Na.^2.*F.^2./(R.*T).*Vm.*(Nai-Nao.*exp(-Vm.*(z_Na.*F)./(R.*T)))./(1-exp(-Vm.*(z_Na.*F)./(R.*T)));
    I_TRPV4Ca = 1e6.*A_m.*PTRPV4Ca.*z_Ca.^2.*F.^2./(R.*T).*Vm.*(Cai-Cao.*exp(-Vm.*(z_Ca.*F)./(R.*T)))./(1-exp(-Vm.*(z_Ca.*F)./(R.*T)));
    I_TRPV4K =  1e6.*A_m.*PTRPV4K .*z_K.^2 .*F.^2./(R.*T).*Vm.*(Ki-Ko.*exp(-Vm.*(z_K.*F)./(R.*T)))./(1-exp(-Vm.*(z_K.*F)./(R.*T)));
else
    I_TRPV4Na = 1e6.*A_m.*PTRPV4Na.*z_Na.*F.*(Nai-Nao);
    I_TRPV4Ca = 1e6.*A_m.*PTRPV4Ca.*z_Ca.*F.*(Cai-Cao);
    I_TRPV4K =  1e6.*A_m.*PTRPV4K .*z_K.*F.*(Ki-Ko);
end

ind = find(tspan>=t,1);
if ~isempty(ind)
    P0_TPRV4 = TRPV_P0(:,ind)';    
else
    P0_TPRV4 = rand(size(I_TRPV4K));
end

I_TRPV4Na = P0_TPRV4.*I_TRPV4Na;
I_TRPV4K =  P0_TPRV4.*I_TRPV4K;
I_TRPV4Ca =  P0_TPRV4.*I_TRPV4Ca;

I_TRPV4 = I_TRPV4K + I_TRPV4Na + I_TRPV4Ca;
I_TRPV4 = fact_bg.*I_TRPV4;

%% Gap junctional current [pA]
I_gj = zeros(1,nnod);
for i = 1:nnod    
    for j = 1:nnz(neighbors(i,:))
            I_gj(i) = I_gj(i) + 1/R_gj*(Vm(i) - Vm(neighbors(i,j)));      
    end    
end

%% differential EqnsNew for moles of ionic spicies; oxygen,  and membrane voltage
dVmdt = -1./Cm.*( I_bg + I_Kir + I_TRPV4 + I_gj - I_stim);

%% all differential equations
ii = 1;
dXdt(ii:tot_cell_number) = dVmdt;      % [mV] membrane voltage

end


function [P0]= TRPV4_P0_threestate1channel(tspan,TRPV4_pars)

i = 1;
meanshutlifetime = TRPV4_pars(i); i = i + 1;
meanopenlifetime = TRPV4_pars(i); i = i + 1;
meanburstlength = TRPV4_pars(i); i = i + 1;
meanblockedlifetime =  TRPV4_pars(i);

k4=1/meanshutlifetime;
k2=1/meanblockedlifetime;
k3andk1=[1 1;meanburstlength -1/k2]\[1/meanopenlifetime;1];
k3=k3andk1(1);
k1=k3andk1(2);

Q = [-(k3+k1) k1 k3; k2 -k2 0; k4 0 -k4];   %transition rate matrix
P = expm(Q);    % transition probability matrix

mc = dtmc(P);
trace = simulate(mc,length(tspan));

P0 = trace;
P0(ismember(trace,[2,3])) = 0;
P0(1,:) = [];
P0 = P0';
end

    
    
