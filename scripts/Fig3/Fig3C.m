clear, clc, close all

%% add necessary folders to the path
addpath(genpath(fullfile(cd ,'..\..')))

%% Simulation time
Tmax = 80;    % [s] Total time of the simulation
dt = 0.0001;  %[s] time increment
tspan = 0:dt*1e3:Tmax*1e3;      % converted to [ms]

%% define models for Kir and TRPV4 currents

tot_cell_number = 20;

parent_CaP = [0:tot_cell_number-1];
n = length(parent_CaP);
m = find(parent_CaP);
Index_CaP = sparse (parent_CaP(m), m, 1, n, n);     % creates an upper triangular matrix with 1 if cells i and j are connected and zero if not (sparse matrix)
Index_CaP = full(Index_CaP);
Index_CaP = Index_CaP+triu(Index_CaP,1).';              % converts the upper triangular matrix to symmetric

%% %%%%%%%%%%%%%%%% Stimulation protocol %%%%%%%%%%%%%%%%%%%%%%%%%
K_stim = true;                % extracellular potassium stimulation
wash_out = false;          % true if you want to remove the potassium stimulus to a new value (K_o_washout) after stim_time
current_stim = false;     % current stimulus

%% Potassium stimulation
K_stim_cell_indx = 1:3;
K_o_Stim = 10;              % [mM] extracellular potassium concentration for stimulated cells
K_o_rest = 3;                 % [mM] extracellular potassium concentration for unstimulated cells
K_stim_onset = 20;         % [s] K Stimulation onset
K_stim_end = 50;            % [s] K Stimulation end
K_i = 150;   % [mM] intracellular potassium concentration

%% potassium washout
K_o_washout = K_o_rest;       % [mM] potassium level after the stimulus is removed (after stim_time)
K_washout_onset = 25;
K_washout_end = Tmax;

%% Current stimulation
I_stim_onset = 5;           % [s] current stimulation onset
I_stim_end = 10;            % [s] current stimulation end
current = -1;        % [pA]   Injected current

%% %%%%%%%%%%%%%%%%%% model parameters %%%%%%%%%%%%%%%%%%%%%%%%

%% constant parameters
R  = 8314.0;	% [mJmol-1K-1]	gas constant
T = 293;   % [K] absolute temperature
F  = 96485.0;	% [Cmol-1] Faraday's constant
RT_F = R*T/F;   % RT/F
z_K  = 1;       % K ion valence

%% Kir channel characteristic
delta_V_kir  = 25*ones(1,tot_cell_number);	% [mV]	voltage diff at half-max. eI_K_i 
G_kir = 0.18*ones(1,tot_cell_number); % [nS/mM^0.5] inward rectifier constant

n_kir = 0.5;	% inward rectifier constant (0.5)
k_kir = 7;	% [mV]	inward rectifier slope factor

%% Background current
E_bg = -30*ones(1,tot_cell_number);        % [mV] resting membrane potential
G_bg = 0.06; % [nS]  lumped background conductance

%% gap junctional resistance
R_gj = 1e-2;        %[Gohm] gap junctional resistance

%% Cell capacitance
Cm = 8*ones(1,tot_cell_number);   % [pF] capillary membrane capacitance

%% -------------------------------------------- TRPV4 current characteristics (stochastic)

% parameters
z_Na = 1	;	%;		Na ion valence
z_Ca = 2	;	%Ca ion valence
A_m = Cm./1e6;	%[cm2]	membrane area
Ca_o = 2.0		;%[mM]	[Ca2+]o extracellular calcium
Na_o = 134.0	; %[mM]	[Na+]o  extracellular sodium
Ca_i = 0.2; 	%[mM]	[Ca2+]i intracellular calcium
Na_i =  10 ;	%[mM]	[Na+]i  intracellular sodium
PCa= 1e-7; %[cm/sec]

PTRPV4Ca= PCa;                       %[cm./sec]
PTRPV4Na=PCa./6.9;                 % PCa:PK:PNa= 6.9:1.35:1
PTRPV4K=1.35.*PCa./6.9;

%% uncomment this section and comment next section to reproduce Fig3C
 s = load('Fig3C_TRPV_P0.mat');
 TRPV_P0 = s.TRPV_P0;

%% uncomment to construct new TRPV4 traces for cells

% disp('Generating TRPV4 traces....')
% tic
% 
% % --------------------------------- %'3 states 1 channel'
% meanshutlifetime= 240e3; %msec
% meanopenlifetime = 90;   %msec
% meanburstlength = 20e3; %msec
% meanblockedlifetime = 33; %msec
% TRPV4_pars = [meanshutlifetime, meanopenlifetime, ...
%     meanburstlength, meanblockedlifetime];
% 
% TRPV_P0 = zeros(tot_cell_number, length(tspan));
% for cells = 1:tot_cell_number
%     TRPV_P0(cells,:) = TRPV4_P0_3state1channel(tspan,TRPV4_pars);
% end
% 
% toc

%% store all parameters
ParNames = who();
for i = 1:numel(ParNames)
    eval(['pars.' ParNames{i} '=' ParNames{i} ';' ]);
end

%% %%%%%%%%%%%%%%%%%% Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial condition
Vm0 = -30*ones(1,tot_cell_number);      %[mV] initial condition for membrane potentials
init = [Vm0];
opt = optimset;
opt.Display = 'off';
X0 = fsolve(@(X) Eqns(0,X,pars), init, opt);

% Solve
disp('solving......')
tic
opts=odeset('OutputFcn',@odeprog);
[t,X] = ode15s(@(t,X) Eqns(t,X,pars),[0 Tmax*1e3],X0,opts);
toc

Vm_post	= X;

%% PLOT Vm vs. Time
fig = figure;
fig.Name = 'VmvsTime';
fig.Renderer = 'painters';
fig.Units = 'inches';
fig.Position = [1 1 4.5 4.5];

cells = 1:tot_cell_number;
cells(K_stim_cell_indx) = [];
unstim_cells = cells;

tplot = t/1000;     % convert to seconds
h2 = plot(tplot, Vm_post(:,unstim_cells),'k','linewidth',2); hold on
h4 = plot(tplot, Vm_post(:,end),'r','linewidth',0.5); hold on
h3 = plot(tplot, Vm_post(:,K_stim_cell_indx),'g','linewidth',0.3); hold on
xlabel('Time(s)'), ylabel('V_m(mV)')


ax = gca;
ax.FontSize = 22;
ax.LineWidth = 3;

ax.Color = 'none';
box off
ax.XTick = [0:20:Tmax];
ax.YTick = [-80:20:-20];
ax.YLim = [-90,-20];
ax.XLim = [0,Tmax];

%% equations

function dXdt = Eqns(t,X,pars)

%% load parameters
tspan = pars.tspan;
tot_cell_number = pars.tot_cell_number;
K_stim = pars.K_stim;
wash_out = pars.wash_out;
current_stim = pars.current_stim;
K_stim_cell_indx = pars.K_stim_cell_indx;
K_o_Stim = pars.K_o_Stim;
K_o_rest = pars.K_o_rest;
K_stim_onset = pars.K_stim_onset;
K_stim_end = pars.K_stim_end;
K_i = pars.K_i;
K_o_washout = pars.K_o_washout;
K_washout_onset = pars.K_washout_onset;
K_washout_end = pars.K_washout_end;
I_stim_onset = pars.I_stim_onset;
I_stim_end = pars.I_stim_end;
current = pars.current;
R = pars.R;
T = pars.T;
F = pars.F;
RT_F = pars.RT_F;
z_K = pars.z_K;
delta_V_kir = pars.delta_V_kir;
G_kir = pars.G_kir;
n_kir = pars.n_kir;
k_kir = pars.k_kir;
E_bg = pars.E_bg;
G_bg = pars.G_bg;
Cm = pars.Cm;
z_Na = pars.z_Na;
z_Ca = pars.z_Ca;
A_m = pars.A_m;
Ca_o = pars.Ca_o;
Na_o = pars.Na_o;
Ca_i = pars.Ca_i;
Na_i = pars.Na_i;
PTRPV4Ca = pars.PTRPV4Ca;
PTRPV4Na = pars.PTRPV4Na;
PTRPV4K = pars.PTRPV4K;
TRPV_P0 = pars.TRPV_P0;
Index_CaP = pars.Index_CaP;
R_gj = pars.R_gj;

%% state variables
X = X(:);
Vm = X';

%% STIMULATION PROTOCOL
K_o = K_o_rest*ones(1,tot_cell_number);
if K_stim
    if t > K_stim_onset*1e3 && t <= K_stim_end*1e3
        K_o(K_stim_cell_indx) = K_o_Stim;
    else
        K_o(K_stim_cell_indx) = K_o_rest;
    end
end
if wash_out
    if t > K_washout_onset*1e3   && t <= K_washout_end*1e3
        K_o = K_o_washout;
    end
end

% Current stimulation
I_stim = zeros(1,tot_cell_number);
if current_stim
    if t >= I_stim_onset*1e3 && t <= I_stim_end*1e3
        I_stim = current;
    else
        I_stim = 0;
    end
end
%% POTASSIUM EQUILLIBRIUM POTENTIAL
E_K(1:tot_cell_number) = RT_F/z_K*log(K_o./K_i);     %[mV]

%% MEMBRANE CURRENTS
I_bg = G_bg.*(Vm - E_bg);        %[pA] lumped background current


I_kir = G_kir.*(K_o).^n_kir .*((Vm - E_K)./(1 + exp((Vm - E_K - delta_V_kir)./k_kir)));   %[pA] whole cell kir current


if Vm ~= 0
    I_TRPV4Na = 1e6.*A_m.*PTRPV4Na.*z_Na.^2.*F.^2./(R.*T).*Vm.*(Na_i-Na_o.*exp(-Vm.*(z_Na.*F)./(R.*T)))./(1-exp(-Vm.*(z_Na.*F)./(R.*T)));
    I_TRPV4Ca = 1e6.*A_m.*PTRPV4Ca.*z_Ca.^2.*F.^2./(R.*T).*Vm.*(Ca_i-Ca_o.*exp(-Vm.*(z_Ca.*F)./(R.*T)))./(1-exp(-Vm.*(z_Ca.*F)./(R.*T)));
    I_TRPV4K =  1e6.*A_m.*PTRPV4K .*z_K.^2 .*F.^2./(R.*T).*Vm.*(K_i-K_o.*exp(-Vm.*(z_K.*F)./(R.*T)))./(1-exp(-Vm.*(z_K.*F)./(R.*T)));
else
    I_TRPV4Na = 1e6.*A_m.*PTRPV4Na.*z_Na.*F.*(Na_i-Na_o);
    I_TRPV4Ca = 1e6.*A_m.*PTRPV4Ca.*z_Ca.*F.*(Ca_i-Ca_o);
    I_TRPV4K =  1e6.*A_m.*PTRPV4K .*z_K.*F.*(K_i-K_o);
end

ind = find(tspan>=t,1);
if ~isempty(ind)
    P0_TPRV4 = TRPV_P0(:,find(tspan>=t,1))';
else
    P0_TPRV4 = rand(size(I_TRPV4K));
end


I_TRPV4 = P0_TPRV4.* (I_TRPV4K + I_TRPV4Na);

% Gap junctional current [pA]
I_gj_temp = zeros(tot_cell_number);
for CaP_i = 1:tot_cell_number
    for CaP_j = 1:tot_cell_number
        if Index_CaP(CaP_i,CaP_j) ~= 0
            I_gj_temp(CaP_i,CaP_j) = 1/R_gj*(Vm(CaP_i)-Vm(CaP_j));
        end
    end
end
I_gj = sum(I_gj_temp, 2);
I_gj = I_gj';

% Total membrane current
I_tot = I_bg + I_kir + I_TRPV4;
% Differential equation
dVmdt = -1./Cm .* (I_tot + I_gj- I_stim);


dXdt = dVmdt(:);

end

function [P0]= TRPV4_P0_3state1channel(tspan,TRPV4_pars)

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







