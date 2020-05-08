clear, clc, close all

%% add necessary folders to the path
addpath(genpath(fullfile(cd ,'..\..')))

%% Simulation time
Tmax = 800; % [s] Total time of the simulation
dt = 0.0001;  %[s] time increment
tspan = 0:dt*1e3:Tmax*1e3;      % converted to [ms]

%% %%%%%%%%%%%%%%%% Stimulation protocol %%%%%%%%%%%%%%%%%%%%%%%%%
K_stim = true;                % extracellular potassium stimulation
current_stim = false;     % current stimulus

%% Potassium stimulation
K_o_rest = 3;                 % [mM] extracellular potassium concentration for unstimulated cells
K_stim_onset = 70;      % [s] K Stimulation onset
K_stim_end = 230;     % [s] K Stimulation end
K_i = 150;                   % [mM] intracellular potassium concentration

K_out = @(t) 3*heaviside(t+1e-3) + 2*heaviside(t-200e3) + ...
    3*heaviside(t-400e3) - 5*heaviside(t-600e3);

%% Current stimulation
I_stim_onset = 5;           % [s] current stimulation onset
I_stim_end = 10;            % [s] current stimulation end
current = -1;        % [pA]   Injected current

%% %%%%%%%%%%%%%%%%%% model parameters %%%%%%%%%%%%%%%%%%%%%%%%

%% constant parameters
R  = 8314.0;	% [mJmol-1K-1]	gas constant
T = 293;  % [K] absolute temperature
F  = 96485.0;	% [Cmol-1] Faraday's constant
RT_F = R*T/F;   % RT/F
z_K  = 1;       % K ion valence

%% Kir channel characteristic
delta_V_kir  = 25;% [mV]	voltage diff at half-max. eI_K_i
G_kirbar = 0.18;   % [nS/mM^0.5] inward rectifier constant

n_kir = 0.5;	% inward rectifier constant (0.5)
k_kir = 7;	% [mV]	inward rectifier slope factor

%% Background current
E_bg = -30;        % [mV] resting membrane potential
G_bg = 0.06; % [nS]  lumped background conductance

%% Cell capacitance
Cm = 8;   % [pF] capillary membrane capacitance

%% TRPV4 current characteristics (stochastic)

% parameters
z_Na = 1	;	%;		Na ion valence
z_Ca = 2	;	%Ca ion valence
A_m = Cm./1e6;	%[cm2]	membrane area
Ca_o = 2.0		;%[mM]	[Ca2+]o extracellular calcium
Na_o = 134.0	; %[mM]	[Na+]o  extracellular sodium
Ca_i = 0.2; 	%[mM]	[Ca2+]i intracellular calcium
Na_i =  10 ;	%[mM]	[Na+]i  intracellular sodium
PCa= 1e-7; %[cm./sec]
PTRPV4Ca= PCa;   %[cm./sec]
PTRPV4Na=PCa./6.9;
PTRPV4K=1.35.*PCa./6.9;

%% TRPV4 current

% --------------------------------- %'3 states 1 channel'
disp('Generating TRPV4 traces ....')
meanshutlifetime= 300e3;%msec
meanopenlifetime = 90;%msec
meanburstlength = 20e3;%msec
meanblockedlifetime = 33; %msec
TRPV4_pars = [meanshutlifetime, meanopenlifetime, ...
    meanburstlength, meanblockedlifetime];

TRPV_P0 = TRPV4_P0_3state1channel(tspan,TRPV4_pars);

%% store all parameters
ParNames = who();
for i = 1:numel(ParNames)
    eval(['pars.' ParNames{i} '=' ParNames{i} ';' ]);
end

%% Solution
% initial condition
Vm0 = -30;      %[mV] initial condition for membrane potentials
init = [Vm0];
opt = optimset;
opt.Display = 'off';
X0 = fsolve(@(X) Equations(0,X,pars), init, opt);

% Solve

tic
disp('Solving....')
opts=odeset('OutputFcn',@odeprog);
[t,X] = ode113(@(t,X) Equations(t,X,pars),[0 Tmax*1e3],X0,opts);
toc

Vm	= X;

%% Plot Vm vs. Time
Fig = figure;
Fig.Renderer = 'painters';

tplot = t/1000;     % convert to seconds
h1 = plot(tplot/60, Vm,'k','linewidth',2);

ax = gca;
ax.FontSize = 22;
ax.FontName = 'arial';
ax.LineWidth = 3;
xlabel('t(min)'), ylabel('V_m (mV)')
ax = gca; ax.Color = 'none';
ylim([-90,-10])
ax.YTick = [-80:20:-20];
ax.XTick = [0:5:15];
box off

%% Equations
function dXdt = Equations(t,X,pars)

tspan = pars.tspan;
K_i = pars.K_i;
current_stim = pars.current_stim;
I_stim_onset = pars.I_stim_onset;
I_stim_end = pars.I_stim_end;
current = pars.current;
R = pars.R;
T = pars.T;
F = pars.F;
RT_F = pars.RT_F;
z_K = pars.z_K;
delta_V_kir = pars.delta_V_kir;
G_kirbar = pars.G_kirbar;
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
K_out = pars.K_out;

%% state variables
X = X(:);
Vm = X;

%% STIMULATION PROTOCOL
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
%% POTASSIUM EQUILLIBRIUM POTENTIAL
E_K = RT_F/z_K*log(K_o./K_i);     %[mV]

%% MEMBRANE CURRENTS
I_bg = G_bg.*(Vm - E_bg);        %[pA] lumped background current
I_kir = G_kirbar.*(K_o).^n_kir .*((Vm - E_K)./(1 + exp((Vm - E_K - delta_V_kir)./k_kir)));   %[pA] whole cell kir current


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

% Total membrane current
I_tot = I_bg + I_kir + I_TRPV4;
% Differential equation
dVmdt = -1./Cm.* (I_tot - I_stim);

dXdt= dVmdt;

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


