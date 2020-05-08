%% stimulation protocol
Options.SensitizeBarrel = 'off';     % use different ratio for barrel nodes
Options.StimulateArtNodes = StimulateArtNodesCond;
Options.StimulateCapillaries = StimulateCapillariesCond;

seedForBarrel = 8404 ; 
xbarrel=100; ybarrel=100;  zbarrel= 1000;

%% Parameters

Rgj_cap = CapRgj;   % [Gohm] cEC Rgj
Rgj_PA=  10e-3; % [Gohm] art/venul single cell Rgj

E_bg = -30;        % [mV] resting membrane potential
G_bg_cap = 0.06;    %nS
G_bg_PA = 0.2;  %nS

G_kirbar_cap   =   cECGkirbarRatio*G_bg_cap ;
G_kirbar_PA    =  PAGkirbarRatio*0.6/sqrt(3);

ratioBarrel = 14;    % used if SensitizeBarrel is on for GKirbar cap

%% Set stim node indices

stimnode = [];
for i = 1:numel(seed)
    iseed = seed(i);
    inodeList = getnodefromsource(G, iseed, xspan, yspan, zspan);
    stimnode = [stimnode inodeList];
end

if strcmp(Options.SensitizeBarrel,'on')
    barrelInd = getnodefromsource(G, seedForBarrel, xbarrel, ybarrel, zbarrel);
else
    barrelInd = [];
end

stimnode = stimnode(randperm(numel(stimnode), round(PercentNodes*numel(stimnode)/100)));


if strcmp(Options.StimulateArtNodes,'off')
    stimnode(ismember(stimnode, artInd)) = [];
end

if strcmp(Options.StimulateCapillaries,'off')
    stimnode(ismember(stimnode, capInd)) = [];
end

stimnode(ismember(stimnode, venInd)) = [];

G.Nodes.Stim(:) = 0;
G.Nodes.Stim(stimnode) = 1;


%% background current
G_bg = G_bg_cap *ones(1,nnod);  %nS
G_bg(artInd) = G_bg_PA; 
G_bg(venInd) = G_bg_PA;

%% Kir channel characteristic
K_o_rest = 3*ones(1,nnod);                 % [mM] extracellular potassium concentration for unstimulated cells
delta_V_kir  = 25;    % [mV]	voltage diff at half-max. eI_K_i 
G_kirbar = zeros(1,nnod);   % [nS/mM^0.5] inward rectifier constant
G_kirbar(capInd) = G_kirbar_cap;      

if strcmp(Options.SensitizeBarrel,'on')
    G_kirbar(barrelInd) = ratioBarrel*G_bg_cap;
end

G_kirbar(artInd) = G_kirbar_PA;   % [nS/mM^0.5] inward rectifier constant
G_bg = G_bg.*LengthFactor;

G_kirbar(venInd) = 0;
G_kirbar = G_kirbar.*LengthFactor;

n_kir = 0.5;	% inward rectifier constant
k_kir = 7.0;	% [mV]	inward rectifier slope factor

%% find Rgj

R_gj = Adj;  % [Gohm]  gap junctional resistance

% disp('Calculating Rgj...')
edgelist = G.Edges{:,1};

for k = 1:size(edgelist,1)
    nodes = edgelist(k,:);
    i = nodes(1);
    j = nodes(2);
    edgeInd = k;
    ilength = G.Edges.L(edgeInd)/desiredLCap;
    itype = G.Edges.Type(edgeInd);
    
    if itype == 2
        fact = ilength/circArtECs;
        R_gj(i,j) = Rgj_PA*fact;
    elseif itype == 1
        fact = ilength/circVenECs;
        R_gj(i,j) = Rgj_PA*fact;
    else
        fact = ilength;
        R_gj(i,j) = CapRgj*fact;
    end
    
    cond1 = ismember(i,capInd) & ismember(j,artInd);
    cond2 = ismember(i,artInd) & ismember(j,capInd);
    if cond1 | cond2
        R_gj(i,j) = Rgj_PA *ilength;
    end
    R_gj(j,i) = R_gj(i,j);
end

Param.R_gj = R_gj;

%% implement boundary conditions
% disp('Boundary Conds...')
Gin = zeros(1,nnod);
BoundaryNodes = find(G.degree == 1);
Param.VmBC = -30; %mV voltage of the boundary nodes

if BoundaryConds
    for node = 1:numel(BoundaryNodes)
        i = BoundaryNodes(node);
        if ismember(i, capInd)
            Rin = sqrt(Rgj_cap/G_bg_cap);
        elseif ismember(i, artInd) | ismember(i, venInd)
            Rin = 1/circArtECs*sqrt(Rgj_PA/G_bg_PA);
        elseif ismember(i, circpArtECs) | ismember(i, circpVenECs)
            Rin = 1/circpArtECs*sqrt(Rgj_PA/G_bg_pial);
        end
        Gin(i) = 1/Rin;
    end
end
    
%% pattern of the jacobian matrix
S = Adj;
S = S + eye(size(S));
S = sparse(S);

%% Cell capacitance
Cm = 16*ones(1,nnod);   % [pF] membrane capacitance
Cm(capInd) = 8;
Cm = Cm.*LengthFactor;

%% find neighbors
Param.edgelist = G.Edges{:,1};
neighbors = zeros(nnod, max(G.degree));

for i = 1:nnod
    neighborNodes = G.neighbors(i);
    neighborNodes = neighborNodes';
    neighbors(i,1:numel(neighborNodes)) = neighborNodes;
end

%% store pars
ParNames = who();
for i = 1:numel(ParNames)
    eval(['Param.' ParNames{i} '=' ParNames{i} ';' ]);
end

%% solve
Vm0 = E_bg*ones(1,nnod);      %[mV] initial condition for membrane potentials
opts=odeset('JPattern',S);
[t,Vm_post] = ode15s(@(t,X)equations(t,X,Param,Options),tspan,Vm0,opts);

%% main function
function dVmdt = equations(t,X,Param, Options)
X = X(:);
Vm = X';

%% options
K_stim = Options.K_stim;

%% parameters
RT_F = Param.RT_F;
z_K = Param.z_K;
nnod = Param.nnod;
stimnode = Param.stimnode;
K_stim_onset = Param.K_stim_onset;
K_stim_end = Param.K_stim_end;
Cm = Param.Cm;
G_bg = Param.G_bg;
R_gj = Param.R_gj;
E_bg = Param.E_bg;
delta_V_kir = Param.delta_V_kir;
G_kirbar = Param.G_kirbar;
n_kir = Param.n_kir;
k_kir = Param.k_kir;
K_i = Param.K_i;
K_o_Stim = Param.K_o_Stim;
K_o_rest = Param.K_o_rest;
neighbors = Param.neighbors;
Gin=Param.Gin;
VmBC = Param.VmBC; 

%% Stimulation protocol
% Potassium stimulus
K_o = K_o_rest;

if K_stim
    if t > K_stim_onset*1e3 && t <= K_stim_end*1e3
        K_o(stimnode) = K_o_Stim;
    end
end

%% Potassium equilibrium potential
E_K= RT_F/z_K*log(K_o./K_i);     %[mV]

%% Membrane currents [pA]
I_bg = G_bg.*(Vm - E_bg);        %[pA] lumped background current
I_kir = G_kirbar.*(K_o).^n_kir .*((Vm - E_K)./...
    (1 + exp((Vm - E_K - delta_V_kir)/k_kir)));   %[pA] whole cell kir current

%% Gap junctional current [pA]

I_gj = zeros(1,nnod);
for i = 1:nnod    
    for j = 1:nnz(neighbors(i,:))
            I_gj(i) = I_gj(i) + 1/R_gj(i,neighbors(i,j))*(Vm(i) - Vm(neighbors(i,j))) +  Gin(i)*(Vm(i) - (VmBC));      
    end    
end

%% Total membrane current
I_tot = I_bg + I_kir;

%% Differential equation
dVmdt =  -1./Cm .* (I_tot + I_gj);
dVmdt = dVmdt';

end
