clear, clc, close all

%% add necessary folders to the path
addpath(genpath(fullfile(cd ,'..\..')))

%% Simulation time
Tmax = 20;                % [s] Total time of the simulation
tspan = [0 Tmax*1e3];      % converted to [ms]

splinePlot = true;

%% stimulation protocol
Options.K_stim = true;                % extracellular potassium stimulation
Options.SensitizeBarrel = 'off';     % use different ratio for barrel nodes
Options.StimulateArtNodes = 'off';
Options.StimulateCapNodes = 'on';
Options.BoundaryConds = true;

%% stim region
seed = [300 1200 -950];

xspan = 120;
yspan = 120;
zspan = 150;

PercentNodes = 50;  %percentage of capillaries to be stimulated in the stimulated region

seedForBarrel = 8404 ; 
xbarrel=100; ybarrel=100;  zbarrel= 1000;

%% Parameters

Rgj_cap = 50e-3;  % [Gohm] cEC Rgj
Rgj_PA=  10e-3;  % [Gohm] art/venul  Rgj
Cm_cap = 8;   % [pF] capillary membrane capacitance
Cm_PA= 16;
E_bg = -30;        % [mV] resting membrane potential
G_bg_cap = 0.06; %nS
G_bg_PA = 0.2;   %nS
G_bg_pial = 0.3;   %nS

G_Kir_cap   =   10/sqrt(3)*G_bg_cap ;   
G_Kir_PA    =   3*G_bg_PA;
G_Kir_pial  =   3*G_bg_pial;

%% geometry

disp('Loading Network...')
s = load('Geometry_Fig7A');


Adj = adjacency(s.G);
a = largestcomponent(Adj);
G = subgraph(s.G,a);

Adj = G.adjacency;
Adj = Adj + Adj';
Adj = full(Adj);
Adj(Adj > 1) = 1;

[~,capInd] = subfromedge(G,find(G.Edges.Type == 0));
[~,artInd] = subfromedge(G,find(G.Edges.Type == 1));
[~,venInd] = subfromedge(G,find(G.Edges.Type == 2));
[~,partInd] = subfromedge(G,find(G.Edges.Type == 3));
[~,pvenInd] = subfromedge(G,find(G.Edges.Type == 4));

capInd(ismember(capInd, artInd)) = [];
capInd(ismember(capInd, venInd)) = [];

nnod = numnodes(G);

desiredLCap = 20;
desiredLArt = 33;
desiredL = desiredLArt*ones(1,nnod);
desiredL(capInd) =  desiredLCap;
ratioBarrel = 14;

%% Potassium stimulation
K_o_Stim = 10;              % [mM] extracellular potassium concentration for stimulated cells
K_o_rest = 3*ones(1,nnod);                 % [mM] extracellular potassium concentration for unstimulated cells
K_stim_onset = 1;         % [s] K Stimulation onset
K_stim_end = 15;     % [s] K Stimulation end
K_i = 150;                      % [mM] intracellular potassium concentration

%% Set stim node indices

stimnode = [];
for i = 1:size(seed,1)
    iseed = seed(i,:);
    inodeList = getnodefromsourceXYZ(G, iseed, xspan, yspan, zspan);
    stimnode = [stimnode inodeList];
end
stimregion = stimnode;


if strcmp(Options.SensitizeBarrel,'on')
    barrelInd = getnodefromsource(G, seedForBarrel, xbarrel, ybarrel, zbarrel);
else
    barrelInd = [];
end

stimnode = stimnode(randperm(numel(stimnode), round(PercentNodes*numel(stimnode)/100)));


if strcmp(Options.StimulateArtNodes,'off')
    stimnode(ismember(stimnode, artInd)) = [];
end
if strcmp(Options.StimulateCapNodes,'off')
    stimnode(ismember(stimnode, capInd)) = [];
end

G.Nodes.Stim(:) = 0;
G.Nodes.Stim(stimnode) = 1;

%% plot network
disp('Plotting Network...')

figure, H = plotgraph(G);
if strcmp(Options.SensitizeBarrel,'on')
    highlight(H,barrelInd,'MarkerSize',3,'Marker','o','NodeColor','g')
end
highlight(H,stimnode,'MarkerSize',3,'Marker','o','NodeColor','m')
iview = [ -104.5875   20.0491];
view(iview)

%% find average length
avgL = zeros(1,nnod);

for i = 1:nnod
    
    if isa(G,'graph')
        allEdges = outedges(G,i);
    elseif isa(G,'digraph')
        outEdges = G.outedges(i);
        inEdges = G.inedges(i);
        allEdges = [outEdges; inEdges];
    else
        error('Invalid graph!')
    end
    
    avgL(i) = sum(G.Edges.L(allEdges))/2;
end

circArtECs = 3; %number of Art ECs in the circumferential direction
circVenECs = 3; %number of Vein ECs in the circumferential direction
circpArtECs = 5; %number of pial Art ECs in the circumferential direction
circpVenECs = 5; %number of pial Vein ECs in the circumferential direction
LengthFactor = avgL./desiredL;
LengthFactor(artInd) = circArtECs*LengthFactor(artInd);
LengthFactor(venInd) = circVenECs*LengthFactor(venInd);
LengthFactor(partInd) = circpArtECs*LengthFactor(partInd);
LengthFactor(pvenInd) = circpVenECs*LengthFactor(pvenInd);

stimPercentage = sum(LengthFactor(stimnode));

%% Model parameters
R  = 8314.0;	% [mJmol-1K-1]	gas constant
T = 293.0;      % [K] absolute temperature
F  = 96487.0;	% [Cmol-1] Faraday's constant
RT_F = R*T/F;   % RT/F
z_K  = 1;       % K ion valence

%% background current
G_bg = G_bg_cap *ones(1,nnod);
G_bg(artInd) = G_bg_PA;   
G_bg(venInd) = G_bg_PA;   
G_bg(partInd) = G_bg_pial;  
G_bg(pvenInd) = G_bg_pial;  

%% Kir channel characteristic
delta_V_kir  = 25;    % [mV]	voltage diff at half-max. eI_K_i 
G_kirbar = ones(1,nnod);   % [nS/mM^0.5] inward rectifier constant

G_kirbar(capInd) = G_Kir_cap;      

if strcmp(Options.SensitizeBarrel,'on')
    G_kirbar(Param.barrelInd) = ratioBarrel*G_bg_cap;
end

G_kirbar(artInd) = G_Kir_PA;   
G_kirbar(partInd) = G_Kir_pial;
G_bg = G_bg.*LengthFactor;

G_kirbar(venInd) = 0;
G_kirbar(pvenInd) = 0;
G_kirbar = G_kirbar.*LengthFactor;

n_kir = 0.5;	% inward rectifier constant
k_kir = 7.0;	% [mV]	inward rectifier slope factor

%% find Rgj

R_gj = Rgj_cap *Adj;  % [Gohm]  gap junctional resistance
disp('Calculating Rgj...')
edgelist = G.Edges{:,1};

tic
for k = 1:size(edgelist,1)
    nodes = edgelist(k,:);
    i = nodes(1);
    j = nodes(2);
    edgeInd = k;
    ilength = G.Edges.L(edgeInd);
    itype = G.Edges.Type(edgeInd);
    
    if itype == 2
        fact = ilength/desiredLArt/circArtECs;
        R_gj(i,j) = Rgj_PA;
        R_gj(i,j) = R_gj(i,j)*fact;
    elseif itype == 1
        fact = ilength/desiredLArt/circVenECs;
        R_gj(i,j) = Rgj_PA ;
        R_gj(i,j) = R_gj(i,j)*fact;
    elseif itype == 3
        fact = ilength/desiredLArt/circpArtECs;
        R_gj(i,j) = Rgj_PA ;
        R_gj(i,j) = R_gj(i,j)*fact;
    elseif itype == 4
        fact = ilength/desiredLArt/circpVenECs;
        R_gj(i,j) = Rgj_PA;
        R_gj(i,j) = R_gj(i,j)*fact;
    else
        fact = ilength/desiredLCap;
        R_gj(i,j) = R_gj(i,j)*fact;
    end
    
    cond1 = ismember(i,capInd) & ismember(j,artInd);
    cond2 = ismember(i,artInd) & ismember(j,capInd);
    if cond1 | cond2
        R_gj(i,j) = Rgj_PA *ilength/desiredLCap;
    end
    
    R_gj(j,i) = R_gj(i,j);
end


toc

%% implement boundary conditions
disp('Boundary Conds...')
Gin = zeros(1,nnod);
BoundaryNodes = find(G.degree == 1);
VmBC = -30; %mV voltage of the boundary nodes

tic
if Options.BoundaryConds
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
toc

%% pattern of the jacobian matrix
S = Adj;
S = S + eye(size(S));
S = sparse(S);

%% Cell capacitance
Cm = 16*ones(1,nnod);   % [pF] membrane capacitance
Cm(capInd) = 8;
Cm = Cm.*LengthFactor;

%% find neighbors
edgelist = G.Edges{:,1};
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
disp('Solving...')
tic
opts=odeset('OutputFcn',@odeprog,'JPattern',S);
[t,Vm_post] = ode15s(@(t,X)equations(t,X,Param,Options),tspan,Vm0,opts);
toc

%% Network colormapped

G0 = G;

tplot = t/1000;     % convert to seconds
tInstance = 14;
tInd = find(tplot>=tInstance,1,'first');

nodeval = Vm_post(tInd,:);
edgeval = getedgeval(G0,nodeval);
climits = [-60,-25];

G0.Edges.Vm = edgeval;

G1 = subfromedge(G0,find(G0.Edges.Type > 0));
G0 = rmedge(G0,find(G0.Edges.Type > 0 ));

edgeval = G0.Edges.Vm;
edgeval1 = G1.Edges.Vm;

diameter = normalize(G0.Edges.D,'range',[1 8]);
diameter(G0.Edges.Type == 1) = 5;
diameter(G0.Edges.Type == 2) = 5;
diameter(G0.Edges.Type == 3) = 6;
diameter(G0.Edges.Type == 4) = 6;


diameter1 = normalize(G1.Edges.D,'range',[1 8]);
diameter1(G1.Edges.Type == 1) = 3.5;
diameter(G0.Edges.Type == 0) = 1;
diameter1(G1.Edges.Type == 2) = 2.5;
diameter1(G1.Edges.Type == 3) = 4;
diameter1(G1.Edges.Type == 4) = 4;

figure
hold on
H = plotcolor(G0,edgeval,'V_m(mV)',climits,diameter);
H.EdgeAlpha = 0.4;

H = plotcolor(G1,edgeval1,'V_m(mV)',climits,diameter1);
H.EdgeAlpha = 1;

axis equal
view(iview)

%% Voltage traces
[~,capInd] = subfromedge(G,find(G.Edges.Type == 0));
[~,artInd] = subfromedge(G,find(G.Edges.Type == 1));
[~,venInd] = subfromedge(G,find(G.Edges.Type == 2));
figure
plot(tplot, Vm_post(:,capInd(1:5:end)),'color',[0.7 0.7 0.7],'linewidth',2), hold on
plot(tplot, Vm_post(:,artInd),'r','linewidth',1), hold on
plot(tplot, Vm_post(:,venInd),'b')
xlabel('t(s)'), ylabel('V_m(mV)')
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 16;
ax.FontName = 'arial';
ax.FontWeight = 'bold';
ax.Color = 'none';
box off

%% plot using cubic spline

tpoints = 4;

if splinePlot
    disp('Plotting......')
    for kk = 1:numel(tpoints)
        
        tic
        tpoint = tpoints(kk);
        tInd = find(tplot>=tpoint,1,'first');
        
        nodeval = Vm_post(tInd,:);
        edgeval = getedgeval(G,nodeval);
        Q = edgeval;
        G.Edges.Q = Q;
        val_colored('Q',G,iview);
        caxis([-60,-25])
        
        toc
        set(gcf, 'InvertHardcopy', 'off')
        
    end
    
end

%% main function
function dVmdt = equations(t,X,Param, Options)
X = X(:);
Vm = X';

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

%% --------------------------------------------------- Stimulation protocol
% Potassium stimulus
if K_stim
    if t > K_stim_onset*1e3 && t <= K_stim_end*1e3
        K_o_rest(stimnode) = K_o_Stim;
    end
end

%% Potassium equilibrium potential
E_K= RT_F/z_K*log(K_o_rest./K_i);     %[mV]

%% Membrane currents [pA]
I_bg = G_bg.*(Vm - E_bg);        %[pA] lumped background current
I_kir = G_kirbar.*(K_o_rest).^n_kir .*((Vm - E_K)./...
    (1 + exp((Vm - E_K - delta_V_kir)/k_kir)));   %[pA] whole cell kir current

%% Gap junctional current [pA]

I_gj = zeros(1,nnod);
for i = 1:nnod
    
    for j = 1:nnz(neighbors(i,:))
        I_gj(i) = I_gj(i) + 1/R_gj(i,neighbors(i,j))*(Vm(i) - Vm(neighbors(i,j))) +  Gin(i)*(Vm(i) - VmBC);
    end
    
end

%% Total membrane current
I_tot = I_bg + I_kir;

%% Differential equation
dVmdt =  -1./Cm .* (I_tot + I_gj);
dVmdt = dVmdt';

end






    
    



