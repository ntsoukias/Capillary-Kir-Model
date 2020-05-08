clear, clc, close all

%% add necessary folders to the path
addpath(genpath(fullfile(cd ,'..\..')))

%% Simulation time
Tmax = 5;                % [s] Total time of the simulation
tspan = [0, Tmax*1e3];      % converted to [ms]

%% %%%%%%%%%%%%%%%% Stimulation protocol %%%%%%%%%%%%%%%%%%%%%%%%%
K_stim = true;                % extracellular potassium stimulation

%% geometry
s = load('Geometry_Fig6_Fig7B');
G = s.G;

iview = [15 21];
Adj = adjacency(G);
a = largestcomponent(Adj);
G = subgraph(G,a);

%% Options
Options.SensitizeBarrel = 'on'; 
Options.StimulateArtNodes = 'off';
Options.NoCapillaries = 'off';
Options.BC = true;
Options.SplinePlot = true;

%% Set stim node indices

seed =  3028;

xspan = 140;
yspan = 140;
zspan = 140;

stimnode = [];
for i = 1:numel(seed)
    iseed = seed(i);
    inodeList = getnodefromsource(G, iseed, xspan(i), yspan(i), zspan(i));
    stimnode = [stimnode inodeList];
end

if strcmp(Options.SensitizeBarrel,'on')
    seedForBarrel = 1223;
    barrelInd = getnodefromsource(G, seedForBarrel, 110, 110, 900);
else
    barrelInd = [];
end

stimnodeList = zeros(numnodes(G),1);
stimnodeList(stimnode) = 1;
G.Nodes.stimnode = stimnodeList;

barrellnodeList = zeros(numnodes(G),1);
barrellnodeList(barrelInd) = 1;
G.Nodes.barrellnode = barrellnodeList;

figure, H = plotgraph(G);
if strcmp(Options.SensitizeBarrel,'on')
    highlight(H,barrelInd,'MarkerSize',3,'Marker','o','NodeColor','g')
end
highlight(H,stimnode,'MarkerSize',4,'Marker','o','NodeColor','m')

view(iview)

%% geometry

Adj = G.adjacency;
Adj = Adj + Adj';
Adj = full(Adj);
Adj(Adj > 1) = 1;

[~,capInd] = subfromedge(G,find(G.Edges.Type == 0));
[~,artInd] = subfromedge(G,find(G.Edges.Type == 1));
[~,venInd] = subfromedge(G,find(G.Edges.Type == 2));

ArtDilateInds = intersect(barrelInd,artInd);

if strcmp(Options.StimulateArtNodes,'off')
    stimnode(ismember(stimnode, artInd)) = [];
elseif strcmp(Options.NoCapillaries,'on')
    stimnode(ismember(stimnode, capInd)) = [];
end

capInd(ismember(capInd, artInd)) = [];
capInd(ismember(capInd, venInd)) = [];
nnod = numnodes(G);

%% find average length of the edges connected to a node
desiredLCap = 20;
desiredLArt = 33;
desiredL = desiredLArt*ones(1,nnod);
desiredL(capInd) =  desiredLCap;

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
    
    avgL(i) = mean(G.Edges.L(allEdges));
end

circArtECs = 3; %number of Art ECs in the circumferential direction
circVenECs = 3; %number of Vein ECs in the circumferential direction
LengthFactor = avgL./desiredL;
LengthFactor(artInd) = circArtECs*LengthFactor(artInd);
LengthFactor(venInd) = circVenECs*LengthFactor(venInd);

%% Potassium stimulation

G.Nodes.Stim(:) = 0;
G.Nodes.Stim(stimnode) = 1;

K_o_Stim = 10;              % [mM] extracellular potassium concentration for stimulated cells
K_o_rest = 3*ones(1,nnod);                 % [mM] extracellular potassium concentration for unstimulated cells
K_stim_onset = 1;         % [s] K Stimulation onset
K_stim_end = 5;     % [s] K Stimulation end
K_i = 150;                      % [mM] intracellular potassium concentration

%% %%%%%%%%%%%%%%%%%% model parameters %%%%%%%%%%%%%%%%%%%%%%%%

%% constant parameters
R  = 8314.0;	% [mJmol-1K-1]	gas constant
T = 293.0;      % [K] absolute temperature
F  = 96487.0;	% [Cmol-1] Faraday's constant
RT_F = R*T/F;   % RT/F
z_K  = 1;       % K ion valence

%% background current
E_bg = -30;        % [mV] resting membrane potential
G_bg = 0.06*ones(1,nnod);
G_bg(artInd) = 0.2;
G_bg(venInd) = 0.2;

%% Kir channel characteristic
delta_V_kir  = 25;	% [mV]	voltage diff at half-max. eI_K_i
G_kirbar = ones(1,nnod);   % [nS/mM^0.5] inward rectifier constant

ratio = 5;
factorArt = 1;
ratioBarrel = 20;

G_kirbar(capInd) = ratio/sqrt(K_o_rest(1))*G_bg(capInd(1));

if strcmp(Options.SensitizeBarrel,'on')
    G_kirbar(barrelInd) = ratioBarrel/sqrt(K_o_rest(1))*G_bg(capInd(1));
end
G_kirbar(artInd) = factorArt*0.6/sqrt(3);   % [nS/mM^0.5] inward rectifier constant

G_bg = G_bg.*LengthFactor;

G_kirbar(venInd) = 0;
G_kirbar = G_kirbar.*LengthFactor;

n_kir = 0.5;	% inward rectifier constant
k_kir = 7;	% [mV]	inward rectifier slope factor

%% Gap junctional resistance

Rgj = 400e-3;
R_gj = Rgj*Adj;  % [Gohm]  gap junctional resistance
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
        R_gj(i,j) = 10e-3;
        R_gj(i,j) = R_gj(i,j)*fact;
    elseif itype == 1
        fact = ilength/desiredLArt/circVenECs;
        R_gj(i,j) = 10e-3;
        R_gj(i,j) = R_gj(i,j)*fact;
    elseif itype == 3
        fact = ilength/desiredLArt/circpArtECs;
        R_gj(i,j) = 10e-3;
        R_gj(i,j) = R_gj(i,j)*fact;
    elseif itype == 4
        fact = ilength/desiredLArt/circpVenECs;
        R_gj(i,j) = 10e-3;
        R_gj(i,j) = R_gj(i,j)*fact;
    else
        fact = ilength/desiredLCap;
        R_gj(i,j) = R_gj(i,j)*fact;
    end
    cond1 = ismember(i,capInd) & ismember(j,ArtDilateInds);
    cond2 = ismember(i,ArtDilateInds) & ismember(j,capInd);
    if cond1 | cond2
        R_gj(i,j) = 10e-3 *ilength/desiredLCap;
    end
    R_gj(j,i) = R_gj(i,j);
end

%% boundary conditions
disp('Boundary Conds...')
Gin = zeros(1,nnod);
BoundaryNodes = find(G.degree == 1);
P.VmBC = -35; %mV voltage of the boundary nodes

Rgj_cap = Rgj;
G_bg_cap = G_bg(capInd(1));
Rgj_PA = 10e-3;
G_bg_PA = G_bg(artInd(1));

tic
if Options.BC
    for node = 1:numel(BoundaryNodes)
        i = BoundaryNodes(node);
        if ismember(i, capInd)
            Rin = sqrt(Rgj_cap/G_bg_cap);
            Gin(i) = 1/Rin;
        elseif ismember(i, artInd) | ismember(i, venInd)
            Rin = 1/circArtECs*sqrt(Rgj_PA/G_bg_PA);
            Gin(i) = 1/Rin;
        elseif ismember(i, circpArtECs) | ismember(i, circpVenECs)
            Rin = 1/circpArtECs*sqrt(Rgj_PA/G_bg_pial);
            Gin(i) = 1/Rin;
        end
    end
end
toc
P.Gin=Gin;

%% pattern of the jacobian matrix
S = Adj;
S = S + eye(size(S));
S = sparse(S);

%% Cell capacitance
Cm = 16*ones(1,nnod);   % [pF] membrane capacitance
Cm(capInd) = 8;
Cm = Cm.*LengthFactor;

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
    eval(['P.' ParNames{i} '=' ParNames{i} ';' ]);
end

%% solve
disp('solving...')
Vm0 = -30*ones(1,nnod);      %[mV] initial condition for membrane potentials
opts=odeset('OutputFcn',@odeprog,'JPattern',S);
[t,Vm_post] = ode15s(@(t,X)equations(t,X,P),tspan,Vm0,opts);

%% Network colormapped

G0 = G;

tplot = t/1000;     % convert to seconds
tInstance = 5;
tInd = find(tplot>=tInstance,1,'first');

nodeval = Vm_post(tInd,:);
edgeval = getedgeval(G0,nodeval);
climits = [-80,-20];

G0.Edges.Vm = edgeval;

G1 = subfromedge(G0,find(G0.Edges.Type > 0));
G0 = rmedge(G0,find(G0.Edges.Type > 0 ));

edgeval = G0.Edges.Vm;
edgeval1 = G1.Edges.Vm;

diameter = normalize(G0.Edges.D,'range',[1 8]);
diameter(G0.Edges.Type == 1) = 5;
diameter(G0.Edges.Type == 2) = 5;

diameter1 = normalize(G1.Edges.D,'range',[1 8]);
diameter1(G1.Edges.Type == 1) = 3.5;
diameter(G0.Edges.Type == 0) = 1;
diameter1(G1.Edges.Type == 2) = 2.5;

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
fig = figure;
plot(tplot, Vm_post(:,capInd(1:5:end)),'color',[0.7 0.7 0.7],'linewidth',2), hold on
plot(tplot, Vm_post(:,artInd),'r','linewidth',1), hold on
plot(tplot, Vm_post(:,venInd),'b')
xlabel('t(s)'), ylabel('V_m(mV)')
ax = gca;
ax.LineWidth = 1;
ax.FontSize = 16;
ax.FontName = 'arial';
ax.FontWeight = 'bold';
box off
fig.Color = 'w';

%% plot using cubic spline
if Options.SplinePlot
    disp('Plotting....')
    tpoints = [1:0.5:4];
    iview = [24,-2];
    
    tic
    for kk = 1:numel(tpoints)
        tpoint = tpoints(kk);
        tInd = find(tplot>=tpoint,1,'first');
        
        nodeval = Vm_post(tInd,:);
        edgeval = getedgeval(G,nodeval);
        Q = edgeval;
        eval(sprintf('G.Edges.Q%d = Q;',kk))
    end
    
    
    for kk = 1:numel(tpoints)
        
        tpoint = tpoints(kk);
        
        valname = sprintf('Q%d',kk);
        if kk == 1
            out = val_colored(valname, G, iview);
        else
            
            out.H1.Edges.ToColor = out.H1.Edges.(valname);
            colorcanvas1(out.H1, out.all_hs1, out.allEdgeIDs1)
            
            out.H2.Edges.ToColor = out.H2.Edges.(valname);
            colorcanvas1(out.H2, out.all_hs2, out.allEdgeIDs2)

        end
        
        caxis([-80,-15])
        
        ax = gca;
        ax.Title.String = [sprintf('%1.2f seconds',tpoint - 1)];
        ax.Title.FontSize = 16;
        ax.Title.FontWeight = 'bold';
        
        if kk == numel(tpoints)
            view(0,90)
        end
        
        set(gcf, 'InvertHardcopy', 'off')
        pause(2)
        
    end
    toc
end

%% functions

function dVmdt = equations(t,X,P)

X = X(:);
Vm = X';

%% parameters
RT_F = P.RT_F;
z_K = P.z_K;
nnod = P.nnod;
stimnode = P.stimnode;
K_stim = P.K_stim;
K_stim_onset = P.K_stim_onset;
K_stim_end = P.K_stim_end;
Cm = P.Cm;
G_bg = P.G_bg;
R_gj = P.R_gj;
E_bg = P.E_bg;
delta_V_kir = P.delta_V_kir;
G_kirbar = P.G_kirbar;
n_kir = P.n_kir;
k_kir = P.k_kir;
K_i = P.K_i;
K_o_Stim = P.K_o_Stim;
K_o_rest = P.K_o_rest;
neighbors = P.neighbors;
Gin = P.Gin;
VmBC = P.VmBC;

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
        I_gj(i) = I_gj(i) + 1/R_gj(i,neighbors(i,j))*(Vm(i) - Vm(neighbors(i,j))) +  Gin(i)*(Vm(i) - (VmBC));
    end
end


%% Total membrane current
I_tot = I_bg + I_kir;

%% Differential equation
dVmdt =  -1./Cm .* (I_tot + I_gj);
dVmdt = dVmdt';

end

