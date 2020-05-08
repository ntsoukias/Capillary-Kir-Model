clear, clc, close all

%% add necessary folders to the path
addpath(genpath(fullfile(cd ,'..\..')))

%% Simulation time
Tmax = 8;                % [s] Total time of the simulation
tspan = [0 Tmax*1e3];      % converted to [ms]
BoundaryConds = true;

%% Potassium stimulation
Options.K_stim = true;                % extracellular potassium stimulation
K_o_Stim = 10;              % [mM] extracellular potassium concentration for stimulated cells
K_stim_onset = 1;         % [s] K Stimulation onset
K_stim_end = 5;     % [s] K Stimulation end
K_i = 150;                      % [mM] intracellular potassium concentration

%% geometry

s = load('Geometry_Fig6_Fig7B');
G = s.G;

iview = [15 21];
Adj = adjacency(G);
a = largestcomponent(Adj);
G = subgraph(G,a);

Adj = G.adjacency;
Adj = Adj + Adj';
Adj = full(Adj);
Adj(Adj > 1) = 1;

[~,capInd] = subfromedge(G,find(G.Edges.Type == 0));
[~,artInd] = subfromedge(G,find(G.Edges.Type == 1));
[~,venInd] = subfromedge(G,find(G.Edges.Type == 2));

capInd(ismember(capInd, artInd)) = [];
capInd(ismember(capInd, venInd)) = [];

nnod = numnodes(G);

%% find main art nodes
seed = 3028;
xspan = 100; yspan = 100; zspan = 1500;
NodeList = getnodefromsource(G, seed, xspan, yspan, zspan);
MainArtNodes = intersect(NodeList, artInd);

TopInd = find(G.Nodes.Z(MainArtNodes) == max(G.Nodes.Z(MainArtNodes) ));
TopArtNode = MainArtNodes(TopInd);

%% stimulation region
xspan = 100;    %um
yspan = 100;
zspan = 100;

PercentNodes = 100;  %percentage of capillaries to be stimulated in the stimulated region

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

desiredLCap = 20;   %um desired length of capillaries
desiredLArt = 30;   %um desired length of arterioles
desiredL = desiredLCap*ones(1,nnod);       % desired length of cells (um)
desiredL(artInd) = desiredLArt;
desiredL(venInd) = desiredLArt;
circArtECs = 3; %number of Art ECs in the circumferential direction
circVenECs = 3; %number of Vein ECs in the circumferential direction
LengthFactor = avgL./desiredL;
LengthFactor(artInd) = circArtECs*LengthFactor(artInd);
LengthFactor(venInd) = circVenECs*LengthFactor(venInd);

%% parameters
R  = 8314.0;	% [mJmol-1K-1]	gas constant
T = 293.0;      % [K] absolute temperature
F  = 96487.0;	% [Cmol-1] Faraday's constant
RT_F = R*T/F;   % RT/F
z_K  = 1;       % K ion valence

%% parameters of interest

CapRgjVals = 100e-3;  %cEC Rgj [Gohm]
Seeds = 1359;   % layer IV vs layer VI stimulus
PAGkirbarRatioVals = 4; % increasing PA Gkirbar by a factor
cECGkirbarRatioVals = 5/sqrt(3);    % cEC Gkirbar
StimulateArtNodesConds = {'on'};
StimulateCapillariesConds = {'off'};

Vars = {CapRgjVals, Seeds, PAGkirbarRatioVals, cECGkirbarRatioVals,...
    StimulateArtNodesConds, StimulateCapillariesConds};

AllCombos = combvec(1:numel(CapRgjVals), 1:numel(Seeds),...
    1:numel(PAGkirbarRatioVals), 1:numel(cECGkirbarRatioVals),...
    1:numel(StimulateArtNodesConds), 1:numel(StimulateCapillariesConds));

summaryData = zeros(size(AllCombos,2), size(AllCombos,1) + 2);
all_Qs = zeros(G.numedges, size(AllCombos,2));

for combo = 1:size(AllCombos,2)
    
    fprintf('Solving parameter combination %g....\n',combo)
    ParamSet = combo;
    CapRgj = CapRgjVals(AllCombos(1,combo));
    seed = Seeds(AllCombos(2,combo));
    PAGkirbarRatio = PAGkirbarRatioVals(AllCombos(3,combo));
    cECGkirbarRatio = cECGkirbarRatioVals(AllCombos(4,combo));
    StimulateArtNodesCond = StimulateArtNodesConds{(AllCombos(5,combo))};
    StimulateCapillariesCond = StimulateCapillariesConds{(AllCombos(6,combo))};
    
    conds = AllCombos(:,combo);
    
    % solve for each combination
    scriptFig6A
    
    tInstance = 4;  % [s] time instance for plotting on geometry
    tplot = t/1000;
    tInd = find(tplot>=tInstance,1,'first');
    MeanArtVm = mean(Vm_post(tInd,MainArtNodes));
    TopArtVm = Vm_post(tInd,TopArtNode);
    summaryData(combo,:) = [conds', MeanArtVm, TopArtVm];
    
    % for cubic spline
    nodeval = Vm_post(tInd,:);
    edgeval = getedgeval(G,nodeval);
    Q = edgeval;
    all_Qs(:,combo) = Q;
    eval(sprintf('G.Edges.Q%d = Q;',combo))
    
    close all
end

%% spline network plots
disp('Plotting....')

for kk = 1:size(all_Qs,2)
    
    Q = all_Qs(:,kk);
    valname = sprintf('Q%d',kk);
    if kk == 1
        out = val_colored(valname, G, iview);
    else
        
        out.H1.Edges.ToColor = out.H1.Edges.(valname);
        colorcanvas1(out.H1, out.all_hs1, out.allEdgeIDs1)
        fig.Visible = 'on';
        fig.Color = 'w';
        
        out.H2.Edges.ToColor = out.H2.Edges.(valname);
        colorcanvas1(out.H2, out.all_hs2, out.allEdgeIDs2)
        fig.Visible = 'on';
        fig.Color = 'w';
    end
    
    caxis([-50,-25])
    
    
    set(gcf, 'InvertHardcopy', 'off')
    
end




