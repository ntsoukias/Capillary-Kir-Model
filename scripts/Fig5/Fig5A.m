clear, clc, close all

%% add necessary folders to the path
addpath(genpath(fullfile(cd ,'..\..')))

%% Simulation time
Tmax = 6;                % [s] Total time of the simulation
tspan = [0, Tmax*1e3];      % converted to [ms]

%% geometry
s = load('Geometry_Fig5.mat');

g = s.G;
Adj = g.adjacency;
Adj = Adj + Adj';
Adj = full(Adj);

capEdges = find(g.Edges.Type == 0);
artEdges = find(g.Edges.Type == 2);
venEdges = find(g.Edges.Type == 1);

[~,capInd] = subfromedge(g,capEdges);
[~,artInd] = subfromedge(g,artEdges);
[~,venInd] = subfromedge(g,venEdges);

capInd(ismember(capInd, artInd)) = [];
capInd(ismember(capInd, venInd)) = [];

ArtECsInds = artInd;

%% combined geometry
tot_cell_number = length(Adj);

%% Plot full geometry
EdgeList = g.Edges{:,1};
G = g;

%% find average length of the edges connected to a node
avgL = zeros(1,tot_cell_number);
neighbors = zeros(tot_cell_number, max(G.indegree) + max(G.outdegree));

for i = 1:tot_cell_number
    outEdges = G.outedges(i);
    inEdges = G.inedges(i);
    allEdges = [outEdges; inEdges];
    NeighborNodes = [];
    for j = 1:numel(allEdges)
        NeighborNodes = [NeighborNodes, G.Edges.EndNodes(allEdges(j),:)];
    end
    NeighborNodes(NeighborNodes == i) = [];
    neighbors(i,1:numel(NeighborNodes)) = NeighborNodes;
    avgL(i) = sum(G.Edges.L(allEdges))/2;
end

circArtECs = 3; %number of Art ECs in the circumferential direction
circVenECs = 3;%number of Vein ECs in the circumferential direction
desiredLCap = 20;
desiredLArt = 33;
desiredL = desiredLArt*ones(1,tot_cell_number); %physiological length of cells (um)
desiredL(capInd) = desiredLCap;
LengthFactor = avgL./desiredL;
LengthFactor(artInd) = circArtECs*LengthFactor(artInd);
LengthFactor(venInd) = circVenECs*LengthFactor(venInd);

%% %%%%%%%%%%%%%%%% Stimulation protocol %%%%%%%%%%%%%%%%%%%%%%%%%
K_stim = true;                % extracellular potassium stimulation
PA_indices = ArtECsInds(:); % index of parenchymal arterioles

%% Potassium stimulation
K_o_Stim = 10;              % [mM] extracellular potassium concentration for stimulated cells
K_o_rest = 3*ones(1,tot_cell_number);                 % [mM] extracellular potassium concentration for unstimulated cells
K_stim_onset = 1;         % [s] K Stimulation onset
K_stim_end = Tmax + 1;     % [s] K Stimulation end
K_i = 150;                  % [mM] intracellular potassium concentration

%% %%%%%%%%%%%%%%%%%% model parameters %%%%%%%%%%%%%%%%%%%%%%%%

%% constant parameters
R  = 8314.0;	% [mJmol-1K-1]	gas constant
T = 293.0;      % [K] absolute temperature
F  = 96485.0;	% [Cmol-1] Faraday's constant
RT_F = R*T/F;   % RT/F
z_K  = 1;       % K ion valence

%% Kir channel characteristic
delta_V_kir  = 25;% [mV]	voltage diff at half-max. eI_K_i
n_kir = 0.5;	% inward rectifier constant
k_kir = 7;	% [mV]	inward rectifier slope factor

%% Gap junctional resistance
R_gjCap = 5e-2; %[Gohm]

R_gj = R_gjCap*Adj;  % [Gohm]  gap junctional resistance
edgenode = G.Edges{:,1};

for i = 1:tot_cell_number
    for j = 1:tot_cell_number
        if Adj(i,j) ~= 0
            edgeInd1 = find(ismember(edgenode,[i, j],'rows'));
            edgeInd2 = find(ismember(edgenode,[j, i],'rows'));
            if isempty(edgeInd1)
                edgeInd = edgeInd2;
            else
                edgeInd = edgeInd1;
            end
            ilength = G.Edges.L(edgeInd)/desiredLCap;
            itype = G.Edges.Type(edgeInd);
            
            if itype == 2
                fact = ilength/circArtECs;
            elseif itype == 1
                fact = ilength/circVenECs;
            else
                fact = ilength;
            end
            
            R_gj(i,j) = R_gj(i,j)*fact;
        end
    end
end

%% trials
NumTrials = 10;

%% Cell capacitance
Cm = 16*ones(1,tot_cell_number);   % [pF] membrane capacitance
Cm(capInd) = 8;
Cm = Cm.*LengthFactor;

%%
E_bg = -30;        % [mV] resting membrane potential
GKirbarMeanVals = [0.1, 0.4, 0.5, 0.6]/sqrt(3);% [nS/mM^0.5] inward rectifier constant


GKirbarStd = 0.2/sqrt(3);  %standard deviation of GKir [nS/mM0.5]
GbgMean = 0.06;    %[nS] background conductance
GbgStd = 0;  %[nS] standard deviation of Gbg

Percentage = [10, 12, 14];
stim_number_vals = round(Percentage*numel(capInd)/100);   %number of cells to be stimulated

stimPercentage = zeros(size(stim_number_vals));
AllMeanDeltaVmsArt = zeros(numel(GKirbarMeanVals),numel(stim_number_vals));
AllStdDeltaVmsArt = zeros(numel(GKirbarMeanVals),numel(stim_number_vals));
AllMeanDeltaVmsCap = zeros(numel(GKirbarMeanVals),numel(stim_number_vals));
AllStdDeltaVmsCap = zeros(numel(GKirbarMeanVals),numel(stim_number_vals));

AllVms = zeros(tot_cell_number, numel(GKirbarMeanVals), NumTrials, numel(stim_number_vals));
AllKStimInds = cell(1,numel(stim_number_vals));

%% store pars
ParNames = who();
for i = 1:numel(ParNames)
    eval(['pars.' ParNames{i} '=' ParNames{i} ';' ]);
end

%% solve

disp('Solving...')
for n = 1:numel(stim_number_vals)
    fprintf('Solving for percentage %g....\n',Percentage(n))
    
    number = stim_number_vals(n);    % number of cells to be stimulated
    AllDeltaVmsArt = zeros(numel(GKirbarMeanVals), NumTrials);
    AllDeltaVmsCap = zeros(numel(GKirbarMeanVals), NumTrials);
    IdealPercentage = Percentage(n);
    
    
    for i = 1:numel(GKirbarMeanVals)
        GKirbarMean = GKirbarMeanVals(i);
        
        for j = 1:NumTrials
            
            K_stim_cell_indx = randperm(numel(capInd), number);
            K_stim_cell_indx = capInd(K_stim_cell_indx);
            P = sum(LengthFactor(K_stim_cell_indx))/sum(LengthFactor(capInd(:)))*100;
            
            stimPercentage(n) = P;
            AllKStimInds{n} = K_stim_cell_indx;
            pars.K_stim_cell_indx = K_stim_cell_indx;
            
            G_kirbar = GKirbarMean + GKirbarStd*randn(1,tot_cell_number);
            G_bg = GbgMean + GbgStd*randn(1,tot_cell_number);
            G_kirbar(G_kirbar<=0) = 0;
            G_bg(G_bg<=0) = 0;
            
            % make PA's passive and not heterogeneous
            G_bg(PA_indices) = 0.28;
            G_kirbar(PA_indices) = 0.65/sqrt(3);
            G_kirbar(venInd) = 0;
            G_kirbar = G_kirbar.*LengthFactor;
            G_bg = G_bg.*LengthFactor;
            
            pars.G_kirbar = G_kirbar;
            pars.G_bg = G_bg;
            
            %% solve
            Vm0 = -30*ones(1,tot_cell_number);      %[mV] initial condition for membrane potentials
            [t,Vm] = ode15s(@(t,X)EqnsNew(t,X,pars),tspan,Vm0);
            
            %% find average change in arterial segment
            tplot = t/1000;     % convert to seconds
            t_indx1 = find(tplot<K_stim_onset,1,'last');
            t_indx2 = find(tplot<K_stim_end,1,'last');
            Vm1 = mean(Vm(t_indx1,PA_indices));
            Vm2 = mean(Vm(t_indx2,PA_indices));
            DeltaVmArt = Vm1 - Vm2;
            AllDeltaVmsArt(i,j) = DeltaVmArt;
            
            Vm1 = mean(Vm(t_indx1,capInd));
            Vm2 = mean(Vm(t_indx2,capInd));
            DeltaVmCap = Vm1 - Vm2;
            AllDeltaVmsCap(i,j) = DeltaVmCap;
            
            AllVms(:,i,j,n) = Vm(end,:);
            
        end
    end
    
    MeanDeltaVmsArt = -mean(AllDeltaVmsArt,2);
    StdDeltaVmsArt = std(AllDeltaVmsArt,0,2);
    
    AllMeanDeltaVmsArt(:,n) = MeanDeltaVmsArt;
    AllStdDeltaVmsArt(:,n) = StdDeltaVmsArt;
    
    MeanDeltaVmsCap = -mean(AllDeltaVmsCap,2);
    StdDeltaVmsCap = std(AllDeltaVmsCap,0,2);
    
    AllMeanDeltaVmsCap(:,n) = MeanDeltaVmsCap;
    AllStdDeltaVmsCap(:,n) = StdDeltaVmsCap;
    
end

%% grouped bar charts
DistInds = 1:numel(GKirbarMeanVals);

MeanDeltaVms = AllMeanDeltaVmsArt;
StdDeltaVms = AllStdDeltaVmsArt;

y = MeanDeltaVms;
e = StdDeltaVms;
fig = figure;
fig.Renderer = 'painters';
hB = bar(y); hold on

x=[];                      % collect the X values with offset for each
for i=1:length(hB)
    x=[x;hB(i).XData+hB(i).XOffset];
    hB(i).FaceColor = [1 1 1] - i/length(hB) ;
end
hold on
hEB = errorbar(x.',y,e,'.k','linewidth',2);  % add error bars

ylabel('\DeltaV_m(mV)')
ax = gca;
ax.LineWidth = 3;
ax.FontSize = 22;
ax.FontName = 'arial';
ax.Color = 'none';
ax.XAxisLocation = 'top';
ax.XTickLabel = '';
box off
ax.YLim = [-50,0];
ax.YTick = -50:10:0;

leg = legend([hB(1),hB(2),hB(3)],...
    {'10%','12%','14%'},'location','southwest');
leg.Box = 'off';

%% plot average on geometry
disp('Plotting on geometry...')
aveVms = mean(AllVms, 3);
counter = 1;
for p = 1:numel(Percentage)
    for d = 1:numel(DistInds)
        nodeval = aveVms(:,DistInds(d),1,p);
        edgeval = getedgeval(G,nodeval);
        Q = edgeval;
        climits = [-80,-20];
        KStimInd = AllKStimInds{p};
        CubicSplineNetworkNew(G,Q, KStimInd)
        
        set(gcf, 'InvertHardcopy', 'off')
        
    end
end


%% Equations
function dVmdt = EqnsNew(t,X,pars)

X = X(:);
Vm = X';

%% parameters
RT_F = pars.RT_F;
z_K= pars.z_K;
tot_cell_number= pars.tot_cell_number;
K_stim_cell_indx= pars.K_stim_cell_indx;
K_stim = pars.K_stim;
K_stim_onset = pars.K_stim_onset;
K_stim_end= pars.K_stim_end;
Cm = pars.Cm;
G_bg = pars.G_bg;
R_gj = pars.R_gj;
E_bg = pars.E_bg;
delta_V_kir = pars.delta_V_kir;
G_kirbar = pars.G_kirbar;
n_kir = pars.n_kir;
k_kir = pars.k_kir;
K_i = pars.K_i;
K_o_Stim  = pars.K_o_Stim;
K_o_rest = pars.K_o_rest;
neighbors = pars.neighbors;

%% --------------------------------------------------- Stimulation protocol
% Potassium stimulus
K_o = K_o_rest;

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
I_kir = G_kirbar.*(K_o).^n_kir .*((Vm - E_K)./(1 + exp((Vm - E_K - delta_V_kir)/k_kir)));   %[pA] whole cell kir current

%% Gap junctional current [pA]

I_gj = zeros(1,tot_cell_number);
for i = 1:tot_cell_number
    for j = 1:nnz(neighbors(i,:))
        I_gj(i) = I_gj(i) + 1/R_gj(i,neighbors(i,j))*(Vm(i) - Vm(neighbors(i,j)));
    end
end


%% Total membrane current
I_tot = I_bg + I_kir;

%% Differential equation
dVmdt =  -1./Cm .* (I_tot + I_gj - I_stim);
dVmdt = dVmdt';

end


function CubicSplineNetworkNew(G,Q, KStimInd)

A = G.adjacency;
D = normalize(G.Edges.D,'range',[1 8]);
D(G.Edges.Type == 1) = 9;
D(G.Edges.Type == 2) = 9;
D(G.Edges.Type == 0) = 3;
REdges = D/2;

npts = length(A);
nodenames = cell(1,npts);

for ii = 1:npts
    nodenames{ii} = num2str(ii);
end

AttributeList = ...
    {'EndNodes' 'Weight' 'REdges', 'Q'};
edgenod = G.Edges{:,1};
nseg = numedges(G);
EdgeTable = table([edgenod(:,1) edgenod(:,2)],ones(nseg,1),REdges, Q,...
    'VariableNames',AttributeList);

X = G.Nodes.X; Y = G.Nodes.Y; Z = G.Nodes.Z;
NodeTable = table(X,Y,Z,nodenames', 'VariableNames',{'X' 'Y' 'Z' 'nodenames'});
g = graph(EdgeTable,NodeTable);
xyz = [g.Nodes.X, g.Nodes.Y, g.Nodes.Z];


%% remove edges of the nodes with degree higher than 2
% find all nodes with degrees higher than 2
Bif_nodeIDs = find(g.degree>2);
g_modified = g;
for ii = 1:numel(Bif_nodeIDs)
    g_modified = g_modified.rmedge(g_modified.outedges(Bif_nodeIDs(ii)));
end

%%
[bins,binsizes] = g_modified.conncomp('outputform','cell');
strands = bins;

%%

remove_indx = [];
badstrands = [];

for ii = 1:numel(strands)
    strand = strands{ii};
    if isequal(numel(strand), 1)
        continue
    end
    
    End_nodes = strand((ismember(g_modified.degree(strand),[0,1])));
    full_strand = strand;
    if ~isempty(End_nodes)
        for jj = 1:numel(End_nodes)
            Add_Neighbor = intersect(g.neighbors(End_nodes(jj)), Bif_nodeIDs);
            full_strand = [Add_Neighbor', full_strand];
        end
    end
    full_strand = unique(full_strand,'stable');
    strands{ii} = full_strand;
    
    % check if the order is okay
    check = 1;
    istrand = strands{ii};
    for i = 1:numel(istrand) - 1
        nod1 = istrand(i);
        nod2 = istrand(i + 1);
        iedge = findedge(g, nod1, nod2);
        if isequal(iedge, 0)
            check = 0;
            break
        end
    end
    
    if check
        continue
    end
    
    sg = subgraph(g,full_strand);
    
    % find degree ones of the strand
    NodeDegrees = sg.degree(1:numel(full_strand));
    Deg1Nodes = find(NodeDegrees == 1);
    if isempty(Deg1Nodes)
        % loop is happening
        % remove one egde
        nodeIDs = sg.Edges{1,1};
        node1 = nodeIDs(1);
        node2 = nodeIDs(2);
        sg1 = sg;
        sg1 = rmedge(sg1,node1,node2);
        full_strand = full_strand(shortestpath(sg1,node1,node2));
        strands{ii} = full_strand;
        continue
    end
    startNode = Deg1Nodes(1);
    new_strand = [startNode];
    sg_modified = sg;
    
    while length(new_strand)<length(full_strand)
        % find neighbors of start node
        neighbor_nodes = sg_modified.neighbors(startNode);
        new_strand = [new_strand, neighbor_nodes];
        sg_modified = sg_modified.rmedge(sg_modified.outedges(startNode));
        startNode = new_strand(end);
    end
    new_strand = unique(new_strand,'stable');
    strands{ii} = full_strand(new_strand);
    
end

%%
ind = find(binsizes == 1);
count = numel(strands);
for i = 1:numel(ind)
    istrand = strands{ind(i)};
    if isstring(istrand)
        strand = str2double(istrand);
    else
        strand = istrand;
    end
    Add_Neighbor = g.neighbors(strand);
    for j = 1:numel(Add_Neighbor)
        count = count + 1;
        strands{count} = [strand, Add_Neighbor(j)];
    end
end

%% plot with cubic spline

disp('plotting.......')
tic
CubicSplineFig = figure;
set(gcf, 'visible','off');
hold all

n = 30;  % number of faces for tubes
inner_points = 3;    % number of points in the axial direction for each segment

errcount = 0;
errs = [];
count = 1;

for s = 1:length(strands)
    strand = strands{s};
    if numel(strand) == 1
        continue
    end
    
    XYZ = xyz(strand,:)';
    cs = cscvn(XYZ);
    for ii = 1:length(cs.breaks) - 1
        nod1 = strand(ii);
        nod2 = strand(ii + 1);
        iedge = findedge(g, nod1, nod2);
        if iedge == 0
            errcount = errcount + 1;
            errs = [errs, s];
            continue
        else
            rval = g.Edges.REdges(iedge);
            strandQ = g.Edges.Q(iedge);
        end
        xtest = linspace(cs.breaks(ii),cs.breaks(ii+1),inner_points);
        ytest = linspace(cs.breaks(ii),cs.breaks(ii+1),inner_points);
        ztest = linspace(cs.breaks(ii),cs.breaks(ii+1),inner_points);
        
        ctestx = 0;
        ctesty = 0;
        ctestz = 0;
        
        for jj = 1:cs.order
            ctestx = ctestx + cs.coefs((ii-1)*3+1,jj)*(xtest - cs.breaks(ii)).^(cs.order - jj);
            ctesty = ctesty + cs.coefs((ii-1)*3+2,jj)*(ytest - cs.breaks(ii)).^(cs.order - jj);
            ctestz = ctestz + cs.coefs((ii-1)*3+3,jj)*(ztest - cs.breaks(ii)).^(cs.order - jj);
        end
        
        [X,Y,Z] = tubeplot([ctestx;ctesty;ctestz],rval,n,0);
        h = surface(X,Y,Z,'edgecolor','none');
        all_hs(count) = h;
        allEdgeIDs(count) = iedge;
        count = count + 1;
        
        
    end
end

fig = figure(CubicSplineFig);
fig.Color = 'w';

shading interp
view(44,6)
brighten(1)


map = jet;
colormap(flipud(map))
caxis([-80 -20])

ax = gca;
ax.Color = 'k';
ax.LineWidth = 1;
ax.FontSize = 16;
ax.FontName = 'arial';
ax.FontWeight = 'bold';
ax.XTick = [];
ax.YTick = [];
ax.ZTick = [];
axis image;
axis equal
toc

%% assign color for edges
for i = 1:numel(all_hs)
    h = all_hs(i);
    iedge = allEdgeIDs(i);
    strandQ = g.Edges.Q(iedge);
    h.CData = strandQ*ones(size(h.ZData));
end

%% add spheres at stimulated nodes
stimNodes = KStimInd;
for i = 1:length(stimNodes)
    nodeID = stimNodes(i);
    rs = 3*mean(g.Edges.REdges);          % radius for the sphere representing the cells
    [Xs,Ys,Zs] = sphere(20);
    Xs = rs*Xs + g.Nodes.X(nodeID);
    Ys = rs*Ys + g.Nodes.Y(nodeID);
    Zs = rs*Zs + g.Nodes.Z(nodeID);
    
    h = surf(Xs,Ys,Zs,'linestyle','none');
    h.FaceColor = 'w';
end
end

function [x,y,z]=tubeplot(curve,r,n,ct)
% Usage: [x,y,z]=tubeplot(curve,r,n,ct)
%
% Tubeplot constructs a tube, or warped cylinder, along
% any 3D curve, much like the build in cylinder function.
% If no output are requested, the tube is plotted.
% Otherwise, you can plot by using surf(x,y,z);
%
% Example of use:
% t=linspace(0,2*pi,50);
% tubeplot([cos(t);sin(t);0.2*(t-pi).^2],0.1);
% daspect([1,1,1]); camlight;
%
% Arguments:
% curve: [3,N] vector of curve data
% r      the radius of the tube
% n      number of points to use on circumference. Defaults to 8
% ct     threshold for collapsing points. Defaults to r/2
%
% The algorithms fails if you have bends beyond 90 degrees.
% Janus H. Wesenberg, july 2004

if nargin<3 || isempty(n), n=8;
    if nargin<2, error('Give at least curve and radius');
    end;
end;
if size(curve,1)~=3
    error('Malformed curve: should be [3,N]');
end;
if nargin<4 || isempty(ct)
    ct=0.5*r;
end


%Collapse points within 0.5 r of each other
npoints=1;
for k=2:(size(curve,2)-1)
    if norm(curve(:,k)-curve(:,npoints))>ct;
        npoints=npoints+1;
        curve(:,npoints)=curve(:,k);
    end
end
%Always include endpoint
if norm(curve(:,end)-curve(:,npoints))>0
    npoints=npoints+1;
    curve(:,npoints)=curve(:,end);
end

%deltavecs: average for internal points.
%           first strecth for endpoitns.
dv=curve(:,[2:end,end])-curve(:,[1,1:end-1]);

%make nvec not parallel to dv(:,1)
nvec=zeros(3,1);
[buf,idx]=min(abs(dv(:,1))); nvec(idx)=1;

xyz=repmat([0],[3,n+1,npoints+2]);

%precalculate cos and sing factors:
cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);

%Main loop: propagate the normal (nvec) along the tube
for k=1:npoints
    convec=cross(nvec,dv(:,k));
    convec=convec./norm(convec);
    nvec=cross(dv(:,k),convec);
    nvec=nvec./norm(nvec);
    %update xyz:
    xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1])+...
        cfact.*repmat(r*nvec,[1,n+1])...
        +sfact.*repmat(r*convec,[1,n+1]);
end;

%finally, cap the ends:
xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);

%,extract results:
x=squeeze(xyz(1,:,:));
y=squeeze(xyz(2,:,:));
z=squeeze(xyz(3,:,:));

%... and plot:
if nargout<3, surf(x,y,z); end;
end



