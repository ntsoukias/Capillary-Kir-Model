clear, clc, close all

%% add necessary folders to the path
addpath(genpath(fullfile(cd ,'..\..')))

%% Constant parameters
R  = 8314.0	;	%[mJmol-1K-1]	gas constant
T = 293; %[K]			absolute temperature
F  = 96485.0;	%[Cmol-1]		Faraday's constant
RT_F = R*T/F;
z_K  = 1;    %		K ion valence

%% Kir channel characteristics
delta_V_kir  = 25; %[mV]	voltage diff at half-max. eI_K_i
G_kirbar = 0.18; %[nS/mM^0.5] inward rectifier constant
n_kir = 0.5;	%		inward rectifier constant
k_kir = 7;	% [mV]	inward rectifier slope factor
K_i = 150;   %[mM] intracullular potassium level

%% Background equlibrium potential
E_bg = -30;  %[mV] resting membrane potential

%% Ko and Gbg values
K_o_values = linspace(1,10,50);                    % Ko values for Kir current [mM]
Gbg_values = linspace(0,0.7,3000);              % Gbg values [nS]

%% Preallocate space for all_intersections, lower and upper Gbg and ratios, and critical Gbg and ratio
all_intersections = zeros(3,length(Gbg_values),length(K_o_values));         % All the roots for all Ko values,
cellState = zeros(length(K_o_values), length(Gbg_values));

lower_Gbg = zeros(1,numel(K_o_values));                 % lower bound for Gbg which result in bistability for each Ko value
lower_ratio = zeros(1,numel(K_o_values));                 % lower bound for Gbg/Gkirmax which result in bistability for each Ko value

upper_Gbg = zeros(1,numel(K_o_values));                % upper bound for Gbg which result in bistability for each Ko value
upper_ratio = zeros(1,numel(K_o_values));                % upper bound for Gbg/Gkirmax which result in bistability for each Ko value

Gbg_crit = zeros(1,numel(K_o_values));                     % Critical Gbg value for each Ko value (border between hyper fav., and depo. fav.)
crit_ratio = zeros(1,numel(K_o_values));                     % Critical Gbg/Gkirmax for each Ko values (border between hyper fav., and depo. fav.)

%% Find bistability region for each Ko value

lower_Vm = zeros(size(K_o_values));
upper_Vm = zeros(size(K_o_values));
crit_Vm = zeros(size(K_o_values));
dVmDep = zeros(size(K_o_values));
dVmHyp = zeros(size(K_o_values));
all_Iblue = zeros(numel(Gbg_values),numel(K_o_values));
all_Ired = zeros(numel(Gbg_values),numel(K_o_values));
all_dVmblue = zeros(numel(Gbg_values),numel(K_o_values));
all_dVmred = zeros(numel(Gbg_values),numel(K_o_values));

tic

f = waitbar(0, 'Finding the bistable region ...');
for j = 1:numel(K_o_values)
    
    K_o = K_o_values(j);
    E_K = RT_F/z_K*log(K_o/K_i);
    Gmax_kir = G_kirbar * (K_o^n_kir);                 % max kir channel conductance [nS]
    V_m = linspace(-200,0,500);                     % range of Vm for finding the roots [mV]
    
    for i = 1:numel(Gbg_values)
        
        % Find all the possible roots (intersections) for each Gbg value
        
        Gbg = Gbg_values(i);
        total_I = @(V_m) Gmax_kir *...
            (V_m - E_K)./ (1 + exp((V_m - E_K - delta_V_kir)/k_kir)) + Gbg*(V_m - E_bg);
        % find the intersections (roots) of total_I
        zV_m = V_m(total_I(V_m).*circshift(total_I(V_m),[0 -1]) <= 0);  % Estimate zero crossings
        zV_m = zV_m(1:end-1);                           % Eliminate any due to ‘wrap-around’ effect
        for k1 = 1:length(zV_m)
            roots(k1) = fzero(total_I, zV_m(k1));   % Finds the root(s) of
        end
        
        if exist('roots','var')
            if numel(roots) == 3
                roots = sort(roots);
                all_intersections(:,i,j) = roots;
                VmRegion = linspace(roots(1),roots(3),500);
                Iblue = abs(min(total_I(VmRegion)));
                Ired = abs(max(total_I(VmRegion)));
                dVmblue = roots(3) - roots(2);
                dVmred = roots(2) - roots(1);
                all_Iblue(i,j) = Iblue;
                all_Ired(i,j) = Ired;
                all_dVmblue(i,j) = dVmblue;
                all_dVmred(i,j) = dVmred;
                if Iblue < Ired
                    cellState(j,i) = -1;
                else
                    cellState(j,i) = 1;
                end
                clear roots
            elseif numel(roots) == 1
                all_intersections(1,i,j) = roots;
                if abs(roots - E_K) < abs(roots - E_bg)
                    cellState(j,i) = -2;
                else
                    cellState(j,i) = 2;
                end
                clear roots
            else
                clear roots
            end
        end
    end
    
    %% new
    % Find bistable range
    bistable_range = find(all(all_intersections(:,:,j)) ~= 0);
    
    if ~isempty(bistable_range)
        lowerInd = find(cellState(j,:) == -1,1,'first');
        lower_Gbg(j) = Gbg_values(min(bistable_range));
        lower_ratio(j) = lower_Gbg(j)/G_kirbar;
        
        upperInd = find(cellState(j,:) == 2,1,'first');
        upper_Gbg(j) = Gbg_values(max(bistable_range));
        upper_ratio(j) = upper_Gbg(j)/G_kirbar;
        
        critInd = find(cellState(j,:) == 1,1,'first');
        Gbg_crit(j) = Gbg_values(critInd);
        
        crit_ratio(j) = Gbg_crit(j)/G_kirbar;
        lower_Vm(j) = max(all_intersections(:,min(bistable_range),j));
        upper_Vm(j) = min(all_intersections(:,max(bistable_range),j));
        crit_Vm(j) = all_intersections(2,critInd,j);
        dVmDep(j) = crit_Vm(j) - all_intersections(1,critInd,j);
        dVmHyp(j) = crit_Vm(j) - all_intersections(3,critInd,j);
    else
        continue
    end
    waitbar(j/numel(K_o_values), f);
end
close(f)
indx = find(lower_Gbg ~= 0);

toc

%% plot regions with coloring based on driving force

% ------------------  make custom colormap
n = 10;                % number of colors

R = [1 0 0];   % color first row - red
G = [0 0 0];   % color middle row - white
B = [0 0 1];   % color last row - blue
cmap = [R; G; B];


[X,Y] = meshgrid(1:3,1:n);  % mesh of indices
cmap = interp2(X([1,round(n/2),n],:),Y([1,round(n/2),n],:),cmap,X,Y); % interpolate colormap
cmap1 = cmap(1:round(n/2),:);
cmap2 = cmap(round(n/2):end,:);

% --------------- depolarization
Fig = figure;
Fig.Name = 'Bistability';
Fig.Renderer = 'painters';
Fig.Units = 'inches';
Fig.OuterPosition = [1 1 7 6];
hAxesDep = axes;
hAxesDep.Color = 'none';
axis(hAxesDep,'off');

ylimit = [0,1.63];

colormap(hAxesDep, (cmap1));
x = [K_o_values(indx) fliplr(K_o_values(indx))];
y1 = [crit_ratio(indx) fliplr(upper_ratio(indx))];
c1 = [fliplr(dVmDep(indx)) zeros(size(indx))];


h1 = patch(x, y1,c1); h1.EdgeColor = 'k';
c = colorbar(hAxesDep,'Location','west');
cbar_position = c.Position;
c.Position = [cbar_position(1)+0.05, cbar_position(2)+0.3, 0.5*cbar_position(3), 0.5*cbar_position(4)];
c.Label.String = '\DeltaV_m(mV)';

c.LineWidth = 1;
c.FontSize = 15;
c.FontName = 'arial';
caxis(hAxesDep, [0,15]);

ax = hAxesDep;
xlabel('K_o(mM)'), ylabel('G_{bg}/G_{Kirbar}(mM)^{1/2}')
ax.XLim = [K_o_values(1), K_o_values(end)];
ax.YLim = ylimit;
ax.Position(1) = 0.15;
ax.Position(2) = 0.15;
positions = ax.Position;

y2 = [upper_ratio(indx) ax.YLim(2)*ones(size(indx))];
c2 = [ones(size(indx)) ones(size(indx))];
h1 = patch(x, y2,c2); h1.EdgeColor = 'none';
h1.FaceColor = 'y'; h1.FaceAlpha = 0.4;

% ------------------- hyperpolarization
hAxesHyp = axes;
hAxesHyp.Color = 'none';
colormap(hAxesHyp,(cmap2));
x = [K_o_values(indx) fliplr(K_o_values(indx))];
y2 = [crit_ratio(indx) fliplr(lower_ratio(indx))];

c2 = [dVmHyp(indx) zeros(size(indx)) ];

h2 = patch(x,y2,c2); h2.EdgeColor = 'k';
c = colorbar(hAxesHyp,'Location','east');
caxis(hAxesHyp, [-30, 0]);

cbar_position = c.Position;
c.Position = [cbar_position(1)-0.1, cbar_position(2)+0.1, 0.5*cbar_position(3), 0.5*cbar_position(4)];
c.Label.String = '\DeltaV_m(mV)';
c.LineWidth = 1;
c.FontSize = 15;
c.FontName = 'arial';
hAxesHyp.Position = positions;

ax = hAxesHyp;
ax.XLim = [K_o_values(1), K_o_values(end)];
ax.YLim = [ hAxesDep.YLim(1), hAxesDep.YLim(2)];
xlabel('K_o(mM)'), ylabel('G_{bg}/G_{Kirbar}(mM)^{1/2}')

ax.FontSize = 18;
ax.FontName = 'arial';
ax.LineWidth = 3;
positions = ax.Position;
hold on

x3 =  [K_o_values(indx)   K_o_values(indx(end)) K_o_values(indx(end)+1:end)...
    fliplr(K_o_values)];
y3 = [lower_ratio(indx) ...
    hAxesDep.YLim(2)*ones(size(K_o_values(indx(end):end))) ax.YLim(1)*ones(size(K_o_values))];
c3 = ones(size(y3));
h1 = patch(x3, y3,c3); h1.EdgeColor = 'none'; h1.FaceColor = 'k'; h1.FaceAlpha = 0.4;

% ---------------- link the two overlaying axes so they match at all times to remain accurate
linkaxes([hAxesDep,hAxesHyp]);


