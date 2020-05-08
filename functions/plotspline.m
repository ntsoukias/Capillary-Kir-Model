function out = plotspline(G,n , inner_points, fig)
% Graph object containing edge property to be colored
% edgetocolor must be edge property name as a string

dbstop if error

%% New graph containing edgeproperty to be colored (Q) and diameter (D).
nnod = numnodes(G);
nodenames = cell(1,nnod);

for i = 1:nnod
    nodenames{i} = num2str(i);
end

D = G.Edges.D;
if min(D) > 0.5
    D = normalize(G.Edges.D,'range',[1 8]);
    D(G.Edges.Type == 1) = 19;
    D(G.Edges.Type == 2) = 12;
    D(G.Edges.Type == 0) = 2.5;
    D(G.Edges.Type == 3) = 24;
    D(G.Edges.Type == 4) = 24;
end

G.Edges.D = D(:);
EdgeTable = G.Edges;
NodeTable = G.Nodes;
g = graph(EdgeTable,NodeTable);
g.Nodes.nodenames = nodenames';
xyz = [g.Nodes.X, g.Nodes.Y, g.Nodes.Z];

%% remove edges of the nodes with degree higher than 2
% find all nodes with degrees higher than 2
Bif_nodeIDs = find(g.degree > 2);
g_modified = g;
for i = 1:numel(Bif_nodeIDs)
    g_modified = g_modified.rmedge(g_modified.outedges(Bif_nodeIDs(i)));
end

%%
[bins,binsizes] = g_modified.conncomp('outputform','cell');
strands = bins;

%%
for i = 1:numel(strands)
    strand = strands{i};
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
    strands{i} = full_strand;
    
    % check if the order is okay
    check = 1;
    istrand = strands{i};
    for j = 1:numel(istrand) - 1
        nod1 = istrand(j);
        nod2 = istrand(j + 1);
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
        strands{i} = full_strand;
        continue
    end
    startNode = Deg1Nodes(1);
    new_strand = startNode;
    sg_modified = sg;
    
    while length(new_strand)<length(full_strand)
        % find neighbors of start node
        neighbor_nodes = sg_modified.neighbors(startNode);
        new_strand = [new_strand, neighbor_nodes];
        sg_modified = sg_modified.rmedge(sg_modified.outedges(startNode));
        startNode = new_strand(end);
    end
    new_strand = unique(new_strand,'stable');
    strands{i} = full_strand(new_strand);
    
end

ind = find(binsizes == 1);
count = numel(strands);
% f = waitbar(0, 'Handling degree one strands');
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
%     msg = ['Handling degree one strands...Iteration: ' sprintf('%i', i)...
%         ' of ' num2str(numel(ind))];
end

%% plot with cubic spline

% disp('plotting.......')
CubicSplineFig = fig;
set(gcf, 'visible','off');
hold all

errcount = 0;
errs = [];
count = 1;

radiusList = g.Edges.D/2;

for s = 1:length(strands)
    strand = strands{s};
    if numel(strand) == 1
        continue
    end
    
    XYZ = xyz(strand,:)';
    cs = cscvn(XYZ);
    for i = 1:length(cs.breaks) - 1
        nod1 = strand(i);
        nod2 = strand(i + 1);
        iedge = findedge(g, nod1, nod2);
        if iedge > 0
            rval = radiusList(iedge);
        else
            errcount = errcount + 1;
            errs = [errs, s];
            continue
        end
        xtest = linspace(cs.breaks(i),cs.breaks(i+1),inner_points);
        ytest = linspace(cs.breaks(i),cs.breaks(i+1),inner_points);
        ztest = linspace(cs.breaks(i),cs.breaks(i+1),inner_points);
        
        ctestx = 0;
        ctesty = 0;
        ctestz = 0;
        
        for jj = 1:cs.order
            ctestx = ctestx + cs.coefs((i-1)*3+1,jj)*(xtest - cs.breaks(i)).^(cs.order - jj);
            ctesty = ctesty + cs.coefs((i-1)*3+2,jj)*(ytest - cs.breaks(i)).^(cs.order - jj);
            ctestz = ctestz + cs.coefs((i-1)*3+3,jj)*(ztest - cs.breaks(i)).^(cs.order - jj);
        end
        
        if numel(rval) > 1
            rval = rval(1);
            iedge = iedge(1);
        end
        
        [X,Y,Z] = tubeplot([ctestx;ctesty;ctestz],rval,n,0);
        h = surface(X,Y,Z,'edgecolor','none');
        all_hs(count) = h;
        allEdgeIDs(count) = iedge;
        count = count + 1;
        
    end
end

% fprintf('Error Count during plotting: %i\n', errcount)

out.Graph = g;
out.Handles = all_hs;
out.EdgeIDs = allEdgeIDs;
out.Fig = CubicSplineFig;

end

function [x,y,z]=tubeplot(curve,r,n,ct)

dbstop if error

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
    end
end
if size(curve,1)~=3
    error('Malformed curve: should be [3,N]');
end
if nargin<4 || isempty(ct)
    ct=0.5*r;
end


%Collapse points within 0.5 r of each other
npoints=1;
for k=2:(size(curve,2)-1)
    if norm(curve(:,k)-curve(:,npoints))>ct
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
[~,idx]=min(abs(dv(:,1))); nvec(idx)=1;

xyz=zeros([3,n+1,npoints+2]);

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
end

%finally, cap the ends:
xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);

%,extract results:
x=squeeze(xyz(1,:,:));
y=squeeze(xyz(2,:,:));
z=squeeze(xyz(3,:,:));

%... and plot:
if nargout<3
    surf(x,y,z);
end

end
