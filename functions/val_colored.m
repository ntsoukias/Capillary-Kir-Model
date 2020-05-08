function out = val_colored(valname, G, iview)

%% Define colors based on vessel type
[G00, G12, G34] = splitgraph(G);

fig = figure;

n = 5; % number of faces for tubes
inner_points = 3; % number of points in the axial direction for each segment
out1 = plotspline(G00, n , inner_points, fig);

n = 30; % number of faces for tubes
inner_points = 200; % number of points in the axial direction for each segment
fig = out1.Fig;

out2 = plotspline(G12, n , inner_points, fig);
fig = out2.Fig;

if ~isequal(numedges(G34), 0)
out3 = plotspline(G34, n , inner_points, fig);
fig = out3.Fig;
end


%%
[g1, all_hs1, allEdgeIDs1, ~] = prepfig(out1);
[EdgeTableSorted1, ~] = sortrows(g1.Edges, {'CN'});
NodeTable1 = g1.Nodes;
H1 = graph(EdgeTableSorted1, NodeTable1);
H1.Edges.ToColor = H1.Edges.(valname);
colorcanvas1(H1, all_hs1, allEdgeIDs1)
fig.Visible = 'on';
fig.Color = 'w';

[g2, all_hs2, allEdgeIDs2, ~] = prepfig(out2);
[EdgeTableSorted2, ~] = sortrows(g2.Edges, {'CN'});
NodeTable2 = g2.Nodes;
H2 = graph(EdgeTableSorted2, NodeTable2);
H2.Edges.ToColor = H2.Edges.(valname);
colorcanvas1(H2, all_hs2, allEdgeIDs2)
fig.Visible = 'on';
fig.Color = 'w';

if ~isequal(numedges(G34), 0)
[g3, all_hs3, allEdgeIDs3, ~] = prepfig(out3);
[EdgeTableSorted3, ~] = sortrows(g3.Edges, {'CN'});
NodeTable3 = g3.Nodes;
H3 = graph(EdgeTableSorted3, NodeTable3);
H3.Edges.ToColor = H3.Edges.(valname);
colorcanvas1(H3, all_hs3, allEdgeIDs3)
fig.Visible = 'on';
fig.Color = 'w';
else
    H3 = [];
    all_hs3 = [];
    allEdgeIDs3 = [];
end

view(iview)
ax = gca;
ax.Color = 'k';
axis equal
cmap = jet;
colormap(flipud(cmap));
colorbar

camlight headlight

out.H1 = H1; 
out.H2 = H2;
out.H3 = H3;
out.all_hs1 = all_hs1;
out.all_hs2 = all_hs2;
out.all_hs3 = all_hs3;
out.allEdgeIDs1 = allEdgeIDs1;
out.allEdgeIDs2 = allEdgeIDs2;
out.allEdgeIDs3 = allEdgeIDs3;

