function [g, all_hs, allEdgeIDs, fig] = prepfig(out)

g = out.Graph;
all_hs = out.Handles;
allEdgeIDs = out.EdgeIDs;
fig = out.Fig;

fig.Color = 'w';
zlabel('z(\mum)')
axis tight

ax = gca;
ax.Color = 'w';
ax.LineWidth = 3;
ax.FontSize = 22;
ax.FontName = 'arial';
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
ax.ZAxis.Visible = 'off';
