function H = plotgraph(G)

iart = find(G.Edges.Type == 1);
iven = find(G.Edges.Type == 2);
icap = find(G.Edges.Type == 0);
ipart = find(G.Edges.Type == 3);
ipven = find(G.Edges.Type == 4);

set(gcf,'Color','w')
H = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'ZData',G.Nodes.Z,...
    'Marker','none');
H.ShowArrows = 'off';
highlight(H,'Edges',iart,'EdgeColor','r','LineWidth',3)
highlight(H,'Edges',iven,'EdgeColor','b','LineWidth',3)
highlight(H,'Edges',icap,'EdgeColor',[0.25, 0.25, 0.25],'LineWidth',1)
highlight(H,'Edges',ipart,'EdgeColor','r','LineWidth',5)
highlight(H,'Edges',ipven,'EdgeColor','b','LineWidth',5)
ax = gca;
ax.FontWeight = 'bold';
ax.FontName = 'arial';
ax.FontSize = 16;
ax.AmbientLightColor = 'magenta';
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.ZGrid = 'off';
ax.LineWidth = 2;
ax.XAxis.Label.String = 'x (\mum)';
ax.YAxis.Label.String = 'y (\mum)';
ax.ZAxis.Label.String = 'z (\mum)';

ax.ZLim = [-1200 0];

campos(1e3*[-5.8900   -2.3044    0.0125])
axis equal
box off
