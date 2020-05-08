function h = plotcolor(G, variable, name, climits, diameter)

%% ------- Figure Properties ------- %%
fig = gcf;
fig.Color = 'w';

h = plot(G,'XData',G.Nodes.X,'YData',G.Nodes.Y,'ZData',G.Nodes.Z,...
    'EdgeCData',variable,'LineWidth',diameter,'EdgeAlpha',1,...
    'Marker','none','NodeLabel',[]);
h.ShowArrows = 'off';

%% ------- Axis Properties ------- %%
ax = gca;
ax.FontWeight = 'bold';
ax.FontName = 'arial';
ax.FontSize = 16;
ax.AmbientLightColor = 'magenta';
ax.XGrid = 'off';
ax.YGrid = 'off';
ax.ZGrid = 'off';
ax.XAxis.Label.String = 'x (\mum)';
ax.YAxis.Label.String = 'y (\mum)';
ax.ZAxis.Label.String = 'z (\mum)';
ax.Color = 'k';

%% ------- Colormap Properties ------- %%
c = colorbar;
caxis(climits)
c.Label.String = name;
c.Location = 'eastoutside';

c.FontSize = 12;
c.FontWeight = 'bold';
c.FontName = 'arial';

map = jet;
colormap(flipud(map)) % reverse the colormap
axis equal
