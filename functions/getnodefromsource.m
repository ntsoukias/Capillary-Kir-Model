function nodeList = getnodefromsource(G, sourcenode, xdist, ydist, zdist)
% DESCRIPTION:
% Get all the nodes within a cube centered at the coordinates of
% a source node.

% INPUT:
% G --- graph object
% sourcenode --- node index that will be the center of a cube
% xdist, ydist, zdist --- distance in each direction to consider

X = G.Nodes.X;
Y = G.Nodes.Y;
Z = G.Nodes.Z;

sourceX = X(sourcenode);
sourceY = Y(sourcenode);
sourceZ = Z(sourcenode);

xrange = [sourceX - xdist, sourceX + xdist];
yrange = [sourceY - ydist, sourceY + ydist];
zrange = [sourceZ - zdist, sourceZ + zdist];

X(X <= xrange(1)) = 0;
X(X >= xrange(2)) = 0;

Y(Y <= yrange(1)) = 0;
Y(Y >= yrange(2)) = 0;

Z(Z <= zrange(1)) = 0;
Z(Z >= zrange(2)) = 0;

XYZ = [X, Y, Z];

trashnod = any(XYZ==0,2);
nodeList = 1:numnodes(G);
nodeList(trashnod) = [];
