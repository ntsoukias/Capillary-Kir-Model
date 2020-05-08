function [G0, G12, G34] = splitgraph(G)

% get all type indices

ind0 = find(G.Edges.Type == 0); % capillaries
ind1 = find(G.Edges.Type == 1); % arterioles
ind2 = find(G.Edges.Type == 2); % venules
ind3 = find(G.Edges.Type == 3); % arteries
ind4 = find(G.Edges.Type == 4); % veins

% Create graph of only Type 0
G0 = rmedge(G, [ind1(:); ind2(:); ind3(:); ind4(:)]);
% Create graph of only Type 1 and 2
G12 = rmedge(G, [ind0(:); ind3(:); ind4(:)]);
% Create graph of only Type 3 and 4
G34 = rmedge(G, [ind0(:); ind1(:); ind2(:)]);
