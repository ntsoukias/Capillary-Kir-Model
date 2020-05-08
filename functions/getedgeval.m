function edgeval = getedgeval(G,nodeval)
% DESCRIPTION: Loops through all edges in network. Assigns edge values
% based on average of ends nodes for each edge.

% INPUT:
% G       --- graph object
% nodeval --- value of each node

% OUTPUT:
% edgeval --- value of each edge

nseg = numedges(G);

edgenod = G.Edges{:,1}; % array of from-to nodes of each edge
edgeval = zeros(nseg,1);

for iseg = 1:nseg
    
    nod1 = edgenod(iseg,1);
    nod2 = edgenod(iseg,2);
    edgeval(iseg) = mean([nodeval(nod1) nodeval(nod2)]);
    
end
