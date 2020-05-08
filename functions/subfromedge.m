function [subG,nodeList] = subfromedge(G,edgeInd)
% DESCRIPTION:
% Make a subgraph from edge indicies

% INPUT: 
% G --- graph object
% edgeInd --- index of edges to keep 

edgelist = G.Edges{:,1};

edgeToKeep = edgelist(edgeInd,:);

span = 2*length(edgeToKeep);

nodeList = unique(reshape(edgeToKeep,[span 1]));

subG = subgraph(G,nodeList);

