function colorcanvas3(g, all_hs, allEdgeIDs)

strandQList = g.Edges.ToColor;

%% assign color for edges
for i = 1:numel(all_hs)
    h = all_hs(i);
    iedge = allEdgeIDs(i);
    strandQ = strandQList(iedge, :);
    h.FaceColor = strandQ;
end


