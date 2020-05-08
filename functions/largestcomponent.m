function [B] = largestcomponent(A)
% INPUT
% A ... adjecency matrix

% OUTPUT
% B ... indices of nodes belonging to the largest component

% n ... number of nodes
% mz .. set of neighbors
% x ... for each node: index of the connected component it belongs to
% v ... set of nodes of these current component, whose neighbors have to be checked still
% c ... for each component: number of its nodes

% deterimine number of nodes
n=length(A);

% loop over all nodes to get list of neighbors
for i=1:n
    mz{i}=find(A(i,:));
end

x(1:n)=0;
z=0; % number of components detected so far
k=0; % set of nodes already dealt with
for i=1:n
    % node i does not yet belong to any component?
    if x(i)==0;
        z=z+1; % increase number of detected components
        clear v
        v(1)=i;
        while nnz(v)>0 % while v not empty
            x(v(1))=z; % update component entry for current node v(1)
            k=union(k,v(1)); % flag current node v(1) as done
            b=setdiff(mz{v(1)},k); % neighbors of node v(1) that still have to be taken care of
            v=setdiff(v,v(1)); % delete current node from stack
            v=union(v,b); % add newly detected nodes to interim list of nodes of current component
        end
    end
end

% size of components
c(1:max(x))=0;
% loop over components
for i=1:max(x)
    c(i)=length(find(x==i)); % determine size of components
end

% size of the largest component(s)
cm=find(c==max(c));

% pick largest component
B=find(x==cm(1));
