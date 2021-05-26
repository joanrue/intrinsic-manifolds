function [A] = RMST_from_vecdist(vecdist,gamma,max_neighbors,weighted)
% Relaxed Minimum Spanning Tree Manifold (RMST) Learning Algorithm

if ~exist('edge_type','var') ||  weighted == 0
    weighted = 0;
end

% Create Minimum Spanning Tree
distmat = squareform(vecdist);
MST = minspantree(graph(distmat));
if weighted == 0
    A = adjacency(MST);
elseif weighted == 1
    A = adjacency(MST,'weighted');
end

% Apply the relaxed condition
for i = 1:size(A,1)
    [d_i,INDS] = sort(distmat(i,:),'ascend');
    d_i = d_i(2); % distance to the nearest neighbor of node "i"
    for j = INDS(2:max_neighbors+1)
        % Compute the shortest path in the MST among the nearest "max_neighbors"
        % neighbors
        weights = [];
        nodes = shortestpath(MST,i,j);
        for k = 1:(length(nodes)-1)
            weights = [weights,(distmat(nodes(k),nodes(k+1)))];
        end
        mw = max(weights);
        
        d_j = sort(distmat(j,:),'ascend');
        d_j = d_j(2);
        % Distance between i and j
        d_ij = distmat(i,j);
        
        % If distance between i and j is smaller than the max. weight plus
        % the minimum distances of the implicated nodes add an edge between
        % those nodes.
        
        if ((mw + gamma*(d_i+d_j)) > d_ij)
            if weighted == 0
                A(i,j)=1;
                A(j,i)=1;
            elseif weighted == 1
                A(i,j)=d_ij;
                A(j,i)=d_ij;
            end
        end
    end
end

% Transform distance to weighted adjacency
if weighted
    A(A>0) = max(max(A))./ A(A>0);
end
end
