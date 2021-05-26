function [data] = compute_TMH_ind_new(data,opts,verbose)
%compute_TMH_ind This function computes the temporal manifold harmonics of
%   the preprocessed fMRI data.
%   INPUT:
%           * data: cell(1,num_subjects)


% Compute Coherence Connectivity Dynamics matrix (CCD) for
% individual subjects
Nsubs = length(data);
if nargin==2
    verbose = true;
end

for nsub = 1:Nsubs
    if verbose
        fprintf('Computing TMH for subject %i of %i \n',nsub,Nsubs);    
    end    
   
% Compute adjacency matrix among time points
vecdist=pdist(data{1,nsub}.pattern,opts.distance);
Adjacency = rmst_from_vecdist(vecdist,opts.gamma,opts.maxKNN,opts.weighted);

% Compute Laplacian Eigenmaps

[TMH,TMHeig] = compute_laplacian_eigenmaps(Adjacency,opts.dimension+1,opts.laplacian);
data{1,nsub}.TMH = TMH(:,2:end);
data{1,nsub}.TMHeig = TMHeig(2:end);

end