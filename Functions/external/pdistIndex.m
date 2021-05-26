%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code by Selen Atasoy 31/10/2013
% returns the pdist index/indices for the given squareform indices
% function
%
% INPUT:    dist            pdist array
%           r               row index
%           c               column index
%
% OUTPUT:   idx             indcies in the pdist format
%           vals            values of these indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [idx, vals] = pdistIndex(dist, r_all, c_all)

% check the input
% if isequal(r,c) || (sum(ismember(r,c))>0) || (sum(ismember(c,r))>0)
%     display('diagonal entries are not present in the dist vector');
%     return;
% end;

ctr=1;
for r_ctr=1:length(r_all)
   for c_ctr=1:length(c_all)          
   
    % convert r and c to upper diagonal entries  
    if r_all(r_ctr)>c_all(c_ctr)
        temp    = r_all(r_ctr);
        r       = c_all(c_ctr);
        c       = temp;
    else
        r = r_all(r_ctr);
        c = c_all(c_ctr);
    end;

    % nr of data points
    n = ceil(sqrt(2*length(dist)));

    % number of previous elements
    % pn = n-1;
    % 
    % % number of previous rows
    % pr = (r-1);
    % sum2pr = (1 + pr) .* (pr./2);
    % 
    % % nr of previous row elements
    % pr_elements = (pr .* pn) - sum2pr;
    % idx = pr_elements + c - 1;

    idx(ctr) = (((r-1) * (n-1)) - ((1 + (r-1)) * ((r-1)/2))) + c - 1;
    vals(ctr)  = dist(idx(ctr));
    ctr = ctr +1;
   end;
end;
