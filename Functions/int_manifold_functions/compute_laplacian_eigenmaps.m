%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pearson Lab UNSW
% Code modified by Selen Atasoy
%
% original code by Michail Belkin, Laplacian Eigenmaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% returns the eigenvalues and the eigenvectors of the graphb Laplacian
%
% INPUT:
% A         -   Adjacency matrix
% d         -   number of dimensions
% type      -   type of the Laplacian (combinatorial or normalized)
%
% OUTPUT:
% E         -   eigenvalues
% V         -   eigenvectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Y, V] = computeLaplacianEigenmaps(A, d, type)

if ~exist('type', 'var')
    type = 'combinatorial';
end;

% disp(['Computing ', type, ' Laplacian eigenvectors.']);

W = A;
DegM = sum(W(:,:),2);

%% remove zero entries from the degree matrix
% [val, ind] = find(DegM==0)
% D(val(:,1))=1;
%%%

%% if asked, compute the normalized Laplacian
opts.tol    = 1e-9;
%opts.issym  = 1;
%opts.disp   = 0;% this has been modified, before it was = 5
switch type
    case 'normalized'
        L = spdiags(DegM,0,speye(size(W,1)))-W;
        D = full(spdiags(DegM,0,speye(size(W,1))));
        [E,V] = eig(full(L),D);
        E(:,1) = []; V(1,:)=[]; V(:,1)=[];
    case 'symmetric'
        %
        %     for i=1:size (L)
        %       D(i,i) = sum(L(i,:));
        %       if (D(i,i) ~= 0)
        %          DD(i,i) = 1/sqrt(DegM(i,i));
        %       else disp ('warning 0');
        %          DD(i,i) = 0;
        %       end
        %     end
        %     LL=DD*(DegM-L)*DD;
        %     L = LL;
        DD = spdiags(1./sqrt(DegM),0,speye(size(W,1)));
        L           = speye(size(W,1)) - ( DD * W * DD );
        [E,V] = eigs(L,d+1,'sm',opts);
        
    case 'combinatorial'
        L = spdiags(DegM,0,speye(size(W,1)))-W;
        % eigenvalue decomposition
        [E,V] = eigs(L,d+1,'sm',opts);
   
    case 'levys'
        L = spdiags(DegM,0,speye(size(W,1)))-W;
        L = 0.5*(L + transpose(L));
        [E,V] = eigs(L,d+1,'sm',opts);
end;

%% sort the einemaps and return as the new coordinates
[val, ind] = sort(diag(V));
Y = E(:,1:d);
V = val(1:d);


