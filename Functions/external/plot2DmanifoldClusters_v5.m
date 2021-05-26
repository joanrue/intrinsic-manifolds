%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots 2D of the joint manifold
%
%
% INPUT:    Y           joint manifold
%           X           dataset a x b x n dimensional
%           nr_frames   number of frames of the first dataset
%           dims2plot   dimensions to plot
% 
% OUTPUT:   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot2DmanifoldClusters_v5(Y, dims2plot, num_images,sizeballs)

if (~exist('dims2plot', 'var'))
    dims2plot = [1,2];
end;

xy = [];
M = [];
count = 0; 
clusters = [];
for cmp=1:size(Y,2)
    xy_tmp = Y{cmp}.TMH(:,[dims2plot(1),dims2plot(2)]);
    xy = [xy; xy_tmp];
    for t = 1:size(Y{cmp}.pattern,1)
        count = count+1; 
        M(:,:,count) =  squareform(Y{cmp}.pattern(t,:));
        
    end
    clusters = [clusters ; abs(Y{cmp}.stages)+1];
end;

if (~exist('num_images', 'var'))
    num_images = 30;
end;

k = num_images;

%num_clusters = length(unique(clusters));

cc = colormap(brewermap(9,'set3'));
colorss = cc([2,2,6,7,5],:);
colorss = [0.5732         0    0.1330
    0.5732         0    0.1330
    colorss(3:5,:)];


%color_map = jet(num_clusters);
color_map = colorss;
C(:,1:3) = color_map(clusters(:),:);


%% plot the 1. video %%
if size(xy,2)>size(xy,1)
    xy = xy';
end


a = size(M,1);
b = size(M,2);
n = size(M,3);


e = 1/(2*k);
% plot result
for i=1:2
    xy(:,i) = rescale(xy(:,i), e, 1-e );
end

%A = zeros(k*b,k*a) + 1; % mmax(M(:));
A = zeros(k*b,k*a,3)+1; % mmax(M(:));
for x=1:k
    for y=1:k
        selx = (x-1)*b+1:x*b;
        sely = (y-1)*a+1:y*a;
        % center point
        cx = e + (x-1)/k;
        cy = e + (y-1)/k;
        % distance to center
        d = max( abs( xy(:,1)-cx ), abs( xy(:,2)-cy ) );
        % find closest point
        [v,I] = min(d);
        if v<=e
            %A(selx,sely) = rescale( M(end:-1:1,:,I)' );
            A(sely,selx,:) = 1-repmat(rescale( M(end:-1:1,:,I)')>.75,[1,1,3]).*permute(repmat(1-C(I,:)',[1,90,90]),[2,3,1]);
            
            %A(selx,sely,:) = repmat((rescale( -M(end:-1:1,:,I)')+.75)/1.75,[1,1,3]).*permute(repmat(C(I,:)',[1,90,90]),[2,3,1]);
        end
    end
end


n = size(A,1);
xy = xy*n;

plot_lines = 1;

hold on;
imagesc( [0,n-1]+0.5,[0,n-1]+0.5, A );
%colormap gray(256);
axis tight; axis square; axis off;

%% 1. video 
% scatter(xy(part1,1),xy(part1,2),10, C(part1,:), '*');
% for i=1:size(xy,1)
%     plot(xy(i,1),xy(i,2), 'o', ...
%         'MarkerSize',sizeballs, 'MarkerEdgeColor', C(i,:), 'MarkerFaceColor', C(i,:));
%     hold on;
% end;

%%
axis xy;
if plot_lines
    for x=0:a:n
        line([x x], [0 n], 'Color',[0.9 0.9 0.9]);
        line([0 n], [x x], 'Color',[0.9 0.9 0.9]);
    end
end
hold off;



