%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%    TMH PIPELINE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%  Add external functions to path %%%%%%%%%%%%%%%%%%%%%%%%

addpath('external')
addpath('int_manifold_functions')
addpath('svm_functions')
addpath('surr_functions')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Load Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path2data = ['..' filesep 'Data'];

load([path2data filesep 'sleep_scoring_N0_subs'])
load([path2data filesep 'TC_N0_subs_REV'])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Define Parameters  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Relaxed Minimum Spanning Tree params.
opts.gamma = 3;
opts.maxKNN = 5;
opts.distance = 'cosine'; 
% Laplacian Eigenmaps params.
opts.weighted = 0;
opts.laplacian = 'combinatorial'; % Graph Laplacian can be combinatorial, normalized

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Subsampling  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select only those subjects that have a minjmum number of time points
% (TRs) per each condition, and select only a maximum number of time to
% have same number of samples per condition (for classification purposes)

opts.TSmax = 87;    % Max num of TRs per each condition.
opts.TSmin = 87;    % Min num of TRs per each condition (18 subjects with
% a minimum of 87 samples per condition. 12 subjects
% with a minimum of 144 samples per condition).

[ind_data,selected_subs,group_data] = preprocess(TC_all_states,inst_states_all,opts.TSmin,opts.TSmax);

Nsub = length(selected_subs);

%%
%%%%%%%%%%%%%%%%%%%  Get individual subject embeddings  %%%%%%%%%%%%%%%%%%%%%

opts.dimension = 90;

for nsub = 1:Nsub    
    [~,ind_data{1,nsub}.PCA,~]=pca(ind_data{1,nsub}.ts','NumComponents',opts.dimension);         
end

%%
%%%%%%%%%%%%%%%  SVM classification varying dimensionality  %%%%%%%%%%%%%%%%%

rng(1)
% load('../Data/SVM_RESULT_PROCUSTRES.mat');
svm_results_params = cell(opts.dimension-1,1);
parpool;
parfor i = 1:(opts.dimension-1)
    fprintf('SVM using %0.2i dimensions\n', i);
    
    %%%%%%%%%%%%%%%%%%%  Get all-subjects group embeddings  %%%%%%%%%%%%%%%%%%%%%
    
    opts_svm = opts;
    opts_svm.dimension = i;    
    ind_data_svm = compute_ind_manifold(ind_data,opts_svm,false);
    [group_data_svm,~] = compute_group_manifold(group_data,ind_data_svm,opts_svm,Nsub,87);    
    [group_data_svm,~] = compute_group_pca(group_data_svm,ind_data_svm,opts_svm,Nsub,87);            
    
    % 1vs1 SVM individual
    [acc_tmh,acc_pca,~] = svm_func_1vs1(ind_data_svm)
    svm_results_params{i}.acc_1vs1_tmh_ind = acc_tmh;
    svm_results_params{i}.acc_1vs1_pca_ind = acc_pca;
    
    % Multi-class SVM individual
    [acc_tmh,~,C_norm_tmh,acc_pca,~,C_norm_pca,~] = svm_func(ind_data_svm);
    svm_results_params{i}.acc_tmh_ind = acc_tmh;            
    svm_results_params{i}.C_norm_tmh_ind = C_norm_tmh;
    svm_results_params{i}.acc_pca_ind = acc_pca;            
    svm_results_params{i}.C_norm_pca_ind = C_norm_pca;        
    
    % 1vs1 SVM group
    [acc_tmh,acc_pca,conds_1vs1] = svm_func_1vs1(group_data_svm)
    svm_results_params{i}.acc_1vs1_tmh_group = acc_tmh;
    svm_results_params{i}.acc_1vs1_pca_group = acc_pca;
    
    % Multi-class SVM group
    [acc_tmh,~,C_norm_tmh,acc_pca,~,C_norm_pca,conds_multiclass] = svm_func(group_data_svm);
    svm_results_params{i}.acc_tmh_group = acc_tmh;            
    svm_results_params{i}.C_norm_tmh_group = C_norm_tmh;
    svm_results_params{i}.acc_pca_group = acc_pca;            
    svm_results_params{i}.C_norm_pca_group = C_norm_pca;    
end

[acc_ts,~] = svm_func_1vs1_ts(ind_data);
svm_results_hd.acc_1vs1_ind = acc_ts;            

[acc_ts,~,C_norm_ts] = svm_func_ts(ind_data);
svm_results_hd.acc_ind = acc_ts;            
svm_results_hd.C_norm_ind = C_norm_ts;

[acc_ts,~] = svm_func_1vs1_ts(group_data);
svm_results_hd.acc_1vs1_group = acc_ts;            

[acc_ts,~,C_norm_ts] = svm_func_ts(group_data);
svm_results_hd.acc_group = acc_ts;            
svm_results_hd.C_norm_group = C_norm_ts;

%% Test low-dimensionality

p_vals_ind = zeros(88,1); for i = 1:88, p_vals_ind(i) = ranksum( svm_results_params{i}.acc_tmh_ind, svm_results_params{i+1}.acc_tmh_ind ); end

[~, ~, ~, p_vals]=fdr_bh(p_vals);

%%
%%%%%%%%%%%%%%%  Surrogate analysis (Load data to save time) %%%%%%%%%%%%%%%%%
% 
 
% % Get p-values from surrogate analysis
% run surrogate_analysis in a cluster
% run organize_surrogate_data

load([path2data filesep 'Surrogates' filesep 'svm_results_surr'])

num_sur = 1354;

p_acc_1vs1_tmh_ind  = (1+sum(squeeze(mean(svm_results_surr.acc_1vs1_tmh_ind(1,:,:,:),[2,3]))>repmat(mean(svm_results_params{7}.acc_1vs1_tmh_ind,[1,2]),[num_sur,1])))/(1+num_sur);
p_acc_1vs1_pca_ind = (1+sum(squeeze(mean(svm_results_surr.acc_1vs1_pca_ind(1,:,:,:),[2,3]))>repmat(mean(svm_results_params{7}.acc_1vs1_pca_ind,[1,2]),[num_sur,1])))/(1+num_sur);
p_acc_1vs1_tmh_group = (1+sum(squeeze(mean(svm_results_surr.acc_1vs1_tmh_group(1,:,:,:),[2,3]))>repmat(mean(svm_results_params{7}.acc_1vs1_tmh_group,[1,2]),[num_sur,1])))/(1+num_sur);
p_acc_1vs1_pca_group = (1+sum(squeeze(mean(svm_results_surr.acc_1vs1_pca_group(1,:,:,:),[2,3]))>repmat(mean(svm_results_params{7}.acc_1vs1_pca_group,[1,2]),[num_sur,1])))/(1+num_sur);

p_acc_1vsall_tmh_ind = (1+sum(squeeze(sum(squeeze(mean(svm_results_surr.C_norm_tmh_ind(1,:,:,:,:),4)) .* repmat(eye(4),[1,1,num_sur]),[1,2])) > repmat(sum(mean(svm_results_params{7}.C_norm_tmh_ind,3) .* eye(4),[1,2]),[num_sur,1]))) / (1+num_sur);
p_acc_1vsall_pca_ind = (1+sum(squeeze(sum(squeeze(mean(svm_results_surr.C_norm_pca_ind(1,:,:,:,:),4)) .* repmat(eye(4),[1,1,num_sur]),[1,2])) > repmat(sum(mean(svm_results_params{7}.C_norm_pca_ind,3) .* eye(4),[1,2]),[num_sur,1]))) / (1+num_sur);
p_acc_1vsall_tmh_group = (1+sum(squeeze(sum(squeeze(mean(svm_results_surr.C_norm_tmh_group(1,:,:,:,:),4)) .* repmat(eye(4),[1,1,num_sur]),[1,2])) > repmat(sum(mean(svm_results_params{7}.C_norm_tmh_group,3) .* eye(4),[1,2]),[num_sur,1]))) / (1+num_sur);
p_acc_1vsall_pca_group = (1+sum(squeeze(sum(squeeze(mean(svm_results_surr.C_norm_pca_group(1,:,:,:,:),4)) .* repmat(eye(4),[1,1,num_sur]),[1,2])) > repmat(sum(mean(svm_results_params{7}.C_norm_pca_group,3) .* eye(4),[1,2]),[num_sur,1]))) / (1+num_sur);

p_acc_1vs1_ind_ranksum = ranksum(svm_results_params{7}.acc_1vs1_tmh_ind(:),svm_results_params{7}.acc_1vs1_pca_ind(:));
p_acc_1vs1_group_ranksum = ranksum(svm_results_params{7}.acc_1vs1_tmh_group(:),svm_results_params{7}.acc_1vs1_pca_group(:));
p_acc_1vsall_ind_ranksum = ranksum(svm_results_params{7}.C_norm_tmh_ind(repmat(eye(4),[1,1,18])==1),svm_results_params{7}.C_norm_pca_ind(repmat(eye(4),[1,1,18])==1));
p_acc_1vsall_group_ranksum = ranksum(svm_results_params{7}.C_norm_tmh_group(repmat(eye(4),[1,1,18])==1),svm_results_params{7}.C_norm_pca_group(repmat(eye(4),[1,1,18])==1));

[~,~,~, p] = fdr_bh([p_acc_1vs1_tmh_ind,p_acc_1vs1_pca_ind,p_acc_1vs1_tmh_group,p_acc_1vs1_pca_group,p_acc_1vsall_tmh_ind,p_acc_1vsall_pca_ind,p_acc_1vsall_tmh_group,p_acc_1vsall_pca_group]);
[~,~,~, p2] = fdr_bh([p_acc_1vs1_ind_ranksum,p_acc_1vs1_group_ranksum,p_acc_1vsall_ind_ranksum,p_acc_1vsall_group_ranksum]);

fprintf('\n******************************************\n')
fprintf('\nResults SVM -- Table S1\n')


fprintf('\nInd. TMH Average accuracy SVM 1vs1 = %0.2f +- %0.2f, pvalue against surrogates = %0.4f \n', mean(svm_results_params{7}.acc_1vs1_tmh_ind,[1,2]),std(svm_results_params{7}.acc_1vs1_tmh_ind,0,[1,2]), p(1))
fprintf('\nInd. PCA Average accuracy SVM 1vs1 = %0.2f +- %0.2f, pvalue against surrogates = %0.4f \n', mean(svm_results_params{7}.acc_1vs1_pca_ind,[1,2]),std(svm_results_params{7}.acc_1vs1_pca_ind,0,[1,2]), p(2))    
fprintf('\nInd. SVM 1vs1 ranksum = %0.4f \n', p2(1))
fprintf('\nGroup TMH Average accuracy SVM 1vs1 = %0.2f +- %0.2f, pvalue against surrogates = %0.4f \n', mean(svm_results_params{7}.acc_1vs1_tmh_group,[1,2]),std(svm_results_params{7}.acc_1vs1_tmh_group,0,[1,2]), p(3))
fprintf('\nGroup PCA Average accuracy SVM 1vs1 = %0.2f +- %0.2f, pvalue against surrogates = %0.4f \n', mean(svm_results_params{7}.acc_1vs1_pca_group,[1,2]),std(svm_results_params{7}.acc_1vs1_pca_group,0,[1,2]), p(4))
fprintf('\nGroup. SVM 1vs1 ranksum = %0.4f \n', p2(2))


fprintf('\nInd. TMH Average accuracy SVM 1vsAll = %0.2f +- %0.2f, pvalue against surrogates = %0.4f \n', mean(diag(mean(svm_results_params{7}.C_norm_tmh_ind,3) .* eye(4))),std(diag(mean(svm_results_params{7}.C_norm_tmh_ind,3) .* eye(4))), p(5))
fprintf('\nInd. PCA Average accuracy SVM 1vsAll = %0.2f +- %0.2f, pvalue against surrogates = %0.4f \n', mean(diag(mean(svm_results_params{7}.C_norm_pca_ind,3) .* eye(4))),std(diag(mean(svm_results_params{7}.C_norm_pca_ind,3) .* eye(4))), p(6))    
fprintf('\nInd. SVM 1vsall ranksum = %0.4f \n', p2(3))
fprintf('\nGroup TMH Average accuracy SVM 1vsAll = %0.2f +- %0.2f, pvalue against surrogates = %0.4f \n', mean(diag(mean(svm_results_params{7}.C_norm_tmh_group,3) .* eye(4))),std(diag(mean(svm_results_params{7}.C_norm_tmh_group,3) .* eye(4))), p(7))
fprintf('\nGroup PCA Average accuracy SVM 1vsAll = %0.2f +- %0.2f, pvalue against surrogates = %0.4f \n', mean(diag(mean(svm_results_params{7}.C_norm_pca_group,3) .* eye(4))),std(diag(mean(svm_results_params{7}.C_norm_pca_group,3) .* eye(4))), p(8))
fprintf('\nGroup. SVM 1vsall ranksum = %0.4f \n', p2(4))

fprintf('\n******************************************\n')
%%  
% 
% % Ranksum pvalues
comparisons = {'Awake-N1', 'Awake-N2', 'Awake-N3', 'N1-N2','N1-N3','N2-N3'};
states = {'Awake', 'N1', 'N2', 'N3'};

ord_1v1 = [6,5,3,4,2,1];
ord_1vsall = [4,3,2,1];

p_rksum_1vs1 = zeros(2,length(comparisons));
p_rksum_1vs_all = zeros(2,length(states));
 
for i = 1:6
    p_rksum_1vs1(1,i) = ranksum(svm_results_params{7}.acc_1vs1_tmh_ind(:,i),svm_results_params{7}.acc_1vs1_pca_ind(:,i));
    p_rksum_1vs1(2,i) = ranksum(svm_results_params{7}.acc_1vs1_tmh_group(:,i),svm_results_params{7}.acc_1vs1_pca_group(:,i));    
end

for i = 1:4
    p_rksum_1vs_all(1,i) = ranksum(squeeze(svm_results_params{7}.C_norm_tmh_ind(i,i,:)),squeeze(svm_results_params{7}.C_norm_pca_ind(i,i,:)));
    p_rksum_1vs_all(2,i) = ranksum(squeeze(svm_results_params{7}.C_norm_tmh_group(i,i,:)),squeeze(svm_results_params{7}.C_norm_pca_group(i,i,:)));
end

[~,~,~, p_fdr_1vs1] = fdr_bh(p_rksum_1vs1(:));
[~,~,~, p_fdr_1vsall] = fdr_bh(p_rksum_1vs_all(:));
p_rksum_1vs1 = reshape(p_fdr_1vs1,size(p_rksum_1vs1));
p_rksum_1vs_all = reshape(p_fdr_1vsall,size(p_rksum_1vs_all));
fprintf('\n******************************************\n')
fprintf('\nResults SVM 1vs1 -- Table S2\n')

for i = 1:length(comparisons)    
    fprintf('\n\nInd TMH  comparison %s,  Mean acc. %0.4f +- %0.4f  \n', comparisons{i}, mean(svm_results_params{7}.acc_1vs1_tmh_ind(:,ord_1v1(i))),std(svm_results_params{7}.acc_1vs1_tmh_ind(:,ord_1v1(i))))
    fprintf('Ind PCA  comparison %s,  Mean acc. %0.4f +- %0.4f  \n', comparisons{i}, mean(svm_results_params{7}.acc_1vs1_pca_ind(:,ord_1v1(i))),std(svm_results_params{7}.acc_1vs1_pca_ind(:,ord_1v1(i))))
    fprintf('Ranksum median comparison Ind. TMH vs PCA, %s,  p-value = %0.4f \n',comparisons{i}, p_rksum_1vs1(1,ord_1v1(i)))    
 
    fprintf('\nGroup TMH  comparison %s,  Mean acc. %0.4f +- %0.4f  \n', comparisons{i}, mean(svm_results_params{7}.acc_1vs1_tmh_group(:,ord_1v1(i))),std(svm_results_params{7}.acc_1vs1_tmh_group(:,ord_1v1(i))))
    fprintf('Group PCA  comparison %s,  Mean acc. %0.4f +- %0.4f  \n', comparisons{i}, mean(svm_results_params{7}.acc_1vs1_pca_group(:,ord_1v1(i))),std(svm_results_params{7}.acc_1vs1_pca_group(:,ord_1v1(i))))     
    fprintf('Ranksum median comparison Group TMH vs PCA, %s, p-value = %0.4f \n',comparisons{i}, p_rksum_1vs1(2,ord_1v1(i)))
end
fprintf('\n******************************************\n')


fprintf('\nResults SVM 1vsAll -- Table S2\n')
means_tmh_ind = diag(mean(svm_results_params{7}.C_norm_tmh_ind,3));
stds_tmh_ind = diag(std(svm_results_params{7}.C_norm_tmh_ind,0,3));
means_pca_ind = diag(mean(svm_results_params{7}.C_norm_pca_ind,3));
stds_pca_ind = diag(std(svm_results_params{7}.C_norm_pca_ind,0,3));
means_tmh_group = diag(mean(svm_results_params{7}.C_norm_tmh_group,3));
stds_tmh_group = diag(std(svm_results_params{7}.C_norm_tmh_group,0,3));
means_pca_group = diag(mean(svm_results_params{7}.C_norm_pca_group,3));
stds_pca_group = diag(std(svm_results_params{7}.C_norm_pca_group,0,3));

for i = 1:length(states)
    
    fprintf('\n\nInd TMH  comparison %s,  Mean acc. %0.4f +- %0.4f  \n', states{i}, means_tmh_ind(ord_1vsall(i)),stds_tmh_ind(ord_1vsall(i)))
    fprintf('Ind PCA  comparison %s,  Mean acc. %0.4f +- %0.4f  \n', states{i}, means_pca_ind(ord_1vsall(i)),stds_pca_ind(ord_1vsall(i)))
    fprintf('Ranksum median comparison Ind. TMH vs PCA, %s, p-value = %0.4f \n',states{i}, p_rksum_1vs_all(1,ord_1vsall(i)))    
    
    fprintf('\nGroup TMH  comparison %s,  Mean acc. %0.4f +- %0.4f  \n', states{i}, means_tmh_group(ord_1vsall(i)),stds_tmh_group(ord_1vsall(i)))
    fprintf('Group PCA  comparison %s,  Mean acc. %0.4f +- %0.4f  \n', states{i}, means_pca_group(ord_1vsall(i)),stds_pca_group(ord_1vsall(i)))
    fprintf('Ranksum median comparison Group TMH vs PCA, %s, p-value = %0.4f \n',states{i}, p_rksum_1vs_all(2,ord_1vsall(i)))
end
fprintf('\n******************************************\n')
%% Check avg acc for group of stages 

tmp = []; 
tmp2 = []; 

tmp3 = []; 
tmp4 = []; 

for i = 1:18, 
    tmp = [tmp,diag(svm_results_params{7}.C_norm_tmh_ind(ord_1vsall(3:4),ord_1vsall(3:4),i))];
    tmp2 = [tmp2,diag(svm_results_params{7}.C_norm_tmh_group(ord_1vsall(3:4),ord_1vsall(3:4),i))];
    
    tmp3 = [tmp3,diag(svm_results_params{7}.C_norm_tmh_ind(ord_1vsall(1:2),ord_1vsall(1:2),i))];
    tmp4 = [tmp4,diag(svm_results_params{7}.C_norm_tmh_group(ord_1vsall(1:2),ord_1vsall(1:2),i))];
end
fprintf('\n Ind TMH comparison stages 3&4,  Mean acc. %0.4f +- %0.4f  \n', mean(tmp(:)),std(tmp(:)))
fprintf('\n Group TMH comparison stages 3&4,  Mean acc. %0.4f +- %0.4f  \n', mean(tmp2(:)),std(tmp2(:)))


fprintf('\n Ind TMH comparison stages A&1,  Mean acc. %0.4f +- %0.4f  \n', mean(tmp3(:)),std(tmp3(:)))
fprintf('\n Group TMH comparison stages A&1,  Mean acc. %0.4f +- %0.4f  \n', mean(tmp4(:)),std(tmp4(:)))
%%
%%%%%%%%%%%%%%%%%%%%%%%  Area under the ROC Curve  %%%%%%%%%%%%%%%%%%%%%%%%%%
opts.dimension = 90;    
ind_data = compute_ind_manifold(ind_data,opts,false);
opts.dimension = 7;    
[group_data,~] = compute_group_manifold(group_data,ind_data,opts,Nsub,87);    
[group_data,~] = compute_group_pca(group_data,ind_data,opts,Nsub,87);    

[TMH_all, PCA_all] = roc_auc_manifolds(ind_data,opts,'median');

fprintf('\n******************************************\n')
fprintf('\nResults ROC analsysis -- Table S3\n')
for i = 1:length(comparisons)
    fprintf('\nTMH AUC-ROC %s = %0.2f +- %0.2f \n',comparisons{i}, mean(eval(sprintf('TMH_all.auc%0.1i',ord_1v1(i)))),std(eval(sprintf('TMH_all.auc%0.1i',ord_1v1(i)))))
    fprintf('PCA AUC-ROC %s = %0.2f +- %0.2f \n',comparisons{i}, mean(eval(sprintf('PCA_all.auc%0.1i',ord_1v1(i)))),std(eval(sprintf('PCA_all.auc%0.1i',ord_1v1(i)))))
    fprintf('Ranksum median comparison TMH vs PCA, p-value = %0.4f \n',TMH_all.pvals(i))
end
dimensions_used = [];
for i = 1:Nsub
    dimensions_used = [dimensions_used,length(unique(TMH_all.bestdim(i,:)))];
end
fprintf('\n Dimensions used --> Max = %0.01i,   Min = %0.01i,   Median = %0.01f \n',max(dimensions_used),min(dimensions_used),median(dimensions_used));
fprintf('\n******************************************\n')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%    PLOT RESULTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

size_ = 80;
alpha_ = 0.5;
cc = colormap(brewermap(9,'set3'));
colorss = cc([2,2,6,7,5],:);
colorss = [0.5732         0    0.1330
    0.5732         0    0.1330
    colorss(3:5,:)];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               FIGURE 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure of 2-D Individual TMH for figure 1

h1 = figure('Units','inches','Position',[0 0 10 10],'PaperPositionMode','auto'); hold on;

% Get all time points (for plotting manifold)

TSmax_tmp = 155; % Max num of TRs per each condition.
TSmin_tmp = 144; % Min num of TRs per each condition.
subj2plot = 52;

[Y,~,~] = preprocess(TC_all_states(subj2plot),inst_states_all(subj2plot),TSmin_tmp,TSmax_tmp);
Y = compute_ind_manifold(Y,opts);

dims2plot = [2,1];
num_images = 15;
size_balls = 8;
plot2DmanifoldClusters_v5(Y, dims2plot, num_images,size_balls)
save2pdf('../Figures/Fig1/CCD_graph',h1,300)

h2 = figure(2); clf
imagesc(squareform(Y{1}.pattern(10,:)),'AlphaData',0.4)
colormap gray(256);
xticks([]),yticks([])
daspect([1,1,1])
ax = gca;
ax.LineWidth = 20;
ax.XColor = colorss(2,:);
ax.YColor = colorss(2,:);
save2pdf('../Figures/Fig1/PhaseCoh1',h2,300)

h3 = figure(3); clf
imagesc(squareform(Y{1}.pattern(20,:)),'AlphaData',0.4)
colormap gray(256);
xticks([]),yticks([])
daspect([1,1,1])
ax = gca;
ax.LineWidth = 20;
ax.XColor = colorss(2,:);
ax.YColor = colorss(2,:);
save2pdf('../Figures/Fig1/PhaseCoh2',h3,300)

h4 = figure(4); clf
imagesc(squareform(Y{1}.pattern(30,:)),'AlphaData',0.4)
colormap gray(256);
xticks([]),yticks([])
daspect([1,1,1])
ax = gca;
ax.LineWidth = 20;
ax.XColor = colorss(2,:);
ax.YColor = colorss(2,:);
save2pdf('../Figures/Fig1/PhaseCoh3',h4,300)

h5 = figure(5); clf
imagesc(squareform(Y{1}.pattern(40,:)),'AlphaData',0.4)
colormap gray(256);
xticks([]),yticks([])
daspect([1,1,1])
ax = gca;
ax.LineWidth = 20;
ax.XColor = colorss(5,:);
ax.YColor = colorss(5,:);
save2pdf('../Figures/Fig1/PhaseCoh4',h5,300)

h6 = figure('Units','inches','Position',[0 0 4.5 4.5],'PaperPositionMode','auto'); hold on;
scatter3(Y{1}.TMH(:,1),Y{1}.TMH(:,2),Y{1}.TMH(:,3),...
    size_,colorss(abs(Y{1}.stages)+1,:),'o','filled','MarkerFaceAlpha',0.3)
view(3)

grid on
ax = gca;
box on
ax.LineWidth = 1.5;
save2pdf('../Figures/Fig1/Intrinsic_manifold',h6,300)

h7 = figure('Units','inches','Position',[0 0 4.5 4.5],'PaperPositionMode','auto'); hold on;

vecdist=pdist(Y{1}.pattern(:,[5,2,4]),'euclidean');
Adjacency = RMST_from_vecdist(vecdist,opts.gamma,opts.maxKNN,opts.weighted);
edges = find(Adjacency);
for e = 1:length(edges)
    [i,j] = ind2sub(size(Adjacency),edges(e));
    p = plot3(Y{1}.pattern([i,j],5),Y{1}.pattern([i,j],2),Y{1}.pattern([i,j],4),'-k','LineWidth',1);hold on
    p.Color(4) = 0.1;
end

scatter3(Y{1}.pattern(:,5),Y{1}.pattern(:,2),Y{1}.pattern(:,4), ...
    size_,colorss(abs(Y{1}.stages)+1,:),'o','filled','MarkerFaceAlpha',0.7)

view(3)
xlim([-1.2,1.2])
ylim([-1.2,1.2])
zlim([-1.2,1.2])
grid on
ax = gca;
box on
ax.LineWidth = 1.5;
save2pdf('../Figures/Fig1/graph_brain_dynamics',h7,300)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               FIGURE 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure of Individual low-dim intrinsic manifolds and PCA, and the corresponding 
% group embeedings
subj2plot = [3,4,6,10,11,14,16,18];

opt.FontSize = 20;

for i = subj2plot
    h8 = figure('Units','inches','Position',[0 0 4 4],'PaperPositionMode','auto'); hold on;
    
    ids_sub = (1+(87*4)*(i-1)):(87*4)*i;
    scatter3(group_data.TMH(ids_sub,2),group_data.TMH(ids_sub,1),group_data.TMH(ids_sub,3),...
        size_,colorss(abs(group_data.stages(ids_sub))+1,:),'o','filled','MarkerFaceAlpha',alpha_)
 
    xticks([-0.1, 0, 0.1]),yticks([-0.1, 0, 0.1]), zticks([-0.1, 0, 0.1])
    xticklabels([])
    yticklabels([])
    zticklabels([])
    
    grid on
    ax = gca;
    box on    
    ax.LineWidth = 1.5;
    h8 = fancy_figure(h8, opt);
    view(3)   
    save2pdf(sprintf('../Figures/Fig2/Embeddings_TMH_%i',i),h8,300)    
    close
    
    % PCA
    h8 = figure('Units','inches','Position',[0 0 4 4],'PaperPositionMode','auto'); hold on;
    ids_sub = (1+(87*4)*(i-1)):(87*4)*i;
    scatter3(group_data.PCA(ids_sub,2),group_data.PCA(ids_sub,1),group_data.PCA(ids_sub,3),...
        size_*2,colorss(abs(group_data.stages(ids_sub))+1,:),'o','filled','MarkerFaceAlpha',alpha_)
 
    xticks([-0.1, 0, 0.1]),yticks([-0.1, 0, 0.1]), zticks([-0.1, 0, 0.1])
    xticklabels([])
    yticklabels([])
    zticklabels([])
    
    grid on
    ax = gca;
    box on    
    ax.LineWidth = 1.5;
    h8 = fancy_figure(h8, opt);
    view(3)   
    save2pdf(sprintf('../Figures/Fig2/Embeddings_PCA_%i',i),h8,300)
    
    close
    
end

h8 = figure('Units','inches','Position',[0 0 6 6],'PaperPositionMode','auto'); hold on;

scatter3(group_data.TMH(:,2),group_data.TMH(:,1),group_data.TMH(:,3),...
    size_,colorss(abs(group_data.stages)+1,:),'o','filled','MarkerFaceAlpha',alpha_)

xlim([-0.10,0.10]),ylim([-0.1,0.1]),zlim([-0.1,0.1])
xticks([-0.10, 0, 0.1]),yticks([-0.1, 0, 0.1]), zticks([-0.1, 0, 0.1])
grid on,ax = gca;box on,ax.LineWidth = 1.5;
h8 = fancy_figure(h8, opt);
xticklabels([])
yticklabels([])
zticklabels([])
view(3)
save2pdf(sprintf('../Figures/Fig2/Embeddings_TMH_all',i),h8,300)
close

h8 = figure('Units','inches','Position',[0 0 2 2],'PaperPositionMode','auto'); hold on;

scatter3(group_data.PCA(:,1),group_data.PCA(:,2),group_data.PCA(:,3),...
    size_,colorss(abs(group_data.stages)+1,:),'o','filled','MarkerFaceAlpha',alpha_)

xticks([-20, 0, 20]),yticks([-10, 0, 10]), zticks([-10, 0, 10])
grid on,ax = gca;box on,ax.LineWidth = 1.5;
h8 = fancy_figure(h8, opt);
view(3)
xticklabels([])
yticklabels([])
zticklabels([])
save2pdf(sprintf('../Figures/Fig2/Embeddings_PCA_all',i),h8,300)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               FIGURE 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure of SVM accuracies with different number of dimensions

%% remove after final run
opts.dimension = 90;
%%

for i = 1:(opts.dimension-1)
    accuracies_TMH_ind(:,i) = svm_results_params{i}.acc_1vs1_tmh_ind(:);
    accuracies_PCA_ind(:,i) = svm_results_params{i}.acc_1vs1_pca_ind(:);
    accuracies_TMH_group(:,i) = svm_results_params{i}.acc_1vs1_tmh_group(:);
    accuracies_PCA_group(:,i) = svm_results_params{i}.acc_1vs1_pca_group(:);
    
    accuracies_TMH_ind_ci(1,i,:) = abs(prctile(svm_results_params{i}.acc_1vs1_tmh_ind(:),[75,25])-median(svm_results_params{i}.acc_1vs1_tmh_ind(:)));
    accuracies_PCA_ind_ci(1,i,:) = abs(prctile(svm_results_params{i}.acc_1vs1_pca_ind(:),[75,25])-median(svm_results_params{i}.acc_1vs1_pca_ind(:)));
    accuracies_TMH_group_ci(1,i,:) = abs(prctile(svm_results_params{i}.acc_1vs1_tmh_group(:),[75,25])-median(svm_results_params{i}.acc_1vs1_tmh_group(:)));
    accuracies_PCA_group_ci(1,i,:) = abs(prctile(svm_results_params{i}.acc_1vs1_pca_group(:),[75,25])-median(svm_results_params{i}.acc_1vs1_pca_group(:)));    
    
    accuracies_TMH_ind_1vsall(:,i) = svm_results_params{i}.acc_tmh_ind(:);
    accuracies_PCA_ind_1vsall(:,i) = svm_results_params{i}.acc_pca_ind(:);
    accuracies_TMH_group_1vsall(:,i) = svm_results_params{i}.acc_tmh_group(:);
    accuracies_PCA_group_1vsall(:,i) = svm_results_params{i}.acc_pca_group(:);
    
    accuracies_TMH_ind_ci_1vsall(1,i,:) = abs(prctile(svm_results_params{i}.acc_tmh_ind(:),[75,25])-median(svm_results_params{i}.acc_tmh_ind(:)));
    accuracies_PCA_ind_ci_1vsall(1,i,:) = abs(prctile(svm_results_params{i}.acc_pca_ind(:),[75,25])-median(svm_results_params{i}.acc_pca_ind(:)));
    accuracies_TMH_group_ci_1vsall(1,i,:) = abs(prctile(svm_results_params{i}.acc_tmh_group(:),[75,25])-median(svm_results_params{i}.acc_tmh_group(:)));
    accuracies_PCA_group_ci_1vsall(1,i,:) = abs(prctile(svm_results_params{i}.acc_pca_group(:),[75,25])-median(svm_results_params{i}.acc_pca_group(:)));    
end
%%

cc = colormap(brewermap(28,'Set1'));
c_map(1,:) = [42,138,72]./256;
c_map(2,:) = [235,103,78]./256;
lineProps(1).col{1} = c_map(1,:);
lineProps(2).col{1} = c_map(2,:);


%%

h9 = figure('Units','inches','Position',[3.7778 5.7500 9.5556 2.5833],'PaperPositionMode','auto'); hold on;

subplot(1,2,1)
mseb(1:(opts.dimension-1),median(accuracies_TMH_ind),accuracies_TMH_ind_ci,lineProps(1),1); hold on
mseb(1:(opts.dimension-1),median(accuracies_PCA_ind),accuracies_PCA_ind_ci,lineProps(2),1);
ylabel('Accuracy');
xlim([1,opts.dimension-1])
xlabel('Manifold dimensionality');

subplot(1,2,2)
mseb(1:(opts.dimension-1),median(accuracies_TMH_group),accuracies_TMH_group_ci,lineProps(1),1); hold on
mseb(1:(opts.dimension-1),median(accuracies_PCA_group),accuracies_PCA_group_ci,lineProps(2),1);
xlabel('Manifold dimensionality');
xlim([1,opts.dimension-1])

h9.Children(1).Children(6).FaceAlpha = 0.3;
h9.Children(1).Children(12).FaceAlpha = 0.3;
h9.Children(2).Children(6).FaceAlpha = 0.3;
h9.Children(2).Children(12).FaceAlpha = 0.3;


h9.Children(1).Children(1).LineWidth = 3;
h9.Children(1).Children(6).LineWidth = 3;
h9.Children(1).Children(4).LineWidth = 2;
h9.Children(1).Children(4).LineWidth = 2;
h9.Children(1).Children(5).LineWidth = 2;
h9.Children(1).Children(10).LineWidth = 2;
h9.Children(1).Children(11).LineWidth = 2;

h9.Children(2).Children(1).LineWidth = 3;
h9.Children(2).Children(6).LineWidth = 3;
h9.Children(2).Children(4).LineWidth = 2;
h9.Children(2).Children(4).LineWidth = 2;
h9.Children(2).Children(5).LineWidth = 2;
h9.Children(2).Children(10).LineWidth = 2;
h9.Children(2).Children(11).LineWidth = 2;

h9.Children(1).YLim = [0.38 1.02];
h9.Children(2).YLim = [0.38 1.02];
opt.FontSize = 20;
h9 = fancy_figure(h9, opt);

save2pdf('../Figures/Fig3/Manifold_dimensionality',h9,300)





%%

h9_2 = figure('Units','inches','Position',[3.7778 5.7500 9.5556 2.5833],'PaperPositionMode','auto'); hold on;

subplot(1,2,1)
mseb(1:(opts.dimension-1),median(accuracies_TMH_ind_1vsall),accuracies_TMH_ind_ci_1vsall,lineProps(1),1); hold on
mseb(1:(opts.dimension-1),median(accuracies_PCA_ind_1vsall),accuracies_PCA_ind_ci_1vsall,lineProps(2),1);
ylabel('Accuracy');
xlim([1,opts.dimension-1])
xlabel('Manifold dimensionality');

subplot(1,2,2)
mseb(1:(opts.dimension-1),median(accuracies_TMH_group_1vsall),accuracies_TMH_group_ci_1vsall,lineProps(1),1); hold on
mseb(1:(opts.dimension-1),median(accuracies_PCA_group_1vsall),accuracies_PCA_group_ci_1vsall,lineProps(2),1);
xlabel('Manifold dimensionality');
xlim([1,opts.dimension-1])

h9_2.Children(1).Children(6).FaceAlpha = 0.3;
h9_2.Children(1).Children(12).FaceAlpha = 0.3;
h9_2.Children(2).Children(6).FaceAlpha = 0.3;
h9_2.Children(2).Children(12).FaceAlpha = 0.3;


h9_2.Children(1).Children(1).LineWidth = 3;
h9_2.Children(1).Children(6).LineWidth = 3;
h9_2.Children(1).Children(4).LineWidth = 2;
h9_2.Children(1).Children(4).LineWidth = 2;
h9_2.Children(1).Children(5).LineWidth = 2;
h9_2.Children(1).Children(10).LineWidth = 2;
h9_2.Children(1).Children(11).LineWidth = 2;

h9_2.Children(2).Children(1).LineWidth = 3;
h9_2.Children(2).Children(6).LineWidth = 3;
h9_2.Children(2).Children(4).LineWidth = 2;
h9_2.Children(2).Children(4).LineWidth = 2;
h9_2.Children(2).Children(5).LineWidth = 2;
h9_2.Children(2).Children(10).LineWidth = 2;
h9_2.Children(2).Children(11).LineWidth = 2;

h9_2.Children(1).YLim = [0.1 1.02];
h9_2.Children(2).YLim = [0.1 1.02];
opt.FontSize = 20;
h9_2 = fancy_figure(h9_2, opt);

save2pdf('../Figures/Fig3/Manifold_dimensionality_1vs_all',h9_2,300)

%%

rng(0)
h10 = figure('Units','inches','Position',[0 0 8 7],'PaperPositionMode','auto'); hold on; 
methods = {'acc_1vs1_tmh_ind','acc_1vs1_pca_ind','acc_1vs1_tmh_group','acc_1vs1_pca_group'};
fc = [c_map(1:2,:);c_map(1:2,:)];
x = [1,2,3,4];
alphas = [.4,.4,.4,.4];
for i = 1:6
    subplot(2,3,i),hold on
    for j = 1:4
        
        
        data = eval(sprintf('svm_results_params{7}.%s(:,ord_1v1(i))',methods{j}));             
        Violin(data,x(j),'ViolinColor',fc(j,:),'ViolinAlpha',alphas(j),'EdgeColor',fc(j,:),'BoxColor',[0,0,0],'MedianColor',fc(j,:),'ShowMean',true);
    end              
    xlim([0.5,4.5])
    xticks([])    
    ylim([0.4,1.05])
    yticks([0.5,0.75,1])
    plot([mean(x),mean(x)],[0.5,1],'--','Color',[0,0,0.5],'LineWidth',0.5)
    
end
h10 = fancy_figure(h10,opt);

save2pdf('../Figures/Fig3/SVM_accuracies_7d',h10,300)

%%
rng(0)
h10 = figure('Units','inches','Position',[0 0 8 7],'PaperPositionMode','auto'); hold on; 
methods = {'acc_1vs1_tmh_ind','acc_1vs1_pca_ind','acc_1vs1_tmh_group','acc_1vs1_pca_group'};
fc = [c_map(1:2,:);c_map(1:2,:)];
x = [1,2,3,4];
alphas = [.4,.4,.4,.4];
for i = 1:6
    subplot(2,3,i),hold on
    for j = 1:4             
        data = eval(sprintf('svm_results_params{3}.%s(:,ord_1v1(i))',methods{j}));        
        Violin(data,x(j),'ViolinColor',fc(j,:),'ViolinAlpha',alphas(j),'EdgeColor',fc(j,:),'BoxColor',[0,0,0],'MedianColor',fc(j,:),'ShowMean',true);
    end              
    xlim([0.5,4.5])
    xticks([])    
    ylim([0.4,1.05])
    yticks([0.5,0.75,1])
    plot([mean(x),mean(x)],[0.5,1],'--','Color',[0,0,0.5],'LineWidth',0.5)
    
end
h10 = fancy_figure(h10,opt);

save2pdf('../Figures/Fig3/SVM_accuracies_3d',h10,300)
%%

rng(0)
h10 = figure('Units','inches','Position',[0 0 8 7],'PaperPositionMode','auto'); hold on; 
methods = {'C_norm_tmh_ind','C_norm_pca_ind','C_norm_tmh_group','C_norm_pca_group'};
fc = [c_map(1:2,:);c_map(1:2,:)];
x = [1,2,3,4];
alphas = [.4,.4,.4,.4];
for i = 1:4
    subplot(2,3,i),hold on
    
    for j = 1:4        
        data = squeeze(eval(sprintf('svm_results_params{7}.%s(ord_1vsall(i),ord_1vsall(i),:)',methods{j}))); 
        Violin(data,x(j),'ViolinColor',fc(j,:),'ViolinAlpha',alphas(j),'EdgeColor',fc(j,:),'BoxColor',[0,0,0],'MedianColor',fc(j,:),'ShowMean',true);
    end              
    xlim([0.5,4.5])
    xticks([])    
    ylim([-0.05,1.05])
    yticks([0,0.25,0.5,0.75,1.])
    plot([mean(x),mean(x)],[0.,1],'--','Color',[0,0,0.5],'LineWidth',0.5)
    
end
h10 = fancy_figure(h10,opt);

save2pdf('../Figures/Fig3/SVM_accuracies_1vsall_7d',h10,300)

%%
rng(0)
h10 = figure('Units','inches','Position',[0 0 8 7],'PaperPositionMode','auto'); hold on; 
methods = {'C_norm_tmh_ind','C_norm_pca_ind','C_norm_tmh_group','C_norm_pca_group'};
fc = [c_map(1:2,:);c_map(1:2,:)];
x = [1,2,3,4];
alphas = [.4,.4,.4,.4];
for i = 1:4
    subplot(2,3,i),hold on
    
    for j = 1:4        
        
        data = squeeze(eval(sprintf('svm_results_params{3}.%s(ord_1vsall(i),ord_1vsall(i),:)',methods{j}))); 
        Violin(data,x(j),'ViolinColor',fc(j,:),'ViolinAlpha',alphas(j),'EdgeColor',fc(j,:),'BoxColor',[0,0,0],'MedianColor',fc(j,:),'ShowMean',true);
    end              
    xlim([0.5,4.5])
    xticks([])    
    ylim([-0.05,1.05])
    yticks([0,0.25,0.5,0.75,1.])
    plot([mean(x),mean(x)],[0.,1],'--','Color',[0,0,0.5],'LineWidth',0.5)
    
end
h10 = fancy_figure(h10,opt);

save2pdf('../Figures/Fig3/SVM_accuracies_1vsall_3d',h10,300)

%%

rng(0)
h10_2 = figure('Units','inches','Position',[0 0 2 5],'PaperPositionMode','auto'); hold on; 
methods = {'C_norm_tmh_ind','C_norm_pca_ind','C_norm_tmh_group','C_norm_pca_group'};
titles = {'TMH ind','PCA ind','TMH group','PCA group'};
sleepstages = {'A','N1','N2','N3'};
fc = [c_map(1:2,:);c_map(1:2,:)];
x = [1,2,3,4];
alphas = [.4,.4,.4,.4];
for j =1:4
    subplot(4,1,j)
    
    data = squeeze(mean(eval(sprintf('svm_results_params{7}.%s(ord_1vsall,ord_1vsall,:)',methods{j})),3)); 
    h = heatmap(sleepstages,sleepstages,data,'Colormap',gray);
    h.CellLabelFormat = '%.2f';
    caxis(h,[0 1]);
    title(titles{j})    
    
end
h10_2 = fancy_figure(h10_2,opt);

saveas(h10_2,'../Figures/Fig3/SVM_confusion_7d.png')


%%
rng(0)
h10_2 = figure('Units','inches','Position',[0 0 2 5],'PaperPositionMode','auto'); hold on; 
methods = {'C_norm_tmh_ind','C_norm_pca_ind','C_norm_tmh_group','C_norm_pca_group'};
titles = {'TMH ind','PCA ind','TMH group','PCA group'};
sleepstages = {'A','N1','N2','N3'};
fc = [c_map(1:2,:);c_map(1:2,:)];
x = [1,2,3,4];
alphas = [.4,.4,.4,.4];
for j =1:4
    subplot(4,1,j)        
    data = squeeze(mean(eval(sprintf('svm_results_params{3}.%s(ord_1vsall,ord_1vsall,:)',methods{j})),3)); 
    h = heatmap(sleepstages,sleepstages,data,'Colormap',gray);
    h.CellLabelFormat = '%.2f';
    caxis(h,[0 1]);
    title(titles{j})    
end
h10_2 = fancy_figure(h10_2,opt);

saveas(h10_2,'../Figures/Fig3/SVM_confusion_3d.png')

%%
h11 = figure();
scatter(1,1,100,fc(1,:),'filled','MarkerFaceAlpha',0.9)
hold on
scatter(1,2,100,fc(2,:),'filled','MarkerFaceAlpha',0.9)
legend('Intrinsic manifold','Principal components')
h11 = fancy_figure(h11,opt);
save2pdf('../Figures/Fig3/legend',h11,300)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               FIGURE 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLot ROC curves

h12 = figure('Units','inches','Position',[0 0 12 7],'PaperPositionMode','auto'); hold on; 

lineprops1.col{1} = c_map(1,:);

lineprops2.col{1} = c_map(2,:);
lineprops2.style = '-.';

for i = 1:6
subplot(2,3,i);
eval(sprintf('mseb(TMH_all.line%0.1i.XData, TMH_all.line%0.1i.YData,TMH_all.line%0.1i.YData_err./sqrt(18),lineprops1, 1);',i,i,i));hold on;
eval(sprintf('mseb(PCA_all.line%0.1i.XData, PCA_all.line%0.1i.YData,PCA_all.line%0.1i.YData_err./sqrt(18),lineprops2, 1);',i,i,i));hold on;
text(0.45, 0.2, ['IM-AUC=', num2str(mean(eval(sprintf('TMH_all.auc%0.1i',i))), '%5.2f')], 'FontSize', 13, 'Color', c_map(1,:));
text(0.45, 0.1, ['PCA-AUC=', num2str(mean(eval(sprintf('PCA_all.auc%0.1i',i))), '%5.2f')], 'FontSize', 13, 'Color', c_map(2,:));

xlabel({'(1-Specificity)'});
ylabel({'Sensitivity'});
xlim([-0.01,1.01]);
ylim([-0.02,1.03]);
end
for i = 1:6
    h12.Children(i).Children(8).FaceAlpha = 0.3;
    h12.Children(i).Children(14).FaceAlpha = 0.3;
end

opt.FontSize = 20;
h12 = fancy_figure(h12, opt);
save2pdf('../Figures/Fig4/ROC_curves',h12,300)
%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               FIGURE S1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure of Individual TMH low-dim manifold and PCA, and the corresponding 
% group embeedings

opt.FontSize = 20;
for i = 1:Nsub
    h8 = figure('Units','inches','Position',[0 0 4 4],'PaperPositionMode','auto'); hold on;
    ids_sub = (1+(87*4)*(i-1)):(87*4)*i;
    scatter3(group_data.TMH(ids_sub,2),group_data.TMH(ids_sub,1),group_data.TMH(ids_sub,3),...
        size_,colorss(abs(group_data.stages(ids_sub))+1,:),'o','filled','MarkerFaceAlpha',alpha_)
 
    xlim([-0.12,0.12]),ylim([-0.12,0.12]),zlim([-0.12,0.12])
    xticks([-0.1, 0, 0.1]),yticks([-0.1, 0, 0.1]), zticks([-0.1, 0, 0.1])
    xticklabels([])
    yticklabels([])
    zticklabels([])
    
    grid on
    ax = gca;
    box on    
    ax.LineWidth = 1.5;
    h8 = fancy_figure(h8, opt);
    view(3)   
    save2pdf(sprintf('../Figures/FigS1/Embeddings_TMH_%i',i),h8,300)    
    close
    
    % PCA
    h8 = figure('Units','inches','Position',[0 0 4 4],'PaperPositionMode','auto'); hold on;
    ids_sub = (1+(87*4)*(i-1)):(87*4)*i;
    scatter3(group_data.PCA(ids_sub,2),group_data.PCA(ids_sub,1),group_data.PCA(ids_sub,3),...
        size_*2,colorss(abs(group_data.stages(ids_sub))+1,:),'o','filled','MarkerFaceAlpha',alpha_)
 
    xticks([-0.1, 0, 0.1]),yticks([-0.1, 0, 0.1]), zticks([-0.1, 0, 0.1])
    xticklabels([])
    yticklabels([])
    zticklabels([])
    
    grid on
    ax = gca;
    box on    
    ax.LineWidth = 1.5;
    h8 = fancy_figure(h8, opt);
    view(3)   
    save2pdf(sprintf('../Figures/FigS1/Embeddings_PCA_%i',i),h8,300)   
    close
    
end

%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               FIGURE S1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts_plot = opts;
opts_plot.dimension = 7;
[group_data_plot,~] = compute_TMH_group_rotation(group_data,ind_data,opts_plot,Nsub,87);    
[group_data_plot,~] = compute_PCA_group_rotation(group_data_plot,ind_data,opts_plot,Nsub,87);  

opt.FontSize = 20;

for i = 1:2
    h8 = figure('Units','inches','Position',[0 0 4 4],'PaperPositionMode','auto'); hold on;
    scatter3(group_data_plot.TMH(:,i*3-1),group_data_plot.TMH(:,i*3-2),group_data_plot.TMH(:,i*3),...
        size_,colorss(abs(group_data_plot.stages)+1,:),'o','filled','MarkerFaceAlpha',alpha_)
    
    xlim([-0.10,0.10]),ylim([-0.1,0.1]),zlim([-0.1,0.1])
    xticks([-0.10, 0, 0.1]),yticks([-0.1, 0, 0.1]), zticks([-0.1, 0, 0.1])
    grid on,ax = gca;box on,ax.LineWidth = 1.5;
    h8 = fancy_figure(h8, opt);
    xticklabels([])
    yticklabels([])
    zticklabels([])
    view(3)
    save2pdf(sprintf('../Figures/FigS3/Embeddings_TMH_all_dims_%i',i),h8,300)
    close
    
    h8 = figure('Units','inches','Position',[0 0 4 4],'PaperPositionMode','auto'); hold on;      
    scatter3(group_data_plot.PCA(:,i*3-2),group_data_plot.PCA(:,i*3-1),group_data_plot.PCA(:,i*3),...
        size_,colorss(abs(group_data_plot.stages)+1,:),'o','filled','MarkerFaceAlpha',alpha_)
   
    xticks([-20, 0, 20]),yticks([-10, 0, 10]), zticks([-10, 0, 10])
    grid on,ax = gca;box on,ax.LineWidth = 1.5;
    h8 = fancy_figure(h8, opt);
    view(3)
    xticklabels([])
    yticklabels([])
    zticklabels([])
    save2pdf(sprintf('../Figures/FigS3/Embeddings_PCA_all_dims_%i',i),h8,300)
end