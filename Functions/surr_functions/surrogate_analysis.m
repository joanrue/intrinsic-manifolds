
%%  Add external functions to path 

addpath('../external')
addpath('../svm_functions')
addpath('../int_manifold_functions')

realization = str2num(getenv('SLURM_ARRAY_TASK_ID'));
roundidx    = str2num(getenv('ROUND_IDX'));

%% Load Data 

path2data = ['..' filesep '..' filesep 'Data'];

load([path2data filesep 'sleep_scoring_N0_subs'])
load([path2data filesep 'TC_N0_subs_REV'])

%% Define Parameters  

% Relaxed Minimum Spanning Tree params.
opts.gamma = 3;%2;
opts.maxKNN = 5;%5;
opts.distance = 'cosine'; 
% Laplacian Eigenmaps params.
opts.weighted = 0;
opts.laplacian = 'combinatorial'; % Graph Laplacian can be combinatorial, normalized

opts.TSmax = 87;    % Max num of TRs per each condition.
opts.TSmin = 87;    % Min num of TRs per each condition (18 subjects with
% a minimum of 87 samples per condition. 12 subjects
% with a minimum of 144 samples per condition).

Nsub = 18;

%% Surrogate analysis

rng('shuffle');

cond = true;
svm_results_params = cell(10,1);
for i = [3,7]
    fprintf('SVM using %0.2i dimensions\n', i);
    
    opts_svm = opts;
    opts_svm.dimension = i;
    
    %%
    while cond == true
        try
            [ind_data,~,group_data] = preprocess(TC_all_states,inst_states_all,opts_svm.TSmin,opts_svm.TSmax,true);
            for nsub = 1:Nsub                
                [~,ind_data{1,nsub}.PCA,~]=pca(ind_data{1,nsub}.ts','NumComponents',opts_svm.dimension);
            end            
            [~, group_data.PCA, ~]=pca(group_data.ts','NumComponents',opts_svm.dimension);
            ind_data_svm = compute_ind_manifold(ind_data,opts_svm,false);
            [group_data_svm,~] = compute_group_manifold(group_data,ind_data_svm,opts_svm,Nsub,87);
            cond = false;
        catch
            warning('Singular surrogate')
        end
    end
           
    % 1vs1 SVM - TMH
    [acc_tmh,acc_pca,~] = svm_func_1vs1(ind_data_svm);
    svm_results_params{i}.acc_1vs1_tmh_ind = acc_tmh;
    svm_results_params{i}.acc_1vs1_pca_ind = acc_pca;
    
    % Multi-class SVM - TMH
    [acc_tmh,~,C_norm_tmh,acc_pca,~,C_norm_pca,~] = svm_func(ind_data_svm);
    svm_results_params{i}.acc_tmh_ind = acc_tmh;            
    svm_results_params{i}.C_norm_tmh_ind = C_norm_tmh;
    svm_results_params{i}.acc_pca_ind = acc_pca;            
    svm_results_params{i}.C_norm_pca_ind = C_norm_pca;        
    
    % 1vs1 SVM - PCA
    [acc_tmh,acc_pca,conds_1vs1] = svm_func_1vs1(group_data_svm);
    svm_results_params{i}.acc_1vs1_tmh_group = acc_tmh;
    svm_results_params{i}.acc_1vs1_pca_group = acc_pca;
    
    % Multi-class SVM - PCA    
    [acc_tmh,~,C_norm_tmh,acc_pca,~,C_norm_pca,conds_multiclass] = svm_func(group_data_svm);
    svm_results_params{i}.acc_tmh_group = acc_tmh;            
    svm_results_params{i}.C_norm_tmh_group = C_norm_tmh;
    svm_results_params{i}.acc_pca_group = acc_pca;            
    svm_results_params{i}.C_norm_pca_group = C_norm_pca;   
    
    
end

[acc_ts,~] = svm_func_1vs1_ts(ind_data);
svm_results_params{1}.acc_1vs1_ind = acc_ts;            

[acc_ts,~,C_norm_ts] = svm_func_ts(ind_data);
svm_results_params{1}.acc_ind = acc_ts;            
svm_results_params{1}.C_norm_ind = C_norm_ts;

[acc_ts,~] = svm_func_1vs1_ts(group_data);
svm_results_params{1}.acc_1vs1_group = acc_ts;            

[acc_ts,~,C_norm_ts] = svm_func_ts(group_data);
svm_results_params{1}.acc_group = acc_ts;            
svm_results_params{1}.C_norm_group = C_norm_ts;

save(sprintf('../../Data/Surrogates/all_surr/SVM_results_surr_%0.3i',roundidx+realization),'svm_results_params');