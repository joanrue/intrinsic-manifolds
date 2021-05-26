addpath external
path2data = ['..' filesep '..' filesep 'Data'];
sur_dir = ['path2data' filesep 'Surrogates' filesep 'all_surr'];
d = dir_ls(sur_dir);

num_sur = length(d);

svm_results_surr.acc_1vs1_tmh_ind = zeros(2,18,6,num_sur);
svm_results_surr.acc_1vs1_pca_ind = zeros(2,18,6,num_sur);
svm_results_surr.acc_tmh_ind = zeros(2,18,num_sur);
svm_results_surr.C_norm_tmh_ind = zeros(2,4,4,18,num_sur);
svm_results_surr.acc_pca_ind = zeros(2,18,num_sur);
svm_results_surr.C_norm_pca_ind = zeros(2,4,4,18,num_sur);
svm_results_surr.acc_1vs1_tmh_group = zeros(2,18,6,num_sur);
svm_results_surr.acc_1vs1_pca_group = zeros(2,18,6,num_sur);
svm_results_surr.acc_tmh_group = zeros(2,18,num_sur);
svm_results_surr.C_norm_tmh_group = zeros(2,4,4,18,num_sur);
svm_results_surr.acc_pca_group = zeros(2,18,num_sur);
svm_results_surr.C_norm_pca_group = zeros(2,4,4,18,num_sur);

for sur = 1:num_sur
    load([sur_dir filesep d(sur).name])
    if ~isempty(svm_results_params{3})
        svm_results_surr.acc_1vs1_tmh_ind(1,:,:,sur) = svm_results_params{3}.acc_1vs1_tmh_ind;
        svm_results_surr.acc_1vs1_pca_ind(1,:,:,sur) = svm_results_params{3}.acc_1vs1_pca_ind;
        svm_results_surr.acc_tmh_ind(1,:,sur) = svm_results_params{3}.acc_tmh_ind;
        svm_results_surr.C_norm_tmh_ind(1,:,:,:,sur) = svm_results_params{3}.C_norm_tmh_ind;
        svm_results_surr.acc_pca_ind(1,:,sur) = svm_results_params{3}.acc_pca_ind;
        svm_results_surr.C_norm_pca_ind(1,:,:,:,sur) = svm_results_params{3}.C_norm_pca_ind;
        svm_results_surr.acc_1vs1_tmh_group(1,:,:,sur) = svm_results_params{3}.acc_1vs1_tmh_group;
        svm_results_surr.acc_1vs1_pca_group(1,:,:,sur) = svm_results_params{3}.acc_1vs1_pca_group;
        svm_results_surr.acc_tmh_group(1,:,sur) = svm_results_params{3}.acc_tmh_group;
        svm_results_surr.C_norm_tmh_group(1,:,:,:,sur) = svm_results_params{3}.C_norm_tmh_group;
        svm_results_surr.acc_pca_group(1,:,sur) = svm_results_params{3}.acc_pca_group;
        svm_results_surr.C_norm_pca_group(1,:,:,:,sur) = svm_results_params{3}.C_norm_pca_group;
    end
    if ~isempty(svm_results_params{10})
        svm_results_surr.acc_1vs1_tmh_ind(2,:,:,sur) = svm_results_params{10}.acc_1vs1_tmh_ind;
        svm_results_surr.acc_1vs1_pca_ind(2,:,:,sur) = svm_results_params{10}.acc_1vs1_pca_ind;
        svm_results_surr.acc_tmh_ind(2,:,sur) = svm_results_params{10}.acc_tmh_ind;
        svm_results_surr.C_norm_tmh_ind(2,:,:,:,sur) = svm_results_params{10}.C_norm_tmh_ind;
        svm_results_surr.acc_pca_ind(2,:,sur) = svm_results_params{10}.acc_pca_ind;
        svm_results_surr.C_norm_pca_ind(2,:,:,:,sur) = svm_results_params{10}.C_norm_pca_ind;
        svm_results_surr.acc_1vs1_tmh_group(2,:,:,sur) = svm_results_params{10}.acc_1vs1_tmh_group;
        svm_results_surr.acc_1vs1_pca_group(2,:,:,sur) = svm_results_params{10}.acc_1vs1_pca_group;
        svm_results_surr.acc_tmh_group(2,:,sur) = svm_results_params{10}.acc_tmh_group;
        svm_results_surr.C_norm_tmh_group(2,:,:,:,sur) = svm_results_params{10}.C_norm_tmh_group;
        svm_results_surr.acc_pca_group(2,:,sur) = svm_results_params{10}.acc_pca_group;
        svm_results_surr.C_norm_pca_group(2,:,:,:,sur) = svm_results_params{10}.C_norm_pca_group;
    end
    clear svm_results_params
end

save([path2data filesep 'Surrogates' filesep 'svm_results_surr'],'svm_results_surr');