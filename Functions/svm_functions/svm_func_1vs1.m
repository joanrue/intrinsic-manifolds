function [acc_tmh,acc_pca,conds] = svm_func_1vs1(data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    warning('off')
    if ~isfield(data, 'subj_vector')
        % ind analysis - 6-fold cross-validation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare data and normalize        
        n_subjects = length(data);
        
        acc_tmh = zeros(n_subjects,6);       
        acc_pca = zeros(n_subjects,6);
        
        for s = 1:n_subjects
            [C,IA,labels] = unique(data{s}.stages);
            data_tmh = zscore(data{s}.TMH);
            data_pca = zscore(data{s}.PCA);
            num_labels = max(labels);   
            
            
            cond = 0;           
            for k1=1:num_labels
                for k2=k1+1:num_labels
                    cond = cond+1;
                    subset =  find ((labels==k1) | (labels==k2));
                    % train svm
                    model_tmh = fitcsvm(data_tmh(subset,:),labels(subset),'Standardize',true,'CrossVal','on');
                    model_pca = fitcsvm(data_pca(subset,:),labels(subset),'Standardize',true,'CrossVal','on');           

                    % fit posterior distribution 
                    score_model_tmh = fitSVMPosterior(model_tmh);
                    score_model_pca = fitSVMPosterior(model_pca);

                    % get posterior test prob
                    pred_tmh = kfoldPredict(score_model_tmh);
                    pred_pca = kfoldPredict(score_model_pca);
                    
                    acc_tmh(s,cond) = sum(pred_tmh == labels(subset)) ./ numel(subset);
                    acc_pca(s,cond) = sum(pred_pca == labels(subset)) ./ numel(subset);
                    
                end
            end                                    
            
        end
                  
        conds = {};
        cond = 0;
        for k1=1:num_labels
            for k2=k1+1:num_labels
                if ~isempty(k2)
                    cond = cond+1;
                    conds{cond} = [C(k1),C(k2)];
                end
            end
        end
       
    else
    
        % group analysis - leave-one-subject-out cross-validation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare data and normalize
        [~,~,subjects] = unique(data.subj_vector);
        subject_labels = 1:max(subjects);
        [C,~,labels] = unique(data.stages);
        data_tmh = zscore(data.TMH);
        data_pca = zscore(data.PCA);
        num_labels = max(labels);                   

        acc_tmh = zeros(length(subject_labels),6);       
        acc_pca = zeros(length(subject_labels),6);
        
        for s = 1:length(subject_labels)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % split training/testing
            train_data_tmh = data_tmh(subjects~=subject_labels(s),:);  test_data_tmh = data_tmh(subjects==subject_labels(s),:);
            train_data_pca = data_pca(subjects~=subject_labels(s),:);  test_data_pca = data_pca(subjects==subject_labels(s),:);
            train_label    = labels(subjects~=subject_labels(s));      test_label    = labels(subjects==subject_labels(s));                
            cond = 0;           
            for k1=1:num_labels
                for k2=k1+1:num_labels
                    cond = cond+1;
                    subset_train =  find ((train_label==k1) | (train_label==k2));
                    subset_test =  find ((test_label==k1) | (test_label==k2));
                    % train svm
                    svm_model_tmh = compact(fitcsvm(train_data_tmh(subset_train,:),double(train_label(subset_train))));           
                    svm_model_pca = compact(fitcsvm(train_data_pca(subset_train,:),double(train_label(subset_train))));           
                    % fit posterior distribution 
                    model_tmh = fitPosterior(svm_model_tmh ,train_data_tmh(subset_train,:),train_label(subset_train));
                    model_pca = fitPosterior(svm_model_pca ,train_data_pca(subset_train,:),train_label(subset_train));
                    % get posterior test prob                    
                    pred_tmh = predict(model_tmh,test_data_tmh(subset_test,:));
                    pred_pca = predict(model_pca,test_data_pca(subset_test,:));
                    
                    acc_tmh(s,cond) = sum(pred_tmh == test_label(subset_test)) ./ numel(subset_test);
                    acc_pca(s,cond) = sum(pred_pca == test_label(subset_test)) ./ numel(subset_test);
                    
                end
            end
        end
    end
    conds = {};
    cond = 0;
    for k1=1:num_labels
        for k2=k1+1:num_labels
            if ~isempty(k2)
                cond = cond+1;
                conds{cond} = [C(k1),C(k2)];
            end
        end
    end
end