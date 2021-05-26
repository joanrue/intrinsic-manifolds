function [acc_tmh,C_tmh,C_norm_tmh,acc_pca,C_pca,C_norm_pca,conds] = svm_func(data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    warning('off')
    if ~isfield(data, 'subj_vector')
        % ind analysis - 6-fold cross-validation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare data and normalize        
        n_subjects = length(data);
        
        acc_tmh = zeros(n_subjects,1);
        C_tmh = zeros(4,4,n_subjects);
        C_norm_tmh = zeros(4,4,n_subjects);

        acc_pca = zeros(n_subjects,1);
        C_pca = zeros(4,4,n_subjects);
        C_norm_pca = zeros(4,4,n_subjects);

        for s = 1:n_subjects
            [C,~,labels] = unique(data{s}.stages);
            data_tmh = zscore(data{s}.TMH);
            data_pca = zscore(data{s}.PCA);
            num_labels = max(labels);   
            
            prob_tmh = zeros(length(labels),num_labels);
            prob_pca = zeros(length(labels),num_labels);
            
            for k=1:num_labels
                % train svm
                model_tmh = fitcsvm(data_tmh,labels==k,'Standardize',true,'CrossVal','on');
                model_pca = fitcsvm(data_pca,labels==k,'Standardize',true,'CrossVal','on');           
            
                % fit posterior distribution 
                [score_model_tmh] = fitSVMPosterior(model_tmh);
                [score_model_pca] = fitSVMPosterior(model_pca);
                
                % get posterior test prob
                [~, post_prob_tmh] = kfoldPredict(score_model_tmh);
                [~, post_prob_pca] = kfoldPredict(score_model_pca);
                            
                prob_tmh(:,k) = post_prob_tmh(:,2);
                prob_pca(:,k) = post_prob_pca(:,2);                
            end
            
            
            [~,pred_tmh] = max(prob_tmh,[],2);    
            [~,pred_pca] = max(prob_pca,[],2);    
            acc_tmh(s) = sum(pred_tmh == labels) ./ numel(labels);
            acc_pca(s) = sum(pred_pca == labels) ./ numel(labels);
            C_tmh(:,:,s) = confusionmat(labels, pred_tmh);
            C_pca(:,:,s) = confusionmat(labels, pred_pca);
            C_norm_tmh(:,:,s) = C_tmh(:,:,s) ./ sum(C_tmh(:,:,s),2);      
            C_norm_pca(:,:,s) = C_pca(:,:,s) ./ sum(C_pca(:,:,s),2);      
            
        end

        conds = {};
        for k=1:num_labels
            conds{k} = C(k);
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

        acc_tmh = zeros(length(subject_labels),1);
        C_tmh = zeros(num_labels,num_labels,length(subject_labels));
        C_norm_tmh = zeros(num_labels,num_labels,length(subject_labels));

        acc_pca = zeros(length(subject_labels),1);
        C_pca = zeros(num_labels,num_labels,length(subject_labels));
        C_norm_pca = zeros(num_labels,num_labels,length(subject_labels));

        for s = 1:length(subject_labels)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % split training/testing
            train_data_tmh = data_tmh(subjects~=subject_labels(s),:);  test_data_tmh = data_tmh(subjects==subject_labels(s),:);
            train_data_pca = data_pca(subjects~=subject_labels(s),:);  test_data_pca = data_pca(subjects==subject_labels(s),:);
            train_label    = labels(subjects~=subject_labels(s));      test_label    = labels(subjects==subject_labels(s));    
            num_test = length(test_label);        
            prob_tmh = zeros(num_test,num_labels);
            prob_pca = zeros(num_test,num_labels);
            for k=1:num_labels
                % train svm
                svm_model_tmh = compact(fitcsvm(train_data_tmh,double(train_label==k)));           
                svm_model_pca = compact(fitcsvm(train_data_pca,double(train_label==k)));           
                % fit posterior distribution 
                model_tmh = fitPosterior(svm_model_tmh ,train_data_tmh,double(train_label==k));
                model_pca = fitPosterior(svm_model_pca ,train_data_pca,double(train_label==k));
                % get posterior test prob
                [~,post_prob_tmh] = predict(model_tmh,test_data_tmh);
                [~,post_prob_pca] = predict(model_pca,test_data_pca);
                prob_tmh(:,k) = post_prob_tmh(:,2);
                prob_pca(:,k) = post_prob_pca(:,2);
            end

            [~,pred_tmh] = max(prob_tmh,[],2);    
            [~,pred_pca] = max(prob_pca,[],2);    
            acc_tmh(s) = sum(pred_tmh == test_label) ./ numel(test_label);
            acc_pca(s) = sum(pred_pca == test_label) ./ numel(test_label);
            C_tmh(:,:,s) = confusionmat(test_label, pred_tmh);
            C_pca(:,:,s) = confusionmat(test_label, pred_pca);
            C_norm_tmh(:,:,s) = C_tmh(:,:,s) ./ sum(C_tmh(:,:,s),2);      
            C_norm_pca(:,:,s) = C_pca(:,:,s) ./ sum(C_pca(:,:,s),2);      
        end
    end
    conds = {};
    for k=1:num_labels
        conds{k} = C(k);
    end
end