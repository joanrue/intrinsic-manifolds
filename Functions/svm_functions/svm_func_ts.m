function [acc_ts,C_ts,C_norm_ts] = svm_func_ts(data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if ~isfield(data, 'subj_vector')
        % ind analysis - 6-fold cross-validation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare data and normalize        
        n_subjects = length(data);
        
        acc_ts = zeros(n_subjects);
        C_ts = zeros(4,4,n_subjects);
        C_norm_ts = zeros(4,4,n_subjects);
       
        for s = 1:n_subjects
            [~,~,labels] = unique(data{s}.stages);
            data_ts = zscore(data{s}.ts');            
            num_labels = max(labels);   
            
            prob_ts = zeros(length(labels),num_labels);
            
            for k=1:num_labels
                % train svm
                model_ts = fitcsvm(data_ts,labels==k,'Standardize',true,'CrossVal','on');                
            
                % fit posterior distribution 
                [score_model_ts] = fitSVMPosterior(model_ts);
                                
                % get posterior test prob
                [~, post_prob_ts] = kfoldPredict(score_model_ts);
                
                prob_ts(:,k) = post_prob_ts(:,2);                
            end
            
            
            [~,pred_ts] = max(prob_ts,[],2);    
            acc_ts(s) = sum(pred_ts == labels) ./ numel(labels);
            
            C_ts(:,:,s) = confusionmat(labels, pred_ts);            
            C_norm_ts(:,:,s) = C_ts(:,:,s) ./ sum(C_ts(:,:,s),2);                  
            
        end

                  
        

       
    else
    
        % group analysis - leave-one-subject-out cross-validation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare data and normalize
        [~,~,subjects] = unique(data.subj_vector);
        subject_labels = 1:max(subjects);
        [~,~,labels] = unique(data.stages);
        data_ts = zscore(data.ts');        
        num_labels = max(labels);                   

        acc_ts = zeros(length(subject_labels),1);
        C_ts = zeros(num_labels,num_labels,length(subject_labels));
        C_norm_ts = zeros(num_labels,num_labels,length(subject_labels));
       
        for s = 1:length(subject_labels)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % split training/testing
            train_data_ts = data_ts(subjects~=subject_labels(s),:);  test_data_ts = data_ts(subjects==subject_labels(s),:);            
            train_label    = labels(subjects~=subject_labels(s));      test_label    = labels(subjects==subject_labels(s));    
            num_test = length(test_label);        
            prob_ts = zeros(num_test,num_labels);            
            for k=1:num_labels
                % train svm
                svm_model_ts = compact(fitcsvm(train_data_ts,double(train_label==k)));           
                
                % fit posterior distribution 
                model_ts = fitPosterior(svm_model_ts ,train_data_ts,double(train_label==k));
                
                % get posterior test prob
                [~,post_prob_ts] = predict(model_ts,test_data_ts);                
                
                prob_ts(:,k) = post_prob_ts(:,2);
                
            end

            [~,pred_ts] = max(prob_ts,[],2);                
            acc_ts(s) = sum(pred_ts == test_label) ./ numel(test_label);            
            C_ts(:,:,s) = confusionmat(test_label, pred_ts);            
            C_norm_ts(:,:,s) = C_ts(:,:,s) ./ sum(C_ts(:,:,s),2);      
            
        end
    end
end