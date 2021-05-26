function [acc_ts,conds] = svm_func_1vs1_ts(data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    warning('off')
    if ~isfield(data, 'subj_vector')
        % ind analysis - 6-fold cross-validation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prepare data and normalize        
        n_subjects = length(data);
        
        acc_ts = zeros(n_subjects,6);               
        
        for s = 1:n_subjects
            [C,IA,labels] = unique(data{s}.stages);
            data_ts = zscore(data{s}.ts');            
            num_labels = max(labels);   
            
            
            cond = 0;           
            for k1=1:num_labels
                for k2=k1+1:num_labels
                    cond = cond+1;
                    subset =  find ((labels==k1) | (labels==k2));
                    % train svm
                    model_ts = fitcsvm(data_ts(subset,:),labels(subset),'Standardize',true,'CrossVal','on');                    

                    % fit posterior distribution 
                    score_model_ts = fitSVMPosterior(model_ts);                 

                    % get posterior test prob
                    pred_ts = kfoldPredict(score_model_ts);                    
                    
                    acc_ts(s,cond) = sum(pred_ts == labels(subset)) ./ numel(subset);                    
                    
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
        data_ts = zscore(data.ts');        
        num_labels = max(labels);                   

        acc_ts = zeros(length(subject_labels),6);               
        
        for s = 1:length(subject_labels)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % split training/testing
            train_data_ts = data_ts(subjects~=subject_labels(s),:);  test_data_ts = data_ts(subjects==subject_labels(s),:);            
            train_label    = labels(subjects~=subject_labels(s));      test_label    = labels(subjects==subject_labels(s));                
            cond = 0;           
            for k1=1:num_labels
                for k2=k1+1:num_labels
                    cond = cond+1;
                    subset_train =  find ((train_label==k1) | (train_label==k2));
                    subset_test =  find ((test_label==k1) | (test_label==k2));
                    % train svm
                    svm_model_ts = compact(fitcsvm(train_data_ts(subset_train,:),double(train_label(subset_train))));                               
                    % fit posterior distribution 
                    model_ts = fitPosterior(svm_model_ts ,train_data_ts(subset_train,:),train_label(subset_train));                    
                    % get posterior test prob                    
                    pred_ts = predict(model_ts,test_data_ts(subset_test,:));                    
                    
                    acc_ts(s,cond) = sum(pred_ts == test_label(subset_test)) ./ numel(subset_test);                    
                    
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