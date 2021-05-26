function [TMH_all, PCA_all] = ROC_AUC_TMH(ind_data,opts,stat_type)
% Area under the Curve of Receiver operating characteristic (AUC-ROC). 
%   INPUTS: 
%       - ind_data: cell(1,num_subjects)
%            ind_data{1,num_sub}: struct with at least the fields: 
%                    % TMH (time,dimensions)
%                    % PCA (time,dimensions)
%                    % Stages (time,1) -> stage of each TR
%       - group_data: struct with fields:
%                    % TMH (concatenated time for all subjects, dimensions)
%                    % PCA (concatenated time for all subjects, dimensions)
%                    % stages (1, concatenated time for all subjects)
%                    % subj_vector (1, concatenated time for all
%                    subjects) --> to which subject corresponds each TR.
%
%

%%
Nsub = length(ind_data);

TimeA = {};
Time1 = {};
Time2 = {};
Time3 = {};
Stage_vector = [];

for nsub =1:Nsub
    TimeA{nsub}=find(ind_data{1,nsub}.stages==0);
    Time1{nsub}=find(ind_data{1,nsub}.stages==-2);
    Time2{nsub}=find(ind_data{1,nsub}.stages==-3);
    Time3{nsub}=find(ind_data{1,nsub}.stages==-4);
    
    % Save the indices of the selected TRs in a vector
    inds = [ TimeA{nsub}; Time1{nsub}; Time2{nsub}; Time3{nsub}];
    Stage_vector = [Stage_vector,ind_data{1,nsub}.stages(inds)];
end

%%
best_dim = zeros(Nsub,6);
second_best_dim  = zeros(Nsub,6);
bestauc = zeros(Nsub,6);
best_dim_pca = zeros(Nsub,6);
second_best_dim_pca = zeros(Nsub,6);
bestauc_pca = zeros(Nsub,6);
for ndim = 1:opts.dimension
    for nsub = 1:Nsub
        
        Clase(1:length(TimeA{nsub}))=1;
        Clase(length(TimeA{nsub})+1:length(TimeA{nsub})+length(Time1{nsub}))=2;
        Htot=vertcat(ind_data{1,nsub}.TMH(Stage_vector(:,1)==0,ndim),ind_data{1,nsub}.TMH(Stage_vector(:,1)==-2,ndim));
        auc=colAUC(Htot,Clase);
        if bestauc(nsub,1)< abs(auc-0.5)+0.5
            TMH{nsub}.auc1=abs(auc-0.5)+0.5;
            bestauc(nsub,1) = abs(auc-0.5)+0.5;
            second_best_dim(nsub,1) = best_dim(nsub,1);
            best_dim(nsub,1) = ndim;
            
            colAUC(Htot,Clase);
            tmp = findobj(gca,'Type','line');
            
            [x, index] = unique(tmp.XData);
            y = tmp.YData(index);
            TMH{nsub}.line1.XData = linspace(0,1);
            TMH{nsub}.line1.YData = interp1(x,y,linspace(0,1));
        end
        
        Clase(1:length(TimeA{nsub}))=1;
        Clase(length(TimeA{nsub})+1:length(TimeA{nsub})+length(Time2{nsub}))=2;
        Htot=vertcat(ind_data{1,nsub}.TMH(Stage_vector(:,1)==0,ndim),ind_data{1,nsub}.TMH(Stage_vector(:,1)==-3,ndim));
        auc=colAUC(Htot,Clase);
        if bestauc(nsub,2)< abs(auc-0.5)+0.5
            TMH{nsub}.auc2=abs(auc-0.5)+0.5;
            bestauc(nsub,2) = abs(auc-0.5)+0.5;
            second_best_dim(nsub,2) = best_dim(nsub,2);
            best_dim(nsub,2) = ndim;
            colAUC(Htot,Clase);
            tmp = findobj(gca,'Type','line');
            [x, index] = unique(tmp.XData);
            y = tmp.YData(index);
            TMH{nsub}.line2.XData = linspace(0,1);
            TMH{nsub}.line2.YData = interp1(x,y,linspace(0,1));
        end
        
        Clase(1:length(TimeA{nsub}))=1;
        Clase(length(TimeA{nsub})+1:length(TimeA{nsub})+length(Time3{nsub}))=2;
        Htot=vertcat(ind_data{1,nsub}.TMH(Stage_vector(:,1)==0,ndim),ind_data{1,nsub}.TMH(Stage_vector(:,1)==-4,ndim));
        auc=colAUC(Htot,Clase);
        if bestauc(nsub,3)< abs(auc-0.5)+0.5
            TMH{nsub}.auc3=abs(auc-0.5)+0.5;
            bestauc(nsub,3) = abs(auc-0.5)+0.5;
            second_best_dim(nsub,3) = best_dim(nsub,3);
            best_dim(nsub,3) = ndim;
            colAUC(Htot,Clase);
            tmp = findobj(gca,'Type','line');
            
            [x, index] = unique(tmp.XData);
            y = tmp.YData(index);
            TMH{nsub}.line3.XData = linspace(0,1);
            TMH{nsub}.line3.YData = interp1(x,y,linspace(0,1));
        end
        Clase(1:length(Time1{nsub}))=1;
        Clase(length(Time1{nsub})+1:length(Time1{nsub})+length(Time2{nsub}))=2;
        Htot=vertcat(ind_data{1,nsub}.TMH(Stage_vector(:,1)==-2,ndim),ind_data{1,nsub}.TMH(Stage_vector(:,1)==-3,ndim));
        auc=colAUC(Htot,Clase);
        if bestauc(nsub,4)< abs(auc-0.5)+0.5
            TMH{nsub}.auc4=abs(auc-0.5)+0.5;
            bestauc(nsub,4) = abs(auc-0.5)+0.5;
            second_best_dim(nsub,4) = best_dim(nsub,4);
            best_dim(nsub,4) = ndim;
            colAUC(Htot,Clase);
            tmp = findobj(gca,'Type','line');
            
            [x, index] = unique(tmp.XData);
            y = tmp.YData(index);
            TMH{nsub}.line4.XData = linspace(0,1);
            TMH{nsub}.line4.YData = interp1(x,y,linspace(0,1));
        end
        
        Clase(1:length(Time1{nsub}))=1;
        Clase(length(Time1{nsub})+1:length(Time1{nsub})+length(Time3{nsub}))=2;
        Htot=vertcat(ind_data{1,nsub}.TMH(Stage_vector(:,1)==-2,ndim),ind_data{1,nsub}.TMH(Stage_vector(:,1)==-4,ndim));
        auc=colAUC(Htot,Clase);
        if bestauc(nsub,5)< abs(auc-0.5)+0.5
            TMH{nsub}.auc5=abs(auc-0.5)+0.5;
            bestauc(nsub,5) = abs(auc-0.5)+0.5;
            second_best_dim(nsub,5) = best_dim(nsub,5);
            best_dim(nsub,5) = ndim;
            colAUC(Htot,Clase);
            tmp = findobj(gca,'Type','line');
            
            [x, index] = unique(tmp.XData);
            y = tmp.YData(index);
            TMH{nsub}.line5.XData = linspace(0,1);
            TMH{nsub}.line5.YData = interp1(x,y,linspace(0,1));
        end
        
        Clase(1:length(Time2{nsub}))=1;
        Clase(length(Time2{nsub})+1:length(Time2{nsub})+length(Time3{nsub}))=2;
        Htot=vertcat(ind_data{1,nsub}.TMH(Stage_vector(:,1)==-3,ndim),ind_data{1,nsub}.TMH(Stage_vector(:,1)==-4,ndim));
        auc=colAUC(Htot,Clase);
        if bestauc(nsub,6)< abs(auc-0.5)+0.5
            TMH{nsub}.auc6=abs(auc-0.5)+0.5;
            bestauc(nsub,6) = abs(auc-0.5)+0.5;
            second_best_dim(nsub,6) = best_dim(nsub,6);
            best_dim(nsub,6) = ndim;
            colAUC(Htot,Clase);
            tmp = findobj(gca,'Type','line');
            
            [x, index] = unique(tmp.XData);
            y = tmp.YData(index);
            TMH{nsub}.line6.XData = linspace(0,1);
            TMH{nsub}.line6.YData = interp1(x,y,linspace(0,1));
        end
        
        
        Clase(1:length(TimeA{nsub}))=1;
        Clase(length(TimeA{nsub})+1:length(TimeA{nsub})+length(Time1{nsub}))=2;
        
        Htot=vertcat(ind_data{1,nsub}.PCA(Stage_vector(:,1)==0,1),ind_data{1,nsub}.PCA(Stage_vector(:,1)==-2,1));
        auc=colAUC(Htot,Clase);
        if bestauc_pca(nsub,1)< abs(auc-0.5)+0.5
            PCA{nsub}.auc1=abs(auc-0.5)+0.5;
            bestauc_pca(nsub,1) = abs(auc-0.5)+0.5;
            second_best_dim_pca(nsub,1) = best_dim_pca(nsub,1);
            best_dim_pca(nsub,1) = ndim;
            
            colAUC(Htot,Clase);
            tmp = findobj(gca,'Type','line');
            
            [x, index] = unique(tmp.XData);
            y = tmp.YData(index);
            PCA{nsub}.line1.XData = linspace(0,1);
            PCA{nsub}.line1.YData = interp1(x,y,linspace(0,1));
        end
        
        Clase(1:length(TimeA{nsub}))=1;
        Clase(length(TimeA{nsub})+1:length(TimeA{nsub})+length(Time2{nsub}))=2;
        
        Htot=vertcat(ind_data{1,nsub}.PCA(Stage_vector(:,1)==0,1),ind_data{1,nsub}.PCA(Stage_vector(:,1)==-3,1));
        auc=colAUC(Htot,Clase);
        if bestauc_pca(nsub,2)< abs(auc-0.5)+0.5
            PCA{nsub}.auc2=abs(auc-0.5)+0.5;
            bestauc_pca(nsub,2) = abs(auc-0.5)+0.5;
            second_best_dim_pca(nsub,2) = best_dim_pca(nsub,2);
            best_dim_pca(nsub,2) = ndim;
            
            colAUC(Htot,Clase);
            tmp = findobj(gca,'Type','line');
            
            [x, index] = unique(tmp.XData);
            y = tmp.YData(index);
            PCA{nsub}.line2.XData = linspace(0,1);
            PCA{nsub}.line2.YData = interp1(x,y,linspace(0,1));
        end
        
        Clase(1:length(TimeA{nsub}))=1;
        Clase(length(TimeA{nsub})+1:length(TimeA{nsub})+length(Time3{nsub}))=2;
        
        Htot=vertcat(ind_data{1,nsub}.PCA(Stage_vector(:,1)==0,1),ind_data{1,nsub}.PCA(Stage_vector(:,1)==-4,1));
        auc=colAUC(Htot,Clase);
        
        
        if bestauc_pca(nsub,3)< abs(auc-0.5)+0.5
            PCA{nsub}.auc3=abs(auc-0.5)+0.5;
            bestauc_pca(nsub,3) = abs(auc-0.5)+0.5;
            second_best_dim_pca(nsub,3) = best_dim_pca(nsub,3);
            best_dim_pca(nsub,3) = ndim;
            
            %     figure;
            colAUC(Htot,Clase);
            tmp = findobj(gca,'Type','line');
            
            [x, index] = unique(tmp.XData);
            y = tmp.YData(index);
            PCA{nsub}.line3.XData = linspace(0,1);
            PCA{nsub}.line3.YData = interp1(x,y,linspace(0,1));
        end
        
        Clase(1:length(Time1{nsub}))=1;
        Clase(length(Time1{nsub})+1:length(Time1{nsub})+length(Time2{nsub}))=2;
        
        Htot=vertcat(ind_data{1,nsub}.PCA(Stage_vector(:,1)==-2,1),ind_data{1,nsub}.PCA(Stage_vector(:,1)==-3,1));
        auc=colAUC(Htot,Clase);
        
        if bestauc_pca(nsub,4)< abs(auc-0.5)+0.5
            PCA{nsub}.auc4=abs(auc-0.5)+0.5;
            bestauc_pca(nsub,4) = abs(auc-0.5)+0.5;
            second_best_dim_pca(nsub,4) = best_dim_pca(nsub,4);
            best_dim_pca(nsub,4) = ndim;
            %     figure;
            colAUC(Htot,Clase);
            tmp = findobj(gca,'Type','line');
            
            [x, index] = unique(tmp.XData);
            y = tmp.YData(index);
            PCA{nsub}.line4.XData = linspace(0,1);
            PCA{nsub}.line4.YData = interp1(x,y,linspace(0,1));
        end
        
        Clase(1:length(Time1{nsub}))=1;
        Clase(length(Time1{nsub})+1:length(Time1{nsub})+length(Time3{nsub}))=2;
        
        Htot=vertcat(ind_data{1,nsub}.PCA(Stage_vector(:,1)==-2,1),ind_data{1,nsub}.PCA(Stage_vector(:,1)==-4,1));
        auc=colAUC(Htot,Clase);
        if bestauc_pca(nsub,5)< abs(auc-0.5)+0.5
            PCA{nsub}.auc5=abs(auc-0.5)+0.5;
            bestauc_pca(nsub,5) = abs(auc-0.5)+0.5;
            second_best_dim_pca(nsub,5) = best_dim_pca(nsub,5);
            best_dim_pca(nsub,5) = ndim;
            %     figure;
            colAUC(Htot,Clase);
            tmp = findobj(gca,'Type','line');
            
            [x, index] = unique(tmp.XData);
            y = tmp.YData(index);
            PCA{nsub}.line5.XData = linspace(0,1);
            PCA{nsub}.line5.YData = interp1(x,y,linspace(0,1));
        end
        
        Clase(1:length(Time2{nsub}))=1;
        Clase(length(Time2{nsub})+1:length(Time2{nsub})+length(Time3{nsub}))=2;
        
        Htot=vertcat(ind_data{1,nsub}.PCA(Stage_vector(:,1)==-3,1),ind_data{1,nsub}.PCA(Stage_vector(:,1)==-4,1));
        auc=colAUC(Htot,Clase);
        if bestauc_pca(nsub,6)< abs(auc-0.5)+0.5
            PCA{nsub}.auc6=abs(auc-0.5)+0.5;
            bestauc_pca(nsub,6) = abs(auc-0.5)+0.5;
            second_best_dim_pca(nsub,6) = best_dim_pca(nsub,6);
            best_dim_pca(nsub,6) = ndim;
            
            %     figure;
            colAUC(Htot,Clase);
            tmp = findobj(gca,'Type','line');
            
            [x, index] = unique(tmp.XData);
            y = tmp.YData(index);
            PCA{nsub}.line6.XData = linspace(0,1);
            PCA{nsub}.line6.YData = interp1(x,y,linspace(0,1));
        end
    end
end


TMH_all.bestdim = best_dim;
PCA_all.bestdim = best_dim_pca;


%%
close all;
TMH_all.auc1 = [];
TMH_all.auc2 = [];
TMH_all.auc3 = [];
TMH_all.auc4 = [];
TMH_all.auc5 = [];
TMH_all.auc6 = [];
PCA_all.auc1 = [];
PCA_all.auc2 = [];
PCA_all.auc3 = [];
PCA_all.auc4 = [];
PCA_all.auc5 = [];
PCA_all.auc6 = [];

TMH_all.line1.XData = [];
TMH_all.line2.XData = [];
TMH_all.line3.XData = [];
TMH_all.line4.XData = [];
TMH_all.line5.XData = [];
TMH_all.line6.XData = [];
PCA_all.line1.XData = [];
PCA_all.line2.XData = [];
PCA_all.line3.XData = [];
PCA_all.line4.XData = [];
PCA_all.line5.XData = [];
PCA_all.line6.XData = [];

TMH_all.line1.YData = [];
TMH_all.line2.YData = [];
TMH_all.line3.YData = [];
TMH_all.line4.YData = [];
TMH_all.line5.YData = [];
TMH_all.line6.YData = [];
PCA_all.line1.YData = [];
PCA_all.line2.YData = [];
PCA_all.line3.YData = [];
PCA_all.line4.YData = [];
PCA_all.line5.YData = [];
PCA_all.line6.YData = [];

for nsub = 1:Nsub
    TMH_all.auc1 = [TMH_all.auc1,TMH{nsub}.auc1];
    TMH_all.auc2 = [TMH_all.auc2,TMH{nsub}.auc2];
    TMH_all.auc3 = [TMH_all.auc3,TMH{nsub}.auc3];
    TMH_all.auc4 = [TMH_all.auc4,TMH{nsub}.auc4];
    TMH_all.auc5 = [TMH_all.auc5,TMH{nsub}.auc5];
    TMH_all.auc6 = [TMH_all.auc6,TMH{nsub}.auc6];
    PCA_all.auc1 = [PCA_all.auc1, PCA{nsub}.auc1];
    PCA_all.auc2 = [PCA_all.auc2, PCA{nsub}.auc2];
    PCA_all.auc3 = [PCA_all.auc3, PCA{nsub}.auc3];
    PCA_all.auc4 = [PCA_all.auc4, PCA{nsub}.auc4];
    PCA_all.auc5 = [PCA_all.auc5, PCA{nsub}.auc5];
    PCA_all.auc6 = [PCA_all.auc6, PCA{nsub}.auc6];
    
    TMH_all.line1.YData = [TMH_all.line1.YData;TMH{nsub}.line1.YData];
    TMH_all.line2.YData = [TMH_all.line2.YData;TMH{nsub}.line2.YData];
    TMH_all.line3.YData = [TMH_all.line3.YData;TMH{nsub}.line3.YData];
    TMH_all.line4.YData = [TMH_all.line4.YData;TMH{nsub}.line4.YData];
    TMH_all.line5.YData = [TMH_all.line5.YData;TMH{nsub}.line5.YData];
    TMH_all.line6.YData = [TMH_all.line6.YData;TMH{nsub}.line6.YData];
    PCA_all.line1.YData = [PCA_all.line1.YData; PCA{nsub}.line1.YData];
    PCA_all.line2.YData = [PCA_all.line2.YData; PCA{nsub}.line2.YData];
    PCA_all.line3.YData = [PCA_all.line3.YData; PCA{nsub}.line3.YData];
    PCA_all.line4.YData = [PCA_all.line4.YData; PCA{nsub}.line4.YData];
    PCA_all.line5.YData = [PCA_all.line5.YData; PCA{nsub}.line5.YData];
    PCA_all.line6.YData = [PCA_all.line6.YData; PCA{nsub}.line6.YData];
end

TMH_all.pvals(1) = ranksum(TMH_all.auc1,PCA_all.auc1);
TMH_all.pvals(2) = ranksum(TMH_all.auc2,PCA_all.auc2);
TMH_all.pvals(3) = ranksum(TMH_all.auc3,PCA_all.auc3);
TMH_all.pvals(4) = ranksum(TMH_all.auc4,PCA_all.auc4);
TMH_all.pvals(5) = ranksum(TMH_all.auc5,PCA_all.auc5);
TMH_all.pvals(6) = ranksum(TMH_all.auc6,PCA_all.auc6);

[~,~,~,p_fdr] = fdr_bh(TMH_all.pvals);
TMH_all.pvals = reshape(p_fdr,size(TMH_all.pvals));


if strcmp(stat_type,'median')

  
    TMH_all.line1.YData_err = reshape(abs(prctile(TMH_all.line1.YData,[95,5])-median(TMH_all.line1.YData))',[1,100,2]);
    TMH_all.line2.YData_err = reshape(abs(prctile(TMH_all.line2.YData,[95,5])-median(TMH_all.line2.YData))',[1,100,2]);
    TMH_all.line3.YData_err = reshape(abs(prctile(TMH_all.line3.YData,[95,5])-median(TMH_all.line3.YData))',[1,100,2]);
    TMH_all.line4.YData_err = reshape(abs(prctile(TMH_all.line4.YData,[95,5])-median(TMH_all.line4.YData))',[1,100,2]);
    TMH_all.line5.YData_err = reshape(abs(prctile(TMH_all.line5.YData,[95,5])-median(TMH_all.line5.YData))',[1,100,2]);
    TMH_all.line6.YData_err = reshape(abs(prctile(TMH_all.line6.YData,[95,5])-median(TMH_all.line6.YData))',[1,100,2]);
    PCA_all.line1.YData_err = reshape(abs(prctile(PCA_all.line1.YData,[95,5])-median(PCA_all.line1.YData))',[1,100,2]);
    PCA_all.line2.YData_err = reshape(abs(prctile(PCA_all.line2.YData,[95,5])-median(PCA_all.line2.YData))',[1,100,2]);
    PCA_all.line3.YData_err = reshape(abs(prctile(PCA_all.line3.YData,[95,5])-median(PCA_all.line3.YData))',[1,100,2]);
    PCA_all.line4.YData_err = reshape(abs(prctile(PCA_all.line4.YData,[95,5])-median(PCA_all.line4.YData))',[1,100,2]);
    PCA_all.line5.YData_err = reshape(abs(prctile(PCA_all.line5.YData,[95,5])-median(PCA_all.line5.YData))',[1,100,2]);
    PCA_all.line6.YData_err = reshape(abs(prctile(PCA_all.line6.YData,[95,5])-median(PCA_all.line6.YData))',[1,100,2]);


    TMH_all.line1.YData = median(TMH_all.line1.YData);
    TMH_all.line2.YData = median(TMH_all.line2.YData);
    TMH_all.line3.YData = median(TMH_all.line3.YData);
    TMH_all.line4.YData = median(TMH_all.line4.YData);
    TMH_all.line5.YData = median(TMH_all.line5.YData);
    TMH_all.line6.YData = median(TMH_all.line6.YData);

    PCA_all.line1.YData = median(PCA_all.line1.YData);
    PCA_all.line2.YData = median(PCA_all.line2.YData);
    PCA_all.line3.YData = median(PCA_all.line3.YData);
    PCA_all.line4.YData = median(PCA_all.line4.YData);
    PCA_all.line5.YData = median(PCA_all.line5.YData);
    PCA_all.line6.YData = median(PCA_all.line6.YData);

end


if strcmp(stat_type,'mean')
    TMH_all.line1.YData_err = std(TMH_all.line1.YData);
    TMH_all.line2.YData_err = std(TMH_all.line2.YData);
    TMH_all.line3.YData_err = std(TMH_all.line3.YData);
    TMH_all.line4.YData_err = std(TMH_all.line4.YData);
    TMH_all.line5.YData_err = std(TMH_all.line5.YData);
    TMH_all.line6.YData_err = std(TMH_all.line6.YData);
    PCA_all.line1.YData_err = std(PCA_all.line1.YData);
    PCA_all.line2.YData_err = std(PCA_all.line2.YData);
    PCA_all.line3.YData_err = std(PCA_all.line3.YData);
    PCA_all.line4.YData_err = std(PCA_all.line4.YData);
    PCA_all.line5.YData_err = std(PCA_all.line5.YData);
    PCA_all.line6.YData_err = std(PCA_all.line6.YData);


    TMH_all.line1.YData = mean(TMH_all.line1.YData);
    TMH_all.line2.YData = mean(TMH_all.line2.YData);
    TMH_all.line3.YData = mean(TMH_all.line3.YData);
    TMH_all.line4.YData = mean(TMH_all.line4.YData);
    TMH_all.line5.YData = mean(TMH_all.line5.YData);
    TMH_all.line6.YData = mean(TMH_all.line6.YData);

    PCA_all.line1.YData = mean(PCA_all.line1.YData);
    PCA_all.line2.YData = mean(PCA_all.line2.YData);
    PCA_all.line3.YData = mean(PCA_all.line3.YData);
    PCA_all.line4.YData = mean(PCA_all.line4.YData);
    PCA_all.line5.YData = mean(PCA_all.line5.YData);
    PCA_all.line6.YData = mean(PCA_all.line6.YData);

end

TMH_all.line1.XData = linspace(0,1);
TMH_all.line2.XData = linspace(0,1);
TMH_all.line3.XData = linspace(0,1);
TMH_all.line4.XData = linspace(0,1);
TMH_all.line5.XData = linspace(0,1);
TMH_all.line6.XData = linspace(0,1);
PCA_all.line1.XData = linspace(0,1);
PCA_all.line2.XData = linspace(0,1);
PCA_all.line3.XData = linspace(0,1);
PCA_all.line4.XData = linspace(0,1);
PCA_all.line5.XData = linspace(0,1);
PCA_all.line6.XData = linspace(0,1);

end
