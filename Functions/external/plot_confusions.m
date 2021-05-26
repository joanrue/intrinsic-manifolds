dim = 4;
figure(1)
subplot(2,2,1)
mat = squeeze(mean(mean(SVM_results{dim}.confusion_TMH_ind,1),2));
myheatmap(mat)
title('TMH ind')
xticklabels({'A','N1','N2','N3'})
yticklabels({'A','N1','N2','N3'})
xlabel('True state')
ylabel('Predicted state')
subplot(2,2,2)
mat = squeeze(mean(mean(SVM_results{dim}.confusion_TMH_group,1),2));
myheatmap(mat)
title('TMH group')
xticklabels({'A','N1','N2','N3'})
yticklabels({'A','N1','N2','N3'})
xlabel('True state')
ylabel('Predicted state')
subplot(2,2,3)
mat = squeeze(mean(mean(SVM_results{dim}.confusion_PCA_ind,1),2));
myheatmap(mat)
title('PCA ind')
xticklabels({'A','N1','N2','N3'})
yticklabels({'A','N1','N2','N3'})
xlabel('True state')
ylabel('Predicted state')
subplot(2,2,4)
mat = squeeze(mean(mean(SVM_results{dim}.confusion_PCA_group,1),2));
myheatmap(mat)
title('PCA group')
xticklabels({'A','N1','N2','N3'})
yticklabels({'A','N1','N2','N3'})
xlabel('True state')
ylabel('Predicted state')



%%



figure(2)
subplot(2,2,1)
mat = squeeze(mean(mean(SVM_results_high_dim.confusion_ts_ind,1),2));
myheatmap(mat)
title('ts ind')
xticklabels({'A','N1','N2','N3'})
yticklabels({'A','N1','N2','N3'})
xlabel('True state')
ylabel('Predicted state')
subplot(2,2,2)
mat = squeeze(mean(mean(SVM_results_high_dim.confusion_ts_group,1),2));
myheatmap(mat)
title('ts group')
xticklabels({'A','N1','N2','N3'})
yticklabels({'A','N1','N2','N3'})
xlabel('True state')
ylabel('Predicted state')
subplot(2,2,3)
mat = squeeze(mean(mean(SVM_results_high_dim.confusion_phases_ind,1),2));
myheatmap(mat)
title('phases ind')
xticklabels({'A','N1','N2','N3'})
yticklabels({'A','N1','N2','N3'})
xlabel('True state')
ylabel('Predicted state')
subplot(2,2,4)
mat = squeeze(mean(mean(SVM_results_high_dim.confusion_phases_group,1),2));
myheatmap(mat)
title('phases group')
xticklabels({'A','N1','N2','N3'})
yticklabels({'A','N1','N2','N3'})
xlabel('True state')
ylabel('Predicted state')


