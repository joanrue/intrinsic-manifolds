function [group_data,ind_data] = compute_group_manifold(group_data,ind_data,opts,Nsub,tss)

group_data.TMH = zeros(tss*Nsub,opts.dimension);
ori_data = ind_data{1}.TMH(:,1:opts.dimension);
new_data = zeros(Nsub,size(ori_data,1),opts.dimension);
new_data(1,:,:) = ori_data;
new_data = new_data-mean(new_data);
ori_stages = ind_data{1}.stages;
stages = unique(ori_stages);

Y = zeros(length(stages)*2,opts.dimension);
for j = 1:length(stages)
    stage_data = ori_data(ori_stages==stages(j),:);

    norms = sqrt(sum((stage_data.^2),2));
    [~,maxind] = max(norms);
    [~,minind] = min(norms);
    Y(j*2-1,:) = stage_data(maxind,:);
    Y(j*2,:) = stage_data(minind,:);
end

for i = 2:Nsub

    sub_data = ind_data{i}.TMH(:,1:opts.dimension);
    sub_data = sub_data-mean(sub_data);
    sub_stages = ind_data{i}.stages;
    X = zeros(length(stages)*2,opts.dimension);
    for j = 1:length(stages)
        stage_data = sub_data(sub_stages==stages(j),:);       
        norms = sqrt(sum((stage_data.^2),2));
        [~,maxind] = max(norms);
        [~,minind] = min(norms);
        X(j*2-1,:) = stage_data(maxind,:);
        X(j*2,:) = stage_data(minind,:);
    end    
    
    [~,~,transform] = procrustes(Y,X);
    new_data(i,:,:) = transform.b.*sub_data*transform.T+repmat(transform.c(1,:),[size(sub_data,1),1]);
end

for i = 1:Nsub
    ind_data{i}.TMH = squeeze(new_data(i,:,:));
    group_data.TMH((i-1)*tss*length(stages)+1:i*tss*length(stages),:) = squeeze(new_data(i,:,:));
end
end