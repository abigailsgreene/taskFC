% Code to formally relate predictive utility to intersubject consistency
% Written by Abby Greene, 2019

% have to run visualize_intersub first to get *_pred_synch and change those
% variable names to signed_*pred_synch - column 1 = ridge
% weights (ie predictive utility) and column 2 is intersub synch

% Stuff to change, as necessary
nets = 'activation'; % options = full (specific networks) vs bwwin (bw = 0, within net = 1); vs hemisphere (0 = within hemi, 1 = cross) vs distance vs activation vs reliability
taskname = {'gambling','emotion','wm','social','relational'};
sign = 'absolute'; % do we want to take the absolute value of predictive utility and synchrony? 'absolute' = yes, 'signed' = no
standardize = 1; % if standardize==1, zscore predictors and outcome before modeling
% load in node-to-network mapping, with hemisphere labels for each node
load ~/Desktop/map268_subnetwork.mat

if strcmp(sign,'absolute')
    for t = 1:5
        intrinsic_pred_synch{t} = abs(signed_intrinsic_pred_synch{t});
        interaction_pred_synch{t} = abs(signed_interaction_pred_synch{t});
        activation_pred_synch{t} = abs(signed_activation_pred_synch{t});
    end
else
    for t = 1:5
        intrinsic_pred_synch{t} = signed_intrinsic_pred_synch{t};
        interaction_pred_synch{t} = signed_interaction_pred_synch{t};
        activation_pred_synch{t} = signed_activation_pred_synch{t};
    end
end

% label edges by network
labelmat = reshape(find(ones(10,10)),10,10);
within_idx = diag(labelmat);
%between_idx = labelmat(find(tril(ones(10,10),-1)));
netmat = zeros(268,268);
for i = 1:268
    for j = 1:268
        netmat(i,j) = labelmat(table2array(map(map.oldroi==i,3)),table2array(map(map.oldroi==j,3)));
    end
end

lower_idx = find(tril(ones(268,268),-1));

% version where include all edges, regardless of whether they're in
% predictive model or not (vast majority aren't) - not used in Greene et al
% (2019)
% for i = 1:5 % five tasks with intersubj synch
%     betas_intrinsic(:,i) = regress(intrinsic_pred_synch{i}(:,1),[ones(size(intrinsic_pred_synch{i},1),1) intrinsic_pred_synch{i}(:,2) netmat(lower_idx) intrinsic_pred_synch{i}(:,2).*netmat(lower_idx)]);
%     betas_interaction(:,i) = regress(interaction_pred_synch{i}(:,1),[ones(size(interaction_pred_synch{i},1),1) interaction_pred_synch{i}(:,2) netmat(lower_idx) interaction_pred_synch{i}(:,2).*netmat(lower_idx)]);
%     betas_activation(:,i) = regress(activation_pred_synch{i}(:,1),[ones(size(activation_pred_synch{i},1),1) activation_pred_synch{i}(:,2) diag(netmat) activation_pred_synch{i}(:,2).*diag(netmat)]);   
% end

%% version where we only use predictive edges
if strcmp(nets,'full') % use all 55 network pair labels
    netmat_lower = netmat(lower_idx);
    netmat_diag = diag(netmat);
elseif strcmp(nets,'bwwin') % only label edges as between (0) or within (1) network
    netmat_diag = diag(netmat);
    netmat_lower_tmp = netmat(lower_idx);
    netmat_lower = zeros(length(netmat_lower_tmp),1); % between network elements will remain 0
    netmat_lower(find(ismember(netmat_lower_tmp,within_idx))) = 1; % set within network elements to 1
elseif strcmp(nets,'hemisphere') % is edge across or within hemispheres?
    for i = 1:268
        if strcmp(map.hemisphere(find(map.oldroi==i)),'L')
            hemisphere(i) = 1; % R = 0, L = 1
        else hemisphere(i) = 0;
        end
    end
    interhemi = zeros(268,268);
    for i = 1:268
        for j = 1:268
            if hemisphere(i)~=hemisphere(j)
                interhemi(i,j) = 1; % 0 = same hemi, 1 = different
            else
            end
        end
    end
    netmat_lower = interhemi(lower_idx);
elseif strcmp(nets,'distance')
    % label edges based on their nodes' distance from each other (load in
    % file in which Euclidean distance is calculated between every node
    % pair)
    load /data1/software/node_dist.mat
    netmat_lower = node_dist(lower_idx);
elseif strcmp(nets,'activation') % use mean abs activation of nodes on either end of an edge as proxy for how much task perturbs that edge
    checkTaskAct % run code that calculates mean activation per node and loads in HCP effect size/node vectors (can simply load in relevant node x task group-level activation matrix)
elseif strcmp(nets,'reliability') % use reliability of edges as predictors
    % load in reliability for each edge, as calculated in Noble et al
    % (2017)
    load /data17/mri_group/abby_data/ppi_paper/reliabilityAnalyses/ICC_matrix_and_regions_nr_1_ns_1.mat
    netmat_lower = ICC_matrix(lower_idx);
end
for i = 1:5 % iterate through five tasks
    % use mean activation of nodes incident to edge as predictor
    if strcmp(nets,'activation')
        clear netmat_lower
        actmat = zeros(268,268);
        netmat_diag = zeros(268,1);
        for j = 1:268
            for k = 1:268
                if j==k
                    netmat_diag(j) = abs(task_act_GT(j,i));
                else
                actmat(j,k) = (abs(task_act_GT(j,i))+abs(task_act_GT(k,i)))/2; % mean activation of two nodes that edge j,k connects
                end
            end
        end
        netmat_lower = actmat(lower_idx);
    end
    % get indices for predictive features
    intrinsic_idx = find(intrinsic_pred_synch{i}(:,1));
    interaction_idx = find(interaction_pred_synch{i}(:,1));
    activation_idx = find(activation_pred_synch{i}(:,1));
    
    if standardize==1 % nb added zscore around additional netmat term for distance, reliability, and activation (not dummy coded bw/win or hemisphere)
        [betas_intrinsic(:,i),~,~,~,stats_intrinsic{i}] = regress(zscore(intrinsic_pred_synch{i}(intrinsic_idx,1)),[ones(length(intrinsic_idx),1) zscore(intrinsic_pred_synch{i}(intrinsic_idx,2)) zscore(netmat_lower(intrinsic_idx)) (zscore(intrinsic_pred_synch{i}(intrinsic_idx,2).*netmat_lower(intrinsic_idx)))]);
        [betas_interaction(:,i),~,~,~,stats_interaction{i}] = regress(zscore(interaction_pred_synch{i}(interaction_idx,1)),[ones(length(interaction_idx),1) zscore(interaction_pred_synch{i}(interaction_idx,2)) zscore(netmat_lower(interaction_idx)) (zscore(interaction_pred_synch{i}(interaction_idx,2).*netmat_lower(interaction_idx)))]);
        if find(strcmp(nets,{'full','activation'})) % only do regression on activation term predictiveness if we're doing network analysis (meaningless for distance [0] or hemisphere [all same]) or using mean activation as covariate
            [betas_activation(:,i),~,~,~,stats_activation{i}] = regress(zscore(activation_pred_synch{i}(activation_idx,1)),[ones(length(activation_idx),1) zscore(activation_pred_synch{i}(activation_idx,2)) zscore(netmat_diag(activation_idx)) (zscore(activation_pred_synch{i}(activation_idx,2).*netmat_diag(activation_idx)))]);
        end
    else
        [betas_intrinsic(:,i),~,~,~,stats_intrinsic{i}] = regress(intrinsic_pred_synch{i}(intrinsic_idx,1),[ones(length(intrinsic_idx),1) intrinsic_pred_synch{i}(intrinsic_idx,2) netmat_lower(intrinsic_idx) intrinsic_pred_synch{i}(intrinsic_idx,2).*netmat_lower(intrinsic_idx)]);
        [betas_interaction(:,i),~,~,~,stats_interaction{i}] = regress(interaction_pred_synch{i}(interaction_idx,1),[ones(length(interaction_idx),1) interaction_pred_synch{i}(interaction_idx,2) netmat_lower(interaction_idx) interaction_pred_synch{i}(interaction_idx,2).*netmat_lower(interaction_idx)]);
        if find(strcmp(nets,{'full','activation'})) % only do regression on activation term predictiveness if we're doing network analysis (meaningless for distance [0] or hemisphere [all same]) or using mean activation as covariate
            [betas_activation(:,i),~,~,~,stats_activation{i}] = regress(activation_pred_synch{i}(activation_idx,1),[ones(length(activation_idx),1) activation_pred_synch{i}(activation_idx,2) netmat_diag(activation_idx) activation_pred_synch{i}(activation_idx,2).*netmat_diag(activation_idx)]);
        end
    end
end

%save(['/data17/mri_group/abby_data/ppi_paper/hcp/results/synch_vs_predictiveness/' nets '_modeling_zScorePredictors.mat']);

% plot resulting betas
figure; bar(categorical(taskname), betas_interaction','grouped'); title(['Predictive utility vs. synchrony, cdFC, ' nets, ', ' sign]); legend({'Intercept','Synchrony','Network membership','Synchrony x Network'})
figure; bar(categorical(taskname), betas_intrinsic','grouped'); title(['Predictive utility vs. synchrony, intrinsic FC, ' nets, ', ' sign]); legend({'Intercept','Synchrony','Network membership','Synchrony x Network'})
if find(strcmp(nets,{'full','activation'})) % only do regression on activation term predictiveness if we're doing network analysis (meaningless for distance [0] or hemisphere [all same]) or using mean activation as covariate
    figure; bar(categorical(taskname), betas_activation','grouped'); title(['Predictive utility vs. synchrony, activation, ' nets, ', ' sign]); legend({'Intercept','Synchrony','Network membership','Synchrony x Network'})
end

% also curious about correlation between synch and reliability
if strcmp(nets,'reliability')
    for i = 1:5
        [r_intrinsic(i),p_intrinsic(i)] = corr(intrinsic_pred_synch{i}(:,2),ICC_matrix(lower_idx),'type','Spearman');
        [r_interaction(i),p_interaction(i)] = corr(interaction_pred_synch{i}(:,2),ICC_matrix(lower_idx),'type','Spearman');
    end
end
