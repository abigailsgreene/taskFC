% Code to split predictive features by their intersubject synchrony, and
% visualize high- and low-synchrony predictive features by network

% start wtih *pred_synch from visualize_interSubPPI (but call it
% signed_*pred_synch) and map of nodes to networks

signed_intrinsic_pred_synch = intrinsic_pred_synch;
signed_interaction_pred_synch = interaction_pred_synch;
signed_activation_pred_synch = activation_pred_synch;
clear intrinsic_pred_synch interaction_pred_synch activation_pred_synch

%take absolute value of predictive utility and synchrony, if desired
%(approach taken in Greene et al., 2019)
for t = 1:5
    intrinsic_pred_synch{t} = abs(signed_intrinsic_pred_synch{t});
    interaction_pred_synch{t} = abs(signed_interaction_pred_synch{t});
    activation_pred_synch{t} = abs(signed_activation_pred_synch{t});
end

% median split features into high synch and low synch
% get indices of edges >= median for each task/term
intrinsic_synch_bin = zeros(size(intrinsic_pred_synch{1,2},1),5);
interaction_synch_bin = zeros(size(interaction_pred_synch{1,2},1),5);
for i = 1:5
    intrinsic_synch_idx{i} = find(intrinsic_pred_synch{i}(:,2)>=median(intrinsic_pred_synch{i}(find(intrinsic_pred_synch{i}(:,1)),2)));
    interaction_synch_idx{i} = find(interaction_pred_synch{i}(:,2)>=median(interaction_pred_synch{i}(find(interaction_pred_synch{i}(:,1)),2)));
    intrinsic_synch_bin(intrinsic_synch_idx{i},i) = 1; % 1 = upper median, 0 = lower median synch
    interaction_synch_bin(interaction_synch_idx{i},i) = 1;
    intrinsic_synch_binMat(:,:,i) = squareform(intrinsic_synch_bin(:,i));
    interaction_synch_binMat(:,:,i) = squareform(interaction_synch_bin(:,i));
end
% set low synch to -1 instead of 0
intrinsic_synch_binMat(intrinsic_synch_binMat==0) = -1;
interaction_synch_binMat(interaction_synch_binMat==0) = -1;

% now sum up predictive contributions of all edges in given network pair
% for high synchrony and low synch features separately - taken
% directly from visualize_interSubPPI.m with signmat replaced
tasks = {'gambling' 'emotion' 'wm' 'social' 'relational'};
dataset = 'hcp';
homedir = '/data17/mri_group/abby_data/ppi_paper/';
pthresh = '0.1';
for tasknum = 1:length(tasks) % iterate through all tasks
    % first load in the selected edges from ridge prediction with this task,
    % and examine relationship between predictiveness and intersubject PPI
    ridgeResults = dir([homedir dataset '/results/' dataset '_p' pthresh '_100iters*' tasks{tasknum} '*ridge_visualization_noNormBetaContrib_withContribVec.mat'])
    load([ridgeResults(1).folder '/' ridgeResults(1).name]);
    term5_intercept = squeeze(edge_mat{5}(:,:,1));
    term5_intrinsic = squeeze(edge_mat{5}(:,:,2));
    term5_interaction = squeeze(edge_mat{5}(:,:,3));
    
    % threshold predictive contributions by binarized high/low synch (don't
    % really care about intercept)
    for circ = 1:3
        for pos = 1:2 % high synch edges first, then low synch
            if pos == 1 % only going to add up contribution of edges that have high synch
                if circ==1; tmpmat = term5_intercept.*(squeeze(intrinsic_synch_binMat(:,:,tasknum))==1);
                elseif circ==2; tmpmat = term5_intrinsic.*(squeeze(intrinsic_synch_binMat(:,:,tasknum))==1);
                elseif circ==3; tmpmat = term5_interaction.*(squeeze(interaction_synch_binMat(:,:,tasknum))==1);
                end
            else % only going to add up contribution of edges that have low synch
                if circ==1; tmpmat = term5_intercept.*(squeeze(intrinsic_synch_binMat(:,:,tasknum))==-1);
                elseif circ==2; tmpmat = term5_intrinsic.*(squeeze(intrinsic_synch_binMat(:,:,tasknum))==-1);
                elseif circ==3; tmpmat = term5_interaction.*(squeeze(interaction_synch_binMat(:,:,tasknum))==-1);
                end
            end
            for net1 = 1:10
                for net2 = 1:10
                    % so order of dim3 is intercept, intrinsic,
                    % interaction; dim4 is high/low synch
                    
                    tot_circmat_ridge_noNorm(net1,net2,circ,pos,tasknum) = sum(sum(tmpmat(table2array(map(map.category==net1,2)),table2array(map(map.category==net2,2)))));
                    if net1==net2
                        tot_circmat_ridge_noNorm(net1,net2,circ,pos,tasknum) = tot_circmat_ridge_noNorm(net1,net2,circ,pos,tasknum)/2; % to not overcount, since edges will be duplicated (don't have to worry about diagonal, since will never be selected)
                        tot_edges = (length(find(map.category==net1))*(length(find(map.category==net2))-1))/2; % here, don't include diagonal (i.e., node x - node x) since self edges aren't used for prediction
                    else
                        tot_edges = length(find(map.category==net1))*length(find(map.category==net2));
                    end
                    % normalize by number of edges in each network
                    tot_circmat_ridge(net1,net2,circ,pos,tasknum) = tot_circmat_ridge_noNorm(net1,net2,circ,pos,tasknum)/tot_edges;
                    clear tot_edges
                end
            end
            clear tmpmat
        end
    end
    clearvars -except *_pred_synch *binMat tasks dataset homedir pthresh tasknum map tot_circmat*
end

% now can visualize tot_circmat_ridge for different tasks, high/low synch,
% and different terms

% code to compare how similar high and low synch canon net matrices are,
% with theory that they're more similar for interaction than intrinsic
for i = 1:5
    mathigh = squeeze(tot_circmat_ridge_pct(:,:,2,1,i));
    matlow = squeeze(tot_circmat_ridge_pct(:,:,2,2,i));
    [r_intrinsic(i),p_intrinsic(i)] = corr(mathigh(lower_idx),matlow(lower_idx),'type','Spearman');
    clear mathigh matlow
end

for i = 1:5
    mathigh = squeeze(tot_circmat_ridge_pct(:,:,3,1,i));
    matlow = squeeze(tot_circmat_ridge_pct(:,:,3,2,i));
    [r_interaction(i),p_interaction(i)] = corr(mathigh(lower_idx),matlow(lower_idx),'type','Spearman');
    clear mathigh matlow
end