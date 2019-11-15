% Code to load in results of intersubject PPI and calculate group mean
% contrasts + prepare data for visualization in canonical network form
% Written by Abby Greene, 2019

tasks = {'gambling' 'emotion' 'wm' 'social' 'relational'};
dataset = 'hcp';
homedir = '/data17/mri_group/abby_data/ppi_paper/';
symmetrize = 1;
pthresh = '0.1';
load ~/Desktop/map268_subnetwork.mat % mapping of nodes to networks

% initialize variables to store predictiveness and synchrony for each
% task/term
intrinsic_pred_synch = cell(length(tasks),1);
interaction_pred_synch = cell(length(tasks),1);
activation_pred_synch = cell(length(tasks),1);
for tasknum = 1:length(tasks)
    % load intersubject PPI betas for given task
    load([homedir dataset '/results/betas/betas_' tasks{tasknum} '_zscore_withinSub_NodeTimecourses_3models_zeroCenteredTaskReg_convInterxn_condCondition_withCueReg_noGSR_interSubject_group1_all.mat'],'betas_1','betas_2');
    
    % calculate condition contrast betas
    subnum = size(betas_2,3);
    % regardless of task, 1st beta is always intercept and 2nd beta is
    % always timecourse/intrinsic FC
    all_mats(:,:,1) = (mean(squeeze(betas_2(:,:,:,1)),3)+mean(squeeze(betas_1(:,:,:,1)),3))./2; % intercept
    all_mats(:,:,2) = (mean(squeeze(betas_2(:,:,:,2)),3)+mean(squeeze(betas_1(:,:,:,2)),3))./2; % node time course
    
    % all_mats(:,:,:,3) = cdFC, all_vecs = node x sub activation
    % then average given term over all subjects, calculate the contrast
    % between the averages, and average the two runs' contrasts together
    if strcmp(tasks{tasknum},'wm') %wm
        all_vecs = mean((((mean(squeeze(betas_2(:,:,:,4)),3)-mean(squeeze(betas_2(:,:,:,5)),3))+(mean(squeeze(betas_1(:,:,:,4)),3)-mean(squeeze(betas_1(:,:,:,5)),3)))./2),2); % task contrast (2bk-0bk)
        all_mats(:,:,3) = ((mean(squeeze(betas_2(:,:,:,7)),3)-mean(squeeze(betas_2(:,:,:,8)),3))+(mean(squeeze(betas_1(:,:,:,7)),3)-mean(squeeze(betas_1(:,:,:,8)),3)))./2; %interaction contrast (2bk-0bk)
    elseif strcmp(tasks{tasknum},'language') % wrong but it's OK because not using language task for intersubject
        all_vecs = mean((((mean(squeeze(betas_2(:,:,:,3)),3))+(mean(squeeze(betas_1(:,:,:,3)),3)))./2),2); % task contrast (language vs. math [modeling only language because no rest])
        all_mats(:,:,3) = ((mean(squeeze(betas_2(:,:,:,4)),3))+(mean(squeeze(betas_1(:,:,:,4)),3)))./2; %interaction contrast (language vs. math)
    elseif strcmp(tasks{tasknum},'emotion')
        all_vecs = mean(((mean(squeeze(betas_2(:,:,:,4)),3)+mean(squeeze(betas_1(:,:,:,4)),3))./2),2); % just fear (vs neut = baseline)
        all_mats(:,:,3) = (mean(squeeze(betas_2(:,:,:,6)),3)+mean(squeeze(betas_1(:,:,:,6)),3))./2; % just fear interaction (vs neut = baseline)
    elseif strcmp(tasks{tasknum},'gambling') % gambling
        all_vecs = mean((((mean(squeeze(betas_2(:,:,:,4)),3)-mean(squeeze(betas_2(:,:,:,3)),3))+(mean(squeeze(betas_1(:,:,:,4)),3)-mean(squeeze(betas_1(:,:,:,3)),3)))./2),2); % task contrast (gam: win-loss)
        all_mats(:,:,3) = ((mean(squeeze(betas_2(:,:,:,6)),3)-mean(squeeze(betas_2(:,:,:,5)),3))+(mean(squeeze(betas_1(:,:,:,6)),3)-mean(squeeze(betas_1(:,:,:,5)),3)))./2; %interaction contrast (gam: win-loss)
    elseif strcmp(tasks{tasknum},'relational') % relational
        all_vecs = mean((((mean(squeeze(betas_2(:,:,:,5)),3)-mean(squeeze(betas_2(:,:,:,4)),3))+(mean(squeeze(betas_1(:,:,:,5)),3)-mean(squeeze(betas_1(:,:,:,4)),3)))./2),2); % task contrast (relational: rel-match)
        all_mats(:,:,3) = ((mean(squeeze(betas_2(:,:,:,8)),3)-mean(squeeze(betas_2(:,:,:,7)),3))+(mean(squeeze(betas_1(:,:,:,8)),3)-mean(squeeze(betas_1(:,:,:,7)),3)))./2; %interaction contrast (relational: rel-match)
    elseif strcmp(tasks{tasknum},'social') % social
        all_vecs = mean((((mean(squeeze(betas_2(:,:,:,3)),3)-mean(squeeze(betas_2(:,:,:,4)),3))+(mean(squeeze(betas_1(:,:,:,3)),3)-mean(squeeze(betas_1(:,:,:,4)),3)))./2),2); % task contrast (mental-rand)
        all_mats(:,:,3) = ((mean(squeeze(betas_2(:,:,:,5)),3)-mean(squeeze(betas_2(:,:,:,6)),3))+(mean(squeeze(betas_1(:,:,:,5)),3)-mean(squeeze(betas_1(:,:,:,6)),3)))./2; %interaction contrast (mental-rand)
    elseif strcmp(tasks{tasknum},'motor') % for motor, just moving (add all conditions) vs rest (unmodeled)
        all_vecs = mean((((mean(squeeze(betas_2(:,:,:,4)),3)+mean(squeeze(betas_2(:,:,:,5)),3)+mean(squeeze(betas_2(:,:,:,6)),3)+mean(squeeze(betas_2(:,:,:,7)),3)+mean(squeeze(betas_2(:,:,:,8)),3))+(mean(squeeze(betas_1(:,:,:,4)),3)+mean(squeeze(betas_1(:,:,:,5)),3)+mean(squeeze(betas_1(:,:,:,6)),3)+mean(squeeze(betas_1(:,:,:,7)),3)+mean(squeeze(betas_1(:,:,:,8)),3)))./2),2);
        all_mats(:,:,3) = ((mean(squeeze(betas_2(:,:,:,10)),3)+mean(squeeze(betas_2(:,:,:,11)),3)+mean(squeeze(betas_2(:,:,:,12)),3)+mean(squeeze(betas_2(:,:,:,13)),3)+mean(squeeze(betas_2(:,:,:,14)),3))+(mean(squeeze(betas_1(:,:,:,10)),3)+mean(squeeze(betas_1(:,:,:,11)),3)+mean(squeeze(betas_1(:,:,:,12)),3)+mean(squeeze(betas_1(:,:,:,13)),3)+mean(squeeze(betas_1(:,:,:,14)),3)))./2;
    elseif strcmp(tasks{tasknum},'mid') % win - loss
        all_vecs = mean((((mean(squeeze(betas_2(:,:,:,3)),3)-mean(squeeze(betas_2(:,:,:,4)),3))+(mean(squeeze(betas_1(:,:,:,3)),3)-mean(squeeze(betas_1(:,:,:,4)),3)))./2),2); % task
        all_mats(:,:,3) = ((mean(squeeze(betas_2(:,:,:,5)),3)-mean(squeeze(betas_2(:,:,:,6)),3))+(mean(squeeze(betas_1(:,:,:,5)),3)-mean(squeeze(betas_1(:,:,:,6)),3)))./2; % interaction
    end
    
    % symmetrize the beta matrices, except for task activation
    if symmetrize==1
        all_mats(:,:,1) = (squeeze(all_mats(:,:,1))+squeeze(all_mats(:,:,1))')./2;
        all_mats(:,:,2) = (squeeze(all_mats(:,:,2))+squeeze(all_mats(:,:,2))')./2;
        all_mats(:,:,3) = (squeeze(all_mats(:,:,3))+squeeze(all_mats(:,:,3))')./2;
    else
    end
    
    intercept = squeeze(all_mats(:,:,1));
    intrinsicFC = squeeze(all_mats(:,:,2));
    interactionFC = squeeze(all_mats(:,:,3));
    taskAct = diag(interactionFC);
    
    % create network-by-network matrix of isfc for  plotting
    load '/data17/mri_group/abby_data/ppi_paper/hcp/ppiBeta_sign.mat'; % mean (across subjects) signs of each feature; this is in original node order, so if reordering isFC on the fly in tmpmat, then must reorder signmat, too
    for circ = 1:3 % iterate through terms
        for pos = 1:2 % iterate through positive and negative features
            if pos==1
                if circ==1; tmpmat = intercept(table2array(map(:,2)),table2array(map(:,2))).*(squeeze(squeeze(signmat(table2array(map(:,2)),table2array(map(:,2)),2,tasknum)))==1);
                elseif circ==2; tmpmat = intrinsicFC(table2array(map(:,2)),table2array(map(:,2))).*(squeeze(squeeze(signmat(table2array(map(:,2)),table2array(map(:,2)),2,tasknum)))==1);
                elseif circ==3; tmpmat = interactionFC(table2array(map(:,2)),table2array(map(:,2))).*(squeeze(squeeze(signmat(table2array(map(:,2)),table2array(map(:,2)),2,tasknum)))==1);
                end
            else
                if circ==1; tmpmat = intercept(table2array(map(:,2)),table2array(map(:,2))).*(squeeze(squeeze(signmat(table2array(map(:,2)),table2array(map(:,2)),2,tasknum)))==-1);
                elseif circ==2; tmpmat = intrinsicFC(table2array(map(:,2)),table2array(map(:,2))).*(squeeze(squeeze(signmat(table2array(map(:,2)),table2array(map(:,2)),2,tasknum)))==-1);
                elseif circ==3; tmpmat = interactionFC(table2array(map(:,2)),table2array(map(:,2))).*(squeeze(squeeze(signmat(table2array(map(:,2)),table2array(map(:,2)),2,tasknum)))==-1);
                end
            end
            for net1 = 1:10
                for net2 = 1:10
                    if net1==net2 && symmetrize==1
                        tmpmat2 = tmpmat(table2array(map(map.category==net1,1)),table2array(map(map.category==net2,1)));
                        tot_circmat_intersub_noNorm(net1,net2,circ,pos,tasknum) = sum(tmpmat2(find(tril(ones(size(tmpmat2)))))); % to not overcount edges, since matrix is symm but we want to include the diagonal (ie self connections)
                        tot_edges = (length(find(map.category==net1))*(length(find(map.category==net2))+1))/2; % now self edges are meaningful, so can include them
                    else
                        tot_circmat_intersub_noNorm(net1,net2,circ,pos,tasknum) = sum(sum(tmpmat(table2array(map(map.category==net1,1)),table2array(map(map.category==net2,1)))));
                        tot_edges = length(find(map.category==net1))*length(find(map.category==net2));
                    end
                    % nwo normalize by number of edges between nets 1 and 2
                    tot_circmat_intersub(net1,net2,circ,pos,tasknum) = tot_circmat_intersub_noNorm(net1,net2,circ,pos,tasknum)/tot_edges;
                    clear tot_edges tmpmat2
                end
            end
            clear tmpmat*
        end
    end
    
    % now load in the selected edges from ridge prediction with this task,
    % (ie from ridge_visualization.m)
    % and examine relationship between predictiveness and intersubject PPI
    ridgeResults = dir([homedir dataset '/results/' dataset '_p' pthresh '_100iters*' tasks{tasknum} '*ridge_visualization_noNormBetaContrib_withContribVec.mat']);
    load([ridgeResults(1).folder '/' ridgeResults(1).name]);
   
    term5_intercept = squeeze(edge_mat{5}(:,:,1));
    term5_intrinsic = squeeze(edge_mat{5}(:,:,2));
    term5_interaction = squeeze(edge_mat{5}(:,:,3));
    term5_taskAct = diag(squeeze(edge_mat{5}(:,:,4)));
    
    % just select lower triangle for all terms
    lower_idx = find(tril(ones(268,268),-1));
    % rCPM results
    term5_intercept_lower = term5_intercept(lower_idx);
    term5_intrinsic_lower = term5_intrinsic(lower_idx);
    term5_interaction_lower = term5_interaction(lower_idx);
    % intersubj PPI results
    intercept_lower = intercept(lower_idx);
    intrinsicFC_lower = intrinsicFC(lower_idx);
    interactionFC_lower = interactionFC(lower_idx);
    % store predictiveness and synchrony for each term/task for modeling
    intrinsic_pred_synch{tasknum} = [term5_intrinsic_lower intrinsicFC_lower];
    interaction_pred_synch{tasknum} = [term5_interaction_lower interactionFC_lower];
    activation_pred_synch{tasknum} = [term5_taskAct taskAct];
end