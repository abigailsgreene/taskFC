function [] = ridge_wrapper_allTasks_permTest(symmetrize, all_behav, homedir, taskname, famid, dataset, savename, thresh, numiter, beta_idx)

% Inputs
% Symmetrize = 1 (yes) 0 (no)
% all_behav = behavior to predict (n x 1 vector)
% homedir = directory where data live
% taskname = cell array of tasks to run
% famid = 2-column matrix; col1 = subject ID; col2 = family ID; row order=
% same as all_behav
% beta_idx = optional (nonlogical) indexing variable if only using subset of
% subjects
% dataset = string indicating which dataset you're using (hcp or abcd)
% savename = string to include in output filename to indicate predicted
% measure
% thresh = pval threshold for edge selection (correlation with behavior)
% numiter = number of iterations


if ~exist('beta_idx','var')
    beta_idx = 1:length(all_behav);
else
end

% set true behavior and load in allowable permutations from FSL's PALM
% (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM/ExchangeabilityBlocks)
real_behav = all_behav; 
clear all_behav
load('/data17/mri_group/abby_data/private/palmQuickperms_output_ppiHCP_s1200_701subs.mat','Pset'); % Pset columns 2:end are allowable permutations

for tasknum = 1:length(taskname)
    % load betas for given task
    data_tmp = load([homedir dataset '/results/betas/betas_' taskname{tasknum} '_zscore_withinSub_NodeTimecourses_3models_zeroCenteredTaskReg_convInterxn_condCondition_withCueReg_noGSR_allBetas_withCue.mat'],'betas_1','betas_2');
    betas_1 = data_tmp.betas_1(:,:,beta_idx,:);
    betas_2 = data_tmp.betas_2(:,:,beta_idx,:);
    % calculate condition contrast betas
    subnum = size(betas_2,3);
    % regardless of task, 1st beta is always intercept and 2nd beta is
    % always timecourse/intrinsic FC
    for i = 1:subnum
        all_mats(:,:,i,1) = (squeeze(betas_2(:,:,i,1))+squeeze(betas_1(:,:,i,1)))./2; % intercept
        all_mats(:,:,i,2) = (squeeze(betas_2(:,:,i,2))+squeeze(betas_1(:,:,i,2)))./2; % node time course
    end
    % all_mats(:,:,:,3) = cdFC, all_vecs = node x sub activation
    if strcmp(taskname{tasknum},'wm') %wm
        for i = 1:subnum
            all_vecs(:,i) = mean((((squeeze(betas_2(:,:,i,4))-squeeze(betas_2(:,:,i,5)))+(squeeze(betas_1(:,:,i,4))-squeeze(betas_1(:,:,i,5))))./2),2); % task contrast (2bk-0bk)
            all_mats(:,:,i,3) = ((squeeze(betas_2(:,:,i,7))-squeeze(betas_2(:,:,i,8)))+(squeeze(betas_1(:,:,i,7))-squeeze(betas_1(:,:,i,8))))./2; %interaction contrast (2bk-0bk)
        end
    elseif strcmp(taskname{tasknum},'language')
        for i = 1:subnum
            all_vecs(:,i) = mean((((squeeze(betas_2(:,:,i,4)))+(squeeze(betas_1(:,:,i,4))))./2),2); % task contrast (language vs. math [modeling only language because no rest])
            all_mats(:,:,i,3) = ((squeeze(betas_2(:,:,i,6)))+(squeeze(betas_1(:,:,i,6))))./2; %interaction contrast (language vs. math)
        end
    elseif strcmp(taskname{tasknum},'emotion')
        for i = 1:subnum
           all_vecs(:,i) = mean(((squeeze(betas_2(:,:,i,4))+squeeze(betas_1(:,:,i,4)))./2),2); % just fear (vs neut = baseline)
           all_mats(:,:,i,3) = (squeeze(betas_2(:,:,i,6))+squeeze(betas_1(:,:,i,6)))./2; % just fear interaction (vs neut = baseline)
        end
    elseif strcmp(taskname{tasknum},'gambling') % gambling
        for i = 1:subnum
            all_vecs(:,i) = mean((((squeeze(betas_2(:,:,i,4))-squeeze(betas_2(:,:,i,3)))+(squeeze(betas_1(:,:,i,4))-squeeze(betas_1(:,:,i,3))))./2),2); % task contrast (gam: win-loss)
            all_mats(:,:,i,3) = ((squeeze(betas_2(:,:,i,6))-squeeze(betas_2(:,:,i,5)))+(squeeze(betas_1(:,:,i,6))-squeeze(betas_1(:,:,i,5))))./2; %interaction contrast (gam: win-loss)
        end
    elseif strcmp(taskname{tasknum},'relational') % relational
        for i = 1:subnum
            all_vecs(:,i) = mean((((squeeze(betas_2(:,:,i,5))-squeeze(betas_2(:,:,i,4)))+(squeeze(betas_1(:,:,i,5))-squeeze(betas_1(:,:,i,4))))./2),2); % task contrast (relational: rel-match)
            all_mats(:,:,i,3) = ((squeeze(betas_2(:,:,i,8))-squeeze(betas_2(:,:,i,7)))+(squeeze(betas_1(:,:,i,8))-squeeze(betas_1(:,:,i,7))))./2; %interaction contrast (relational: rel-match)
        end
    elseif strcmp(taskname{tasknum},'social') % social
        for i = 1:subnum
            all_vecs(:,i) = mean((((squeeze(betas_2(:,:,i,3))-squeeze(betas_2(:,:,i,4)))+(squeeze(betas_1(:,:,i,3))-squeeze(betas_1(:,:,i,4))))./2),2); % task contrast (mental-rand)
            all_mats(:,:,i,3) = ((squeeze(betas_2(:,:,i,5))-squeeze(betas_2(:,:,i,6)))+(squeeze(betas_1(:,:,i,5))-squeeze(betas_1(:,:,i,6))))./2; %interaction contrast (mental-rand)
        end
    elseif strcmp(taskname{tasknum},'motor') % for motor, just moving (add all conditions) vs rest (unmodeled)
        for i = 1:subnum
           all_vecs(:,i) = mean((((squeeze(betas_2(:,:,i,4))+squeeze(betas_2(:,:,i,5))+squeeze(betas_2(:,:,i,6))+squeeze(betas_2(:,:,i,7))+squeeze(betas_2(:,:,i,8)))+(squeeze(betas_1(:,:,i,4))+squeeze(betas_1(:,:,i,5))+squeeze(betas_1(:,:,i,6))+squeeze(betas_1(:,:,i,7))+squeeze(betas_1(:,:,i,8))))./2),2);
           all_mats(:,:,i,3) = ((squeeze(betas_2(:,:,i,10))+squeeze(betas_2(:,:,i,11))+squeeze(betas_2(:,:,i,12))+squeeze(betas_2(:,:,i,13))+squeeze(betas_2(:,:,i,14)))+(squeeze(betas_1(:,:,i,10))+squeeze(betas_1(:,:,i,11))+squeeze(betas_1(:,:,i,12))+squeeze(betas_1(:,:,i,13))+squeeze(betas_1(:,:,i,14))))./2;
        end
    elseif strcmp(taskname{tasknum},'mid') % win - loss
        for i = 1:subnum
           all_vecs(:,i) = mean((((squeeze(betas_2(:,:,i,3))-squeeze(betas_2(:,:,i,4)))+(squeeze(betas_1(:,:,i,3))-squeeze(betas_1(:,:,i,4))))./2),2); % task
           all_mats(:,:,i,3) = ((squeeze(betas_2(:,:,i,5))-squeeze(betas_2(:,:,i,6)))+(squeeze(betas_1(:,:,i,5))-squeeze(betas_1(:,:,i,6))))./2; % interaction
        end
    end
    
    % symmetrize the beta matrices, except for task activation
    if symmetrize==1
       for i = 1:subnum
           all_mats(:,:,i,1) = (squeeze(all_mats(:,:,i,1))+squeeze(all_mats(:,:,i,1))')./2;
           all_mats(:,:,i,2) = (squeeze(all_mats(:,:,i,2))+squeeze(all_mats(:,:,i,2))')./2;
           all_mats(:,:,i,3) = (squeeze(all_mats(:,:,i,3))+squeeze(all_mats(:,:,i,3))')./2;
       end
    else
    end
    
    beta_order = {'intercept','intrinsic','interaction','task','all'};
    % initialize things for rCPM
    k = 10;
    seeds = randi(1000,1,numiter);
    q_s = zeros(numiter,5);
    r_pearson = zeros(numiter,5);
    r_rank = zeros(numiter,5);
    y = cell(5,numiter);
    coef_total = cell(5,numiter);
    coef0_total = cell(5,numiter);
    lambda_total = cell(5,numiter);
    for iter = 1:numiter
        all_behav = real_behav(Pset(:,iter+1))';
        disp(['size of all_behav is ' num2str(size(all_behav))]);
        for i = 1:5
            try
            if i<=3 % first include each beta individually
                all_mats_tot = squeeze(all_mats(:,:,:,i));
                [q_s(iter,i), r_pearson(iter,i), r_rank(iter,i), y{i,iter}, coef_total{i,iter}, coef0_total{i,iter}, lambda_total{i,iter}] = ridgeCPM_vector_new(all_mats_tot, [], all_behav, thresh, famid, [], [], k, seeds(iter));
            elseif i==4 % now task activation alone
                all_vecs_tot = all_vecs;
                [q_s(iter,i), r_pearson(iter,i), r_rank(iter,i), y{i,iter}, coef_total{i,iter}, coef0_total{i,iter}, lambda_total{i,iter}] = ridgeCPM_vector_new([], all_vecs_tot, all_behav, thresh, famid, [], [], k, seeds(iter));
            elseif i==5 % on last run-through, include all terms
                all_mats_tot = all_mats;
                all_vecs_tot = all_vecs;
                [q_s(iter,i), r_pearson(iter,i), r_rank(iter,i), y{i,iter}, coef_total{i,iter}, coef0_total{i,iter}, lambda_total{i,iter}] = ridgeCPM_vector_new(all_mats_tot, all_vecs_tot, all_behav, thresh, famid, [], [], k, seeds(iter));
            end
            catch % if term fails (ie no edges correlated with gF), move on to next term
                disp(['no edges correlate with gF for term ' beta_order{i} ', task ' taskname{tasknum} '. moving on to next term'])
                disp(['seed = ' num2str(seeds(iter))]);
            end
        end
        all_behav = [];
    end
    save([homedir dataset '/results/',dataset, '_', taskname{tasknum} '_ridge_p',num2str(thresh),'_' num2str(numiter) 'iters_' savename '_prediction_noGSR_permTest_withCue.mat'],'q_s','r_pearson', 'r_rank','y','coef_total','coef0_total','lambda_total','beta_order','seeds','beta_idx','-v7.3');
    clearvars -except homedir taskname tasknum dataset symmetrize all_behav famid beta_idx numiter thresh savename Pset real_behav

end

end
