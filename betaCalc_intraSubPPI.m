% Perform PPI approach adapted from Cole et al. (2013) on HCP data
% Written by Abby Greene, Jan-March 2019
function [] = betaCalc_intraSubPPI(dataset,cond,cue_model,tr,homedir,tasks,taskdirnames)

% input arguments
% 0. dataset = 'hcp' or 'abcd'
% 1. cond = 'all' (just task on/off) or 'condition' (e.g., 2bk vs. 0bk)
% 2. cue_model = include cue regressor as predictor? 0 = no, 1 = yes
% 3. tr = TR in seconds (eg HCP = 0.72)
% 4. homedir = directory where regressors and roimean data are found, and
% where data should be saved
% 5. tasks = cell array list of tasks to perform (options = emotion, gambling,
% motor, relational, social, wm, and language)
% 6. taskdirnames = cell array list of names of task directories (e.g., for
% HCP, tasks = in all caps) in same order as tasks; if doesn't exist, default = tasks
% Note: for HCP, LR = run 2; RL = run 1;

% check if taskdirnames exists
if ~exist('taskdirnames','var')
    taskdirnames = tasks;
else
end

% add SPM path
addpath(genpath('/data1/software/spm8'));

% Generate appropriate HRF (corresponding to used TR; HCP = 0.72)
[hrf,~] = spm_hrf(tr);

% Load subject IDs
load([homedir dataset '/' dataset '_subids.mat']);

% Load IDs of missing nodes (if they exist/we're using them - for now (march 2019), just excluding subs with missing nodes so these files don't exist)
if exist([homedir dataset '/' dataset '_missingNodes.mat'])
    load([homedir dataset '/' dataset '_missingNodes.mat'])
else
end

for t = 1:length(tasks)
    
    % Load boxcar task regressor (to start, just total reg: task on vs. off)
    if cue_model == 1 & ~isempty(find(strcmp(tasks{t},{'emotion','motor','relational','wm', 'language'}))) & strcmp(dataset,'hcp')
        cue_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_cue_lr_regressors.mat']));
        cue_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_cue_rl_regressors.mat']));
    elseif cue_model == 1 & strcmp(tasks{t},'wm') & strcmp(dataset,'abcd')
        load([homedir dataset '/regressors/nback_regressors.mat'],'nback_cue_regressor_ds');
        for i = 1:length(subid)
            cue_2{1,1}{i} = nback_cue_regressor_ds{i,2};
            cue_1{1,1}{i} = nback_cue_regressor_ds{i,1};
        end
    else
    end
    
    if strcmp(cond,'all')
        if strcmp(dataset,'hcp') % this won't really work for emotion or language, since all conditions will be modeled...
            tmp_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_lr_regressor_ds'));
            tmp_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_rl_regressor_ds'));
        elseif strcmp(dataset,'abcd')
            if strcmp(tasks{t},'wm') % 2bk and 0bk modeled together, fixation left unmodeled
                load([homedir dataset '/regressors/nback_regressors.mat'],'nback_regressor_ds');
                for i = 1:length(subid)
                    tmp_2{i} = nback_regressor_ds{i,2};
                    tmp_1{i} = nback_regressor_ds{i,1};
                end
            elseif strcmp(tasks{t},'mid') % win and loss modeled together; neut left unmodeled
                load([homedir dataset '/regressors/mid_regressors.mat'],'mid_regressor_ds');
                for i = 1:length(subid)
                    tmp_2{i} = mid_regressor_ds{i,2};
                    tmp_1{i} = mid_regressor_ds{i,1};
                end
            end
        end
    elseif strcmp(cond,'condition') % 1st number after tmp = run#, 2nd number after tmp = condition
        if strcmp(tasks{t},'emotion') % 1 = fear, 2 = neut
            tmp_2_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*fear*'));
            tmp_1_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*fear*'));
        elseif strcmp(tasks{t},'gambling')
            tmp_2_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*loss*'));
            tmp_2_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*win*'));
            tmp_1_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*loss*'));
            tmp_1_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*win*'));
        elseif strcmp(tasks{t},'motor')
            tmp_2_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_lf_*'));
            tmp_2_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_lh_*'));
            tmp_2_3 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_rf_*'));
            tmp_2_4 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_rh_*'));
            tmp_2_5 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_t_*'));
            tmp_1_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_lf_*'));
            tmp_1_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_lh_*'));
            tmp_1_3 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_rf_*'));
            tmp_1_4 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_rh_*'));
            tmp_1_5 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_t_*'));
        elseif strcmp(tasks{t},'relational')
            tmp_2_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_match_*'));
            tmp_2_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_relation_*'));
            tmp_1_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_match_*'));
            tmp_1_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_relation_*'));
        elseif strcmp(tasks{t},'social')
            tmp_2_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_mental_*'));
            tmp_2_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_rnd_*'));
            tmp_1_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_mental_*'));
            tmp_1_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_rnd_*'));
        elseif strcmp(tasks{t},'wm')
            if strcmp(dataset,'hcp')
                tmp_2_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_2bk_*'));
                tmp_2_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_0bk_*'));
                tmp_1_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_2bk_*'));
                tmp_1_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_0bk_*'));
            elseif strcmp(dataset,'abcd')
                load([homedir dataset '/regressors/nback_regressors.mat'],'nback_2bk_regressor_ds','nback_0bk_regressor_ds');
                for i = 1:length(subid)
                    tmp_2_1{1,1}{i} = nback_2bk_regressor_ds{i,2};
                    tmp_2_2{1,1}{i} = nback_0bk_regressor_ds{i,2};
                    tmp_1_1{1,1}{i} = nback_2bk_regressor_ds{i,1};
                    tmp_1_2{1,1}{i} = nback_0bk_regressor_ds{i,1};
                end
            end
        elseif strcmp(tasks{t},'mid')
            load([homedir dataset '/regressors/mid_regressors.mat'],'mid_win_regressor_ds','mid_loss_regressor_ds');
            for i = 1:length(subid)
                tmp_2_1{1,1}{i} = mid_win_regressor_ds{i,2};
                tmp_2_2{1,1}{i} = mid_loss_regressor_ds{i,2};
                tmp_1_1{1,1}{i} = mid_win_regressor_ds{i,1};
                tmp_1_2{1,1}{i} = mid_loss_regressor_ds{i,1};
            end
        elseif strcmp(tasks{t},'language')
            tmp_2_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_lang_*'));
            tmp_2_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_LR_regressors.mat'],'*_math_*'));
            tmp_1_1 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_lang_*'));
            tmp_1_2 = struct2cell(load([homedir dataset '/regressors/' tasks{t} '_RL_regressors.mat'],'*_math_*'));
        end
    end
    % Now concatenate regressors for given task and subject
    for sub = 1:length(subid) 
        if strcmp(cond,'all')
            if ~isempty(find(strcmp({'gambling','social','language','mid'},tasks{t}))) | cue_model==0 % tasks without cues or not modeling cue
                task_reg_2 = tmp_2{1,1}{sub}'; % Regressors should be column vectors
                task_reg_1 = tmp_1{1,1}{sub}';
            else
                task_reg_2 = [cue_2{1,1}{sub}' tmp_2{1,1}{sub}'];
                task_reg_1 = [cue_1{1,1}{sub}' tmp_1{1,1}{sub}'];
            end
        elseif strcmp(cond,'condition') && strcmp(tasks{t},'motor') % If more than one condition, each column is regressor for one condition (task 3/motor is only one with >2 conditions), yielding time x condition predictor matrix
            if cue_model==0
                task_reg_2 = [tmp_2_1{1,1}{sub}' tmp_2_2{1,1}{sub}' tmp_2_3{1,1}{sub}' tmp_2_4{1,1}{sub}' tmp_2_5{1,1}{sub}'];
                task_reg_1 = [tmp_1_1{1,1}{sub}' tmp_1_2{1,1}{sub}' tmp_1_3{1,1}{sub}' tmp_1_4{1,1}{sub}' tmp_1_5{1,1}{sub}'];
            elseif cue_model==1
                task_reg_2 = [cue_2{1,1}{sub}' tmp_2_1{1,1}{sub}' tmp_2_2{1,1}{sub}' tmp_2_3{1,1}{sub}' tmp_2_4{1,1}{sub}' tmp_2_5{1,1}{sub}'];
                task_reg_1 = [cue_1{1,1}{sub}' tmp_1_1{1,1}{sub}' tmp_1_2{1,1}{sub}' tmp_1_3{1,1}{sub}' tmp_1_4{1,1}{sub}' tmp_1_5{1,1}{sub}'];
            end
        elseif strcmp(cond,'condition') && strcmp(tasks{t},'emotion') % no fixation in emotion so only model one condition
            if cue_model==0
                task_reg_2 = [tmp_2_1{1,1}{sub}'];
                task_reg_1 = [tmp_1_1{1,1}{sub}'];
            else
                task_reg_2 = [cue_2{1,1}{sub}' tmp_2_1{1,1}{sub}'];
                task_reg_1 = [cue_1{1,1}{sub}' tmp_1_1{1,1}{sub}'];
            end
        elseif strcmp(cond,'condition') && strcmp(tasks{t},'mid') % no fixation in MID so only model 2/3 conditions (win and loss) + no cue
            task_reg_2 = [tmp_2_1{1,1}{sub}' tmp_2_2{1,1}{sub}'];
            task_reg_1 = [tmp_1_1{1,1}{sub}' tmp_1_2{1,1}{sub}'];
        elseif strcmp(cond,'condition') && strcmp(tasks{t},'language') % corrected on 9/2/19 - yes cue in language, no fixation - only model one condition
            task_reg_2 = [cue_2{1,1}{sub}' tmp_2_1{1,1}{sub}'];
            task_reg_1 = [cue_1{1,1}{sub}' tmp_1_1{1,1}{sub}'];
        else %it's cond = condition and one of four remaining tasks (all of which have fixation)
            if ~isempty(find(strcmp({'gambling','social'},tasks{t}))) | cue_model==0 %  if it's a task without a cue but with fixation, or we're not modeling cue
                task_reg_2 = [tmp_2_1{1,1}{sub}' tmp_2_2{1,1}{sub}'];
                task_reg_1 = [tmp_1_1{1,1}{sub}' tmp_1_2{1,1}{sub}'];
            else % if it's cond=condition and task = WM or relational (have cue and fixation, and only 2 conditions)
                task_reg_2 = [cue_2{1,1}{sub}' tmp_2_1{1,1}{sub}' tmp_2_2{1,1}{sub}'];
                task_reg_1 = [cue_1{1,1}{sub}' tmp_1_1{1,1}{sub}' tmp_1_2{1,1}{sub}'];
            end
        end
        
        % Binarize the convolved task regressor for the interaction term (per Cole
        % et al) (N.B. convolving non zero-centered regressor, then zero centering as we binarize)
        reg_conv_bin_2 = conv2(task_reg_2,hrf);
        reg_conv_bin_2(size(task_reg_2,1)+1:end,:) = [];
        reg_conv_bin_2(find(reg_conv_bin_2>0)) = 0.5; 
        reg_conv_bin_2(find(reg_conv_bin_2<=0)) = -0.5;
        
        reg_conv_bin_1 = conv2(task_reg_1,hrf);
        reg_conv_bin_1(size(task_reg_1,1)+1:end,:) = [];
        reg_conv_bin_1(find(reg_conv_bin_1>0)) = 0.5;
        reg_conv_bin_1(find(reg_conv_bin_1<=0)) = -0.5;
        
        %         % Alternative to binarizing convolved task regressor for
        %         % interaction term: shift boxcar by 6 seconds (ie append 8 TRs'
        %         % worth of -0.5s to top of each regressor) and binarize regressor
        %         % (subtract off 0.5 from each)
        %         reg_conv_bin_2 = [repmat(-0.5,8,size(task_reg_lr,2)); task_reg_lr-0.5];
        %         reg_conv_bin_1 = [repmat(-0.5,8,size(task_reg_rl,2)); task_reg_rl-0.5];
        %         reg_conv_bin_2(size(task_reg_lr,1)+1:end,:) = []; % chop off those appended time points on the other end
        %         reg_conv_bin_1(size(task_reg_rl,1)+1:end,:) = [];
        
        %Convolve the task regressor with the HRF, then zero center it
        reg_conv_2 = conv2(task_reg_2,hrf)-0.5; % convolve SPM's canonical HRF with HCP task regressor
        reg_conv_1 = conv2(task_reg_1,hrf)-0.5;
        
        reg_conv_2(size(task_reg_2,1)+1:end,:) = []; % get rid of any time points after last TR
        reg_conv_1(size(task_reg_1,1)+1:end,:) = []; % get rid of any time points after last TR
        
        % Load in roimean data (i.e., time x node matrix)
        disp(['running subject ' num2str(sub) ' out of ' num2str(length(subid))]);
        if strcmp(dataset,'hcp')
            roi_data2_tmp = importdata([homedir dataset '/hcp1200/' num2str(subid(sub)) '/' taskdirnames{t} '/LR/matrices/' num2str(subid(sub)) '_' taskdirnames{t} '_LR_NOGSR_roimean.txt']);
            roi_data2_tmp2 = roi_data2_tmp.data(:,2:end); % get rid of first column with time point labels
            % if there are missing nodes, exclude them here
            if exist('missing_nodes','var')
                roi_data2_tmp2(:,missing_nodes) = []; % get rid of missing nodes
            else
            end
            
            roi_data1_tmp = importdata([homedir dataset '/hcp1200/' num2str(subid(sub)) '/' taskdirnames{t} '/RL/matrices/' num2str(subid(sub)) '_' taskdirnames{t} '_RL_NOGSR_roimean.txt']);
            roi_data1_tmp2 = roi_data1_tmp.data(:,2:end); % get rid of first column with time point labels
            % if there are missing nodes, exclude them here
            if exist('missing_nodes','var')
                roi_data1_tmp2(:,missing_nodes) = []; % get rid of missing nodes
            else
            end
        elseif strcmp(dataset,'abcd')
            abcd_files = dir([homedir dataset '/' subid{sub} '/*' taskdirnames{t} '*_roimean.txt']); % all run without GSR
            [~,abcd_order] = sort({abcd_files.name}); % just to be sure that dir is in alphanumeric order to match run # to timestamp
            abcd_files = abcd_files(abcd_order);
            roi_data2_tmp = importdata([abcd_files(2).folder '/' abcd_files(2).name]);
            roi_data2_tmp2 = roi_data2_tmp.data(:,2:end); % get rid of first column with time point labels
            % if there are missing nodes, exclude them here
            if exist('missing_nodes','var')
                roi_data2_tmp2(:,missing_nodes) = []; % get rid of missing nodes
            else
            end
            roi_data1_tmp = importdata([abcd_files(1).folder '/' abcd_files(1).name]);
            roi_data1_tmp2 = roi_data1_tmp.data(:,2:end); % get rid of first column with time point labels
            % if there are missing nodes, exclude them here
            if exist('missing_nodes','var')
                roi_data1_tmp2(:,missing_nodes) = []; % get rid of missing nodes
            else
            end
        end
        
        % zscore the data within subject
        roi_data_2 = zscore(roi_data2_tmp2);
        roi_data_1 = zscore(roi_data1_tmp2);
        
        %% TASK ACTIVATION ONLY
        % just do regression with task timing regressor to compare to
        % previously published task activation maps
        % %         for i = 1:size(roi_data_2,2)
        % %             betas_task_lr(i,sub,:) = regress(roi_data_2(:,i),[ones(size(roi_data_2,1),1) reg_conv_2]);
        % %         end
        % %
        % %         for i = 1:size(roi_data_1,2)
        % %             betas_task_rl(i,sub,:) = regress(roi_data_1(:,i),[ones(size(roi_data_1,1),1) reg_conv_1]);
        % %         end
        % %
        % %         if t==1 % emotion
        % %             task_act_lr = squeeze(betas_task_lr(:,sub,3))-squeeze(betas_task_lr(:,sub,4));
        % %             task_act_rl = squeeze(betas_task_rl(:,sub,3))-squeeze(betas_task_rl(:,sub,4));
        % %             task_act(:,sub) = (task_act_lr+task_act_rl)./2;
        % %         elseif t==2 % gambling
        % %             task_act_lr = squeeze(betas_task_lr(:,sub,3))-squeeze(betas_task_lr(:,sub,2));
        % %             task_act_rl = squeeze(betas_task_rl(:,sub,3))-squeeze(betas_task_rl(:,sub,2));
        % %             task_act(:,sub) = (task_act_lr+task_act_rl)./2;
        % %         elseif t==3 %motor - skip for now, not sure what to do with contrast so just save out all betas
        % %         elseif t==4 % relational
        % %             task_act_lr = squeeze(betas_task_lr(:,sub,4))-squeeze(betas_task_lr(:,sub,3));
        % %             task_act_rl = squeeze(betas_task_rl(:,sub,4))-squeeze(betas_task_rl(:,sub,3));
        % %             task_act(:,sub) = (task_act_lr+task_act_rl)./2;
        % %         elseif t==5 % social
        % %             task_act_lr = squeeze(betas_task_lr(:,sub,2))-squeeze(betas_task_lr(:,sub,3));
        % %             task_act_rl = squeeze(betas_task_rl(:,sub,2))-squeeze(betas_task_rl(:,sub,3));
        % %             task_act(:,sub) = (task_act_lr+task_act_rl)./2;
        % %         elseif t==6 % wm
        % %             task_act_lr = squeeze(betas_task_lr(:,sub,3))-squeeze(betas_task_lr(:,sub,4));
        % %             task_act_rl = squeeze(betas_task_rl(:,sub,3))-squeeze(betas_task_rl(:,sub,4));
        % %             task_act(:,sub) = (task_act_lr+task_act_rl)./2;
        % %         end
        % %         clear task_act_*
        
        %% FULL REGRESSION (with terms dropped to assess collinearity)
        for i = 1:size(roi_data_2,2) % i = predictor node
            % calculate interaction regressor
            reg_int_2 = roi_data_2(:,i).*reg_conv_bin_2;
            for j = 1:size(roi_data_2,2) % j = target/outcome node
                betas_2(j,i,sub,:) = regress(roi_data_2(:,j),[ones(size(roi_data_2,1),1) roi_data_2(:,i) reg_conv_2 reg_int_2]);
                betas_noNodeTimecourse_2(j,i,sub,:) = regress(roi_data_2(:,j),[ones(size(roi_data_2,1),1) reg_conv_2 reg_int_2]);
                betas_noInteraction_2(j,i,sub,:) = regress(roi_data_2(:,j),[ones(size(roi_data_2,1),1) roi_data_2(:,i) reg_conv_2]);
            end
            clear reg_int_2
        end
        % now do same as above but for RL data
        for i = 1:size(roi_data_1,2) % i = predictor node
            % calculate interacton regressor
            reg_int_1 = roi_data_1(:,i).*reg_conv_bin_1;
            for j = 1:size(roi_data_1,2) % j = target/outcome node
                betas_1(j,i,sub,:) = regress(roi_data_1(:,j),[ones(size(roi_data_1,1),1) roi_data_1(:,i) reg_conv_1 reg_int_1]);
                betas_noNodeTimecourse_1(j,i,sub,:) = regress(roi_data_1(:,j),[ones(size(roi_data_1,1),1) reg_conv_1 reg_int_1]);
                betas_noInteraction_1(j,i,sub,:) = regress(roi_data_1(:,j),[ones(size(roi_data_1,1),1) roi_data_1(:,i) reg_conv_1]);
            end
            clear reg_int_1
        end
        
        clear roi_* task_reg* reg_conv* targ_nodes
        
    end
    % so betas(:,:,:,1) = intercept, betas(:,:,:,2) = "intrinsic", betas(:,:,:,3) =
    % task, betas(:,:,:,4) = interaction; first dim = target/outcome node, second dim =
    % predictor node
    
    save([homedir dataset '/results/betas/betas_',tasks{t},'_zscore_withinSub_NodeTimecourses_3models_zeroCenteredTaskReg_convInterxn_condCondition_withCueReg_noGSR_allBetas_withCue.mat'],'betas_2','betas_1','betas_noNodeTimecourse_1','betas_noNodeTimecourse_2','betas_noInteraction_2','betas_noInteraction_1','-v7.3');
    clear tmp* cue_2 cue_1 betas*
    %     clear task_act*
end
end
