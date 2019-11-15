% Code to load in results from rCPM and identify features selected on more
% than npct of iterations + calculate their contributions
% written by Abby Greene, 2019
homedir = '/data17/mri_group/abby_data/ppi_paper/'; % where is the parent directory for inputs and outputs?
dataset = 'hcp'; % dataset?
taskname = {'wm','language','emotion','gambling','relational','social','motor'}; % task names?
gsr = 'noGSR'; % did we do GSR?
symmetrize = 1; % did we symmetrize the beta matrices before inputting them to ridge prediction? 1 = yes, 0 = no
numiters = '1';
pthresh = '0.1';
npct = 0.75; %percentage of the iterations on which we require feature to be selected in order to analyze it

% iterate through all tasks
for t = 1:length(taskname)
    % initialize variable to store indices of selected edges for each task
    selected_edges = cell(1,5);
    selected_edges_betas = cell(1,5);
    edge_mat = cell(1,5);
    coef_total_reshape = cell(1,5);
    
    % load betas for given task
    if strcmp(taskname{t},'language')
        data_tmp = load([homedir dataset '/results/betas/betas_' taskname{t} '_zscore_withinSub_NodeTimecourses_3models_zeroCenteredTaskReg_convInterxn_condCondition_withCueReg_noGSR_allBetas_withCue.mat'],'betas_1','betas_2'); % if language, load in version where we included cue regressor
    else
        data_tmp = load([homedir dataset '/results/betas/betas_' taskname{t} '_zscore_withinSub_NodeTimecourses_3models_zeroCenteredTaskReg_convInterxn_condCondition_withCueReg_noGSR_allBetas.mat'],'betas_1','betas_2');
    end
    if strcmp(dataset,'abcd')
        load([homedir 'abcd/abcdData_noMissingPMAT.mat'],'idx_pmat');
        betas_1 = data_tmp.betas_1(:,:,idx_pmat,:);
        betas_2 = data_tmp.betas_2(:,:,idx_pmat,:);
    else
        betas_1 = data_tmp.betas_1;
        betas_2 = data_tmp.betas_2;
    end
    % calculate condition contrast betas
    subnum = size(betas_2,3);
    % regardless of task, 1st beta is always intercept and 2nd beta is
    % always timecourse/intrinsic FC
    for i = 1:subnum
        all_mats(:,:,i,1) = (squeeze(betas_2(:,:,i,1))+squeeze(betas_1(:,:,i,1)))./2; % intercept
        all_mats(:,:,i,2) = (squeeze(betas_2(:,:,i,2))+squeeze(betas_1(:,:,i,2)))./2; % node time course
    end
    % all_mats(:,:,:,3) = cdFC, all_vecs = node x sub activation
    if strcmp(taskname{t},'wm') %wm
        for i = 1:subnum
            all_vecs(:,i) = mean((((squeeze(betas_2(:,:,i,4))-squeeze(betas_2(:,:,i,5)))+(squeeze(betas_1(:,:,i,4))-squeeze(betas_1(:,:,i,5))))./2),2); % task contrast (2bk-0bk)
            all_mats(:,:,i,3) = ((squeeze(betas_2(:,:,i,7))-squeeze(betas_2(:,:,i,8)))+(squeeze(betas_1(:,:,i,7))-squeeze(betas_1(:,:,i,8))))./2; %interaction contrast (2bk-0bk)
        end
    elseif strcmp(taskname{t},'language')
        for i = 1:subnum
            all_vecs(:,i) = mean((((squeeze(betas_2(:,:,i,4)))+(squeeze(betas_1(:,:,i,4))))./2),2); % task contrast (language vs. math [modeling only language because no rest])
            all_mats(:,:,i,3) = ((squeeze(betas_2(:,:,i,6)))+(squeeze(betas_1(:,:,i,6))))./2; %interaction contrast (language vs. math)
        end
    elseif strcmp(taskname{t},'emotion')
        for i = 1:subnum
            all_vecs(:,i) = mean(((squeeze(betas_2(:,:,i,4))+squeeze(betas_1(:,:,i,4)))./2),2); % just fear (vs neut = baseline)
            all_mats(:,:,i,3) = (squeeze(betas_2(:,:,i,6))+squeeze(betas_1(:,:,i,6)))./2; % just fear interaction (vs neut = baseline)
        end
    elseif strcmp(taskname{t},'gambling') % gambling
        for i = 1:subnum
            all_vecs(:,i) = mean((((squeeze(betas_2(:,:,i,4))-squeeze(betas_2(:,:,i,3)))+(squeeze(betas_1(:,:,i,4))-squeeze(betas_1(:,:,i,3))))./2),2); % task contrast (gam: win-loss)
            all_mats(:,:,i,3) = ((squeeze(betas_2(:,:,i,6))-squeeze(betas_2(:,:,i,5)))+(squeeze(betas_1(:,:,i,6))-squeeze(betas_1(:,:,i,5))))./2; %interaction contrast (gam: win-loss)
        end
    elseif strcmp(taskname{t},'relational') % relational
        for i = 1:subnum
            all_vecs(:,i) = mean((((squeeze(betas_2(:,:,i,5))-squeeze(betas_2(:,:,i,4)))+(squeeze(betas_1(:,:,i,5))-squeeze(betas_1(:,:,i,4))))./2),2); % task contrast (relational: rel-match)
            all_mats(:,:,i,3) = ((squeeze(betas_2(:,:,i,8))-squeeze(betas_2(:,:,i,7)))+(squeeze(betas_1(:,:,i,8))-squeeze(betas_1(:,:,i,7))))./2; %interaction contrast (relational: rel-match)
        end
    elseif strcmp(taskname{t},'social') % social
        for i = 1:subnum
            all_vecs(:,i) = mean((((squeeze(betas_2(:,:,i,3))-squeeze(betas_2(:,:,i,4)))+(squeeze(betas_1(:,:,i,3))-squeeze(betas_1(:,:,i,4))))./2),2); % task contrast (mental-rand)
            all_mats(:,:,i,3) = ((squeeze(betas_2(:,:,i,5))-squeeze(betas_2(:,:,i,6)))+(squeeze(betas_1(:,:,i,5))-squeeze(betas_1(:,:,i,6))))./2; %interaction contrast (mental-rand)
        end
    elseif strcmp(taskname{t},'motor') % for motor, just moving (add all conditions) vs rest (unmodeled)
        for i = 1:subnum
            all_vecs(:,i) = mean((((squeeze(betas_2(:,:,i,4))+squeeze(betas_2(:,:,i,5))+squeeze(betas_2(:,:,i,6))+squeeze(betas_2(:,:,i,7))+squeeze(betas_2(:,:,i,8)))+(squeeze(betas_1(:,:,i,4))+squeeze(betas_1(:,:,i,5))+squeeze(betas_1(:,:,i,6))+squeeze(betas_1(:,:,i,7))+squeeze(betas_1(:,:,i,8))))./2),2);
            all_mats(:,:,i,3) = ((squeeze(betas_2(:,:,i,10))+squeeze(betas_2(:,:,i,11))+squeeze(betas_2(:,:,i,12))+squeeze(betas_2(:,:,i,13))+squeeze(betas_2(:,:,i,14)))+(squeeze(betas_1(:,:,i,10))+squeeze(betas_1(:,:,i,11))+squeeze(betas_1(:,:,i,12))+squeeze(betas_1(:,:,i,13))+squeeze(betas_1(:,:,i,14))))./2;
        end
    elseif strcmp(taskname{t},'mid') % win - loss
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
    
    % vectorize betas, taking only lower triangle (already have node x
    % subject for task activation)
    for i = 1:size(all_mats,4)
        for sub = 1:size(all_mats,3)
            all_mats_vec_tmp = squeeze(all_mats(:,:,sub,i));
            all_mats_vec(:,sub,i) = all_mats_vec_tmp(find(tril(ones(size(all_mats,1),size(all_mats,2)),-1))); % edges x subject x PPI term
            clear *_tmp
        end
    end
    
    num_edges = size(all_mats_vec,1);
    num_nodes = size(all_mats,1);
    
    % load ridge CPM results
%     if strcmp(taskname{t},'language')
%         load([homedir dataset '/results/' dataset '_' taskname{t} '_ridge_p' pthresh '_' numiters 'iters_pmat_prediction_' gsr '_withCue.mat']);
%     else
       load([homedir dataset '/results/motionRegression/' dataset '_' taskname{t} '_ridge_p' pthresh '_' numiters 'iters_10fold_pmat_prediction_' gsr '_withCue_preAllocatedCoeff.mat'],'coef_total_resid'); 
        coef_total = coef_total_resid;
       %     end
    % iterate through 5 terms (intercept, ciFC, cdFC, activation, tot)
    for term = 1:size(coef_total,1)
        clear coef_bin coef_tot
        for it = 1:size(coef_total,2) % and through all iterations
            % now identify edges that are selected in every fold on every iteration
            coef_bin(:,:,it) = coef_total{term,it}; % will also store selected features in binary form to count number of times feature was selected below
            coef_tot(:,:,it) = coef_total{term,it}; % regression coefficients for all selected (non-selected = 0) edges for given iteration
            % get ridge betas into edge x total iter format to plot distribution of betas for selected edges
            coef_total_reshape{term} = [coef_total_reshape{term} coef_total{term,it}];
        end
        coef_bin(abs(coef_bin)>0) = 1; % binarize the betas so it's just selected or not
        for i = 1:size(coef_bin,1) % iterate through every edge
            if sum(sum(coef_bin(i,:,:)))>=npct*(size(coef_bin,2)*size(coef_bin,3)) % if sum for each edge over all folds (dim2) and iterations (dim3) = number of folds* number of iterations, it was selected npct% of time
                selected_edges{term} = [selected_edges{term}; i]; % store indices of selected features for given task and term
                selected_edges_betas{term}(i) = mean(mean(squeeze(coef_tot(i,:,:)))); % vector of all edges (or nodes for task act), with mean (across folds and iterations) betas for selected betas, 0 for not-selected
            else
                selected_edges_betas{term}(i) = 0;
            end
        end
        % now reshape these selected edges' contributions (ridge betas * std of input
        % PPI betas) into nodexnode matrices to get contribution of each
        % feature
        if term < 4 % ie intercept, intrinsic, or interaction
            edge_mat{term} = squareform(selected_edges_betas{term}'.*std(squeeze(all_mats_vec(:,:,term)),0,2));
            edge_vec{term} = sum(edge_mat{term});
        elseif term==4 % ie task activation
            edge_mat{term} = selected_edges_betas{term}'.*std(all_vecs,0,2); % for task activation, just want nodex1 vector with avg betas for selected nodes * std of that node across subjects
            edge_vec{term} = edge_mat{term};
        elseif term==5 % all terms
            edge_mat{term}(:,:,1) = squareform(selected_edges_betas{term}(1:size(coef_total{1,1},1))'.*std(squeeze(all_mats_vec(:,:,1)),0,2));
            edge_mat{term}(:,:,2) = squareform(selected_edges_betas{term}(size(coef_total{1,1},1)+1:2*size(coef_total{1,1},1))'.*std(squeeze(all_mats_vec(:,:,2)),0,2));
            edge_mat{term}(:,:,3) = squareform(selected_edges_betas{term}(2*size(coef_total{1,1},1)+1:3*size(coef_total{1,1},1))'.*std(squeeze(all_mats_vec(:,:,3)),0,2));
            edge_mat{term}(:,:,4) = diag(selected_edges_betas{term}(3*size(coef_total{1,1},1)+1:end)'.*std(all_vecs,0,2)); % for task activation, just put nodes on diagonal of matrix for convenience
            for tt = 1:4 % vectorize this information, if you're interested
               edge_vec{term}(:,tt) = sum(squeeze(edge_mat{term}(:,:,tt))); 
            end
            
            % now calculate term-level contributions
            low_idx = find(tril(ones(268,268))); % only look at contributions from lower triangle, to not duplicate
            for i = 1:3
                % number of selected features from this term normalized
                % by number of features in combined model, so
                % contrib = percentage of features in combined model that
                % are from each term
                contrib(i) = (length(find(edge_mat{term}(:,:,i)))/2)/(length(find(selected_edges_betas{term})));
                
                % now get weighted term contribution
                tmp_mat = squeeze(edge_mat{term}(:,:,i));
                if ~isempty(find(tmp_mat)) % get total contribution of given term
                    beta_contrib(i) = (sum(sum(abs(tmp_mat(low_idx))))); % changed this on 9/9 to only consider lower triangle to not duplicate contributions for given edge
                else
                    beta_contrib(i) = 0;
                end
               
                clear tmp_mat
            end
            % now do same for activation
            contrib(4) = (length(find(edge_mat{term}(:,:,4))))/(length(find(selected_edges_betas{term})));
            tmp_mat = squeeze(edge_mat{term}(:,:,4));
            if ~isempty(find(tmp_mat))
                beta_contrib(4) = (sum(sum(abs(tmp_mat)))); % don't have to worry about duplicating, since values are only on the diagonal
            else
                beta_contrib(4) = 0;
            end
            
            if sum(beta_contrib)~=0 % if no selected edges and beta_contrib is all 0's, skip this normalization (else will get NaNs)
                for ii = 1:4 % make wted contribution into percentages
                    beta_contrib_pct(ii) = beta_contrib(ii)/sum(beta_contrib);
                end
            else
            end
            clear tmp_mat
        end
        % store predictive utility for all features in summary variable
        % across tasks and terms
        selected_edges_betas_tot{t,term} = selected_edges_betas{term};   
    end
    % store overall predictive utility of each term in summary variable
    % across tasks
    beta_contrib_pct_tot{t} = beta_contrib_pct;
    %save([homedir dataset '/results/' dataset '_p' pthresh '_' numiters 'iters_' gsr '_' taskname{t} 'ridge_visualization_noNormBetaContrib_withContribVec_withLangCue_lowerTriangleContribOnly.mat'],'edge_vec','selected_edges','selected_edges_betas','edge_mat','contrib','beta_contrib','beta_contrib_pct','taskname','npct','pthresh','numiters');
    clearvars -except orig parcorr homedir dataset gsr taskname symmetrize numiters pthresh npct beta_contrib_pct_tot* selected_edges_betas_tot*
end

