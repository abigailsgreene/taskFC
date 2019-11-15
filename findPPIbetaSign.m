% Code to find mean (across subjects) sign of PPI betas for given feature

taskname = {'gambling' 'emotion' 'wm' 'social' 'relational' 'language' 'motor'};
signmat = zeros(268,268,3,length(taskname));
signvec = zeros(268,length(taskname));
homedir = '/data17/mri_group/abby_data/ppi_paper/';
dataset = 'hcp';
symmetrize = 1;

for tasknum = 1:length(taskname) % iterate through all tasks
    % load betas for given task
    if strcmp(taskname{tasknum},'language')
        load([homedir dataset '/results/betas/betas_' taskname{tasknum} '_zscore_withinSub_NodeTimecourses_3models_zeroCenteredTaskReg_convInterxn_condCondition_withCueReg_noGSR_allBetas_withCue.mat'],'betas_1','betas_2');
    else
        load([homedir dataset '/results/betas/betas_' taskname{tasknum} '_zscore_withinSub_NodeTimecourses_3models_zeroCenteredTaskReg_convInterxn_condCondition_withCueReg_noGSR_allBetas.mat'],'betas_1','betas_2');
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
   
    % iterate through and set mean positive edges to 1, negative edges to -1
    for term = 1:3
        for e = 1:size(all_mats,1)
            for ee = 1:size(all_mats,2)
                if mean(squeeze(squeeze(all_mats(e,ee,:,term))))>0
                    signmat(e,ee,term,tasknum) = 1;
                else signmat(e,ee,term,tasknum) = -1;
                end
            end
        end
    end
    % now do same for nodes and mean activation
    for e = 1:size(all_mats,1)
        if mean(squeeze(all_vecs(e,:)))>0
            signvec(e,tasknum) = 1;
        else signvec(e,tasknum) = -1;
        end
    end
    
end
save([homedir 'hcp/ppiBeta_sign.mat','signmat','signvec']);
% code to count up number of positive and negative (unique) edges in
% signmat, if you're so inclined
upp_idx = find(triu(ones(268,268),1));
for i = 1:7
    for j = 1:3
        tmpmat = squeeze(signmat(:,:,j,i));
        pos(i,j) = length(find(tmpmat(upp_idx)==1));
        neg(i,j) = length(find(tmpmat(upp_idx)==-1));
        clear tmpmat
    end
end
% and to count up number of positive and negative nodes for each task
for i = 1:7
    pos_task(i) = length(find(signvec(:,i)==1));
    neg_task(i) = length(find(signvec(:,i)==-1));
end
