% This script launch the combination of all choices data from GLIOMA FATIGUE
% It also analyses the proportion of choices
%
% Author: Antonius Wiehler <antonius.wiehler@gmail.com>
% Original: 2018-04-09, copied from TIMECONTROL
% Modified: 2020-02-27, Valentine Facque <facque.valentine@gmail.com>


%% PREPARATION
% =========================================================================
clear all;
close all;
clc;


%% SET CONFIGURATION PARAMETERS
% =========================================================================
% this should be everything that is used in multiple scripts.
cfg.studyname = 'GLIOMA_FATIGUE';

cfg.conditions = {'GLIOMAFATIGUE01control', 'GLIOMAFATIGUE01preOP', 'GLIOMAFATIGUE01postOP'};

task_list(1, 1) = {'TD_ID'};  % time discounting imm-delayed
task_list(1, 2) = {'TD_DD'};  % time discounting delayed-delayed

task_list = reshape(task_list', [], 1);  % reshape into one line of tasks, to undo this use reshape(task_list, 2, [])'

cfg.mainexp.n_blocks   = 23;  % from main testing script
cfg.mainexp.n_sessions = 1;

cfg.CREAexp.n_blocks   = 4;  % from main testing script
cfg.CREAexp.n_sessions = 1;


% directories
cfg.dir.meta         = 'outputs/meta/';
cfg.dir.choice       = 'outputs/choice/mainexp/';
cfg.dir.choiceCREA   = 'outputs/choice/CREAexp/';
cfg.dir.choicecalib  = 'outputs/choice/calibrations/';
cfg.dir.choiceprep   = 'outputs/choice/choiceprep/';
cfg.dir.choice_plots = 'analyze/choice/diagnostic_plots/';
cfg.dir.daystructure = 'outputs/day_structure/';
cfg.dir.rating       = 'outputs/ratings/';

cfg.dir.group        = 'analyze/choice/group/';
cfg.dir.group_plots  = 'analyze/choice/group/plots/';
cfg.dir.group_subject  = 'analyze/choice/group/subject/';

create_missing_directories(cfg.dir);

cfg.choice.indiff_min_distance = 10;  % minimum distance of indifference point from floor and ceiling to be considered as a good delay

% add toolboxes
[~, hostname]=system('hostname');
hostname = deblank(hostname);

if strcmp(hostname, 'UMR-PESSI-WP001')
    cfg.dir.export_fig   = 'C:/Users/student/Documents/Matlab_Toolbox/export_fig-master';
    cfg.dir.gramm = 'C:/Users/student/Documents/Matlab_Toolbox/gramm-master';
    cfg.dir.violin = 'C:/Users/valentine.facque/Documents/Matlab_Toolbox/Violinplot-master';
    addpath(genpath(cfg.dir.export_fig));  % add export_fig for saving plots
    addpath(genpath(cfg.dir.gramm)); % add gramm for easier plots
    addpath(genpath(cfg.dir.violin)); %add violin for violin plots
end

% figure
cfg.fig.width              = 20; % cm
cfg.fig.height             = 10;   % cm
cfg.fig.width_per_subplot  = cfg.fig.width / 4.5; % cm
cfg.fig.height_per_subplot = cfg.fig.height / 2; % cm
cfg.fig.height             = 10;   % cm
cfg.fig.fsiz               = 10;    % font size
cfg.fig.offset             = 0.5;  % offset of x axis to better plot two groups

cfg.fig.colormap     = [ 0    0.4470    0.7410;  % blue
    0.8500    0.3250    0.0980;  % red/orange
    0.9290    0.6940    0.1250;  % yellow
    0.4940    0.1840    0.5560;  % purple
    0.4660    0.6740    0.1880;  % green
    0.3010    0.7450    0.9330;  % light blue
    0.6350    0.0780    0.1840];  % dark red



%% CALIBRATION - COMBINE ALL SUBJECTS
% =========================================================================

for i_t = 1 : length(task_list)
    combine_calibration_data_group_v01(cfg, task_list{i_t});
end  % end for loop tasks

savewhere = strcat(cfg.dir.group, cfg.studyname, '_group_calibration_data');  % whre to save the final table?

if ~exist([savewhere '.mat'], 'file')  % only run, when output does not exist
    calib_data = [];  % to add data later
    
    for i_t = 1 : length(task_list)
        filename = strcat(cfg.dir.group, task_list{i_t}, '_group_calibration_data');  % whre to save the final table?
        tmp = load(filename);
        calib_data = [calib_data; tmp.calib_data];  % add data of current task
    end  % end for loop tasks
    
    calib_data.subject_id (:) = extractBefore(calib_data.subject_id(:),6);
    calib_data.raw_response_time = calib_data.response_time;
    
    save([savewhere '.mat'], 'calib_data');  % save as matlab data
    writetable(calib_data, [savewhere '.csv']);  % save as csv table
    
end  % if file does not exist


%% MAINEXP - COMBINE THE DATA OF ALL SUBJECTS PER TASK
% =========================================================================

savewhere = strcat(cfg.dir.group, cfg.studyname, '_group_choicemain_data');  % whre to save the final table?
if ~exist([savewhere '.mat'], 'file')  % only run, when output does not exist
    
    for i_t = 1 : length(task_list)
        combine_mainexp_choice_data_group_v01(cfg, task_list{i_t}, cfg.conditions);
    end  % end for loop tasks
    
    choicemain_data = [];  % to add data later
    
    for i_t = 1 : length(task_list)
        filename = strcat(cfg.dir.group, task_list{i_t}, '_group_choicemain_data');  % whre to save the final table?
        tmp = load(filename);
        choicemain_data = [choicemain_data; tmp.choicemain_data];  % add data of current task
    end  % end for loop tasks
    
    % Z-score RT per subject
    subjects = unique(choicemain_data.subvisit);
    choicemain_data.raw_response_time = choicemain_data.response_time;
    
    for i_s = 1 : length(subjects)
        idx_subject = strcmp(choicemain_data.subvisit, subjects{i_s});
        choicemain_data.response_time(idx_subject) = zscore(choicemain_data.response_time(idx_subject));
    end
    
    % split blocks for easier plot
    choicemain_data.block_split = zeros(height(choicemain_data), 1);
    split_point = median(unique(choicemain_data.block));
    choicemain_data(choicemain_data.block <= split_point, 'block_split') = {3};
    choicemain_data(choicemain_data.block > split_point, 'block_split') = {4};
    
    
    save([savewhere '.mat'], 'choicemain_data');  % save as matlab data
    writetable(choicemain_data, [savewhere '.csv']);  % save as csv table
    
end  % if file does not exist


%% CREA - COMBINE THE DATA OF ALL SUBJECTS PER TASK
% =========================================================================
savewhere = strcat(cfg.dir.group, cfg.studyname, '_group_choiceCREA_data');  % whre to save the final table?

if ~exist([savewhere '.mat'], 'file')  % only run, when output does not exist
    choiceCREA_data = [];  % to add data later
    
    for i_t = 1 : length(task_list)
        combine_CREAexp_choice_data_group_v01(cfg, task_list{i_t}, cfg.conditions);
    end  % end for loop tasks
    
    for i_t = 1 : length(task_list)
        filename = strcat(cfg.dir.group, task_list{i_t}, '_group_choiceCREA_data');  % whre to save the final table?
        tmp = load(filename);
        choiceCREA_data = [choiceCREA_data; tmp.choiceCREA_data];  % add data of current task
    end  % end for loop tasks
    
    % Z-score RT per subject
    subjects = unique(choiceCREA_data.subvisit);
    choiceCREA_data.raw_response_time = choiceCREA_data.response_time;
    
    
    for i_s = 1 : length(subjects)
        idx_subject = strcmp(choiceCREA_data.subvisit, subjects{i_s});
        choiceCREA_data.response_time(idx_subject) = zscore(choiceCREA_data.response_time(idx_subject));
    end
    
    % split blocks for easier plot
    choiceCREA_data.block_split = zeros(height(choiceCREA_data), 1);
    split_point = median(unique(choiceCREA_data.block));
    choiceCREA_data(choiceCREA_data.block <= split_point, 'block_split') = {1};
    choiceCREA_data(choiceCREA_data.block > split_point, 'block_split') = {2};
    
    save([savewhere '.mat'], 'choiceCREA_data');  % save as matlab data
    writetable(choiceCREA_data, [savewhere '.csv']);  % save as csv table
    
end  % if file does not exist


%% MAINEXP&CREA - COMBINE THE DATA FROM BOTH TASK
% =================================================
savewhere = strcat(cfg.dir.group, cfg.studyname, '_group_choiceALL_data');  % where to save the final table?

if ~exist([savewhere '.mat'], 'file')  % only run, when output does not exist
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceCREA_data.mat');  % where is the bic choice data table
    load(filename);  % load choice_data
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choicemain_data.mat');  % where is the bic choice data table
    load(filename);  % load choice_data
    choiceALL_data = [];
    
    %distinguish part of experiment
    choiceCREA_data.part = cellstr(repmat('1',size(choiceCREA_data.choiceSS)));
    choicemain_data.part = cellstr(repmat('2',size(choicemain_data.choiceSS)));
    
    choiceALL_data = [choiceCREA_data ; choicemain_data];
    
    choiceALL_data.subject_id (:) = extractBefore(choiceALL_data.subject_id(:),6);
    
    
    chtot = choiceALL_data;
    
    % save data for all subjects
    save([savewhere '.mat'],'chtot');  % save as matlab data
    writetable(chtot, [savewhere '.csv']);  % save as csv table
end

%% CALIB&CHOICE - COMBINE THE DATA FROM CALIBRATION & CHOICE ALL
% =================================================
savewhere = strcat(cfg.dir.group, cfg.studyname, '_group_choice&calib_data');  % where to save the final table?

if ~exist([savewhere '.mat'], 'file')  % only run, when output does not exist
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceALL_data.mat');  % load group choice
    chALL = load(filename);  % load choice_data
    chALL = chALL.chtot;
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_calibration_data');  % load group calibration
    tmp_data = load(filename); % load calib_data
    tmp_data = tmp_data.calib_data;
    
    tmp_data.block_split(:) = 0;
    
    chALL.session = [];
    chALL.block= [];
    chALL.within_block= [];
    chALL.trial_type= [];
    chALL.good_delay= [];
    chALL.empirical_indiff_point= [];
    chALL.empirical_DV= [];
    chALL.part= [];
    chtot = [tmp_data;chALL];
    
    for i_l = 1 : height(chtot)
        %TODO SPLIT SPLAT SPLAT
        if strcmp(chtot.version_test(i_l),'old')
            if chtot.block_split(i_l) == 0
                chtot.part(i_l) = 1;
            elseif chtot.block_split(i_l) == 3 | chtot.block_split(i_l) == 4
                chtot.part(i_l) = 2;
            else
            end
            
        elseif strcmp(chtot.version_test(i_l),'new')
            if chtot.block_split(i_l) == 0
                chtot.part(i_l) = 0;
            elseif chtot.block_split(i_l) == 1 | chtot.block_split(i_l) == 2
                chtot.part(i_l) = 1;
            elseif chtot.block_split(i_l) == 3 | chtot.block_split(i_l) == 4
                chtot.part(i_l) = 2;
            else
            end
        end
    end
    
    %rename block_split into parts
    chtot.part_n(chtot.part == 0) = cellstr('Calib');
    chtot.part_n(chtot.part == 1) = cellstr('HOC');
    chtot.part_n(chtot.part == 2) = cellstr('Switch');
    
    %str2double for several
    chtot.chimio = str2double(chtot.chimio);
    chtot.radio = str2double(chtot.radio);
    chtot.FSS_score = str2double(chtot.FSS_score);
    chtot.HAD_Ascore = str2double(chtot.HAD_Ascore);
    chtot.HAD_Dscore = str2double(chtot.HAD_Dscore);
    chtot.STARK_score = str2double(chtot.STARK_score);
    chtot.BIS_score = str2double(chtot.BIS_score);
    chtot.vol_preOP = str2double(chtot.vol_preOP);
    chtot.annoy_calib = str2double(chtot.annoy_calib);
    chtot.annoy_crea = str2double(chtot.annoy_crea);
    chtot.annoy_nswitch = str2double(chtot.annoy_nswitch);
    
    %reencode diploma
    chtot.diploma_score(:) = 0;
    chtot.diploma_score(strcmp(chtot.last_diploma,'Aucun'),:) = 0;
    chtot.diploma_score(strcmp(chtot.last_diploma,'CAP, BEP'),:) = 1;
    chtot.diploma_score(strcmp(chtot.last_diploma,'BAC'),:) = 2;
    chtot.diploma_score(strcmp(chtot.last_diploma,'BAC + 2'),:) = 3;
    chtot.diploma_score(strcmp(chtot.last_diploma,'BAC + 3'),:) = 4;
    chtot.diploma_score(strcmp(chtot.last_diploma,'BAC + 4'),:) = 5;
    chtot.diploma_score(strcmp(chtot.last_diploma,'BAC + 5'),:) = 6;
    chtot.diploma_score(strcmp(chtot.last_diploma,'BAC + 8'),:) = 7;
    
    chtot.real_c = strcat(chtot.new_condition,'_',chtot.condition);
    
    %rename side
    chtot.r_arrow = chtot.choice_side == 37;
    %same press as before
    chtot.diff_press(:) = nan;
    sub = unique(chtot.subvisit);
    for i_s = 1 : length(sub)
        idx_sub = strcmp(chtot.subvisit,sub{i_s}); %indexing subject
        
        tsk = unique(chtot(idx_sub,:).task);
        for i_t = 1 : length(tsk)
            idx_tsk = strcmp(chtot.task,tsk{i_t}); %indexing task
            idx_s_t = idx_sub & idx_tsk; %combine idx subject and task
            
            blk = unique(chtot(idx_s_t,:).block_split);
            for i_b = 1 : length(blk)
                idx_blk = chtot.block_split == blk(i_b); %indexing block
                idx_s_t_b = idx_sub & idx_tsk & idx_blk; %combine idx subject and task
                
                max_trial = max(chtot(idx_s_t_b,:).trial_number);
                min_trial = min(chtot(idx_s_t_b,:).trial_number);
                for i_tn = (min_trial + 1) : max_trial
                    idx_tn = chtot.trial_number == i_tn ;
                    idx_tn_1 = chtot.trial_number == (i_tn-1);
                    idx_s_t_b_tn = idx_sub & idx_tsk & idx_blk & idx_tn;
                    idx_s_t_b_tn_1 = idx_sub & idx_tsk & idx_blk & idx_tn_1;
                    chtot.diff_press(idx_s_t_b_tn) = chtot.r_arrow(idx_s_t_b_tn) == chtot.r_arrow(idx_s_t_b_tn_1);
                    
                end %loop trial number
                
            end %loop block
        end %loop task
    end %loop subject
    
    
    
    % save data for all subjects
    save([savewhere '.mat'],'chtot');  % save as matlab data
    writetable(chtot, [savewhere '.csv']);  % save as csv table
    
    chtot_new = chtot(strcmp(chtot.version_test,'new'),:);
    chtot_old = chtot(strcmp(chtot.version_test,'old'),:);
    save([savewhere '_new.mat'],'chtot_new');  % save as matlab data
    save([savewhere '_old.mat'],'chtot_old');  % save as matlab data
    
    %checking outliers and removing them
    save1 = strcat(cfg.dir.group, cfg.studyname, '_group_choice&calib_data_outliers');  % where to save the final table?
    
    if ~exist([save1 '.txt'], 'file')
        chtot_out = chtot;
        chtot_out(chtot_out.block_split == 0,:) = [];
        
        outliers = varfun(@mean ,chtot_out,'InputVariables',{'choiceSS'},'GroupingVariables',{'subvisit','new_condition'});
        out = outliers.mean_choiceSS > 0.9 | outliers.mean_choiceSS < 0.1; %outliers if more than 90% or less than 10% just after calib
        outliers = outliers(out,:);
        outliers = sortrows(outliers,'mean_choiceSS','ascend');
        
        % write result to text file
        diary([save1 '.txt']);
        disp(outliers);
        diary off;
        save(strcat(cfg.dir.group, cfg.studyname, '_outliers_list.mat'),'outliers');  % save as matlab data
    end
    
    
    %removing outliers
    for i_o = 1 : length(outliers.subvisit)
        chtot(strcmp(chtot.subvisit,outliers.subvisit{i_o}),:)=[];
    end
    
    
    chtot_new = chtot(strcmp(chtot.version_test,'new'),:);
    chtot_old = chtot(strcmp(chtot.version_test,'old'),:);
    save([save1 '_new.mat'],'chtot_new');  % save as matlab data
    save([save1 '_old.mat'],'chtot_old');  % save as matlab data
end

%% CALIBRATION - IMPULSIVE CHOICE x CONDITION
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_all_01_calib_impulsive_condition', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_calibration_data.mat');  % where is the bic choice data table
    load(filename);  % load choice_data
    
    g = gramm('x', calib_data.costLL_step, 'y', calib_data.choiceSS, 'color', calib_data.new_condition, 'line', calib_data.task);  % define data
    g.stat_summary('geom',{'point', 'line'});  % this does plot the mean point , 'errorbar'
    g.axe_property('YLim', [0 1],'XTick',[1:5], 'XTickLabel', {'3 days', '1 week', '1 month', '3 months', '1 year'});
    g.set_names('column', '', 'row', '', 'x', 'delay', 'y', 'impulsive choice rate', 'color', 'condition', 'line', 'task');
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    
end % if files does not exist

%% CHECK SENSIBILITY OF CHOICES
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_all_02_choice_sensitive', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 8;  % 5 plus legend
    n_y_subplots = 4;
    
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choice&calib_data_outliers_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    
    % gramm
    g = gramm('x', chtot.block_split, 'y', chtot.choiceSS, 'color', chtot.condition);  % define data
    g.stat_summary('geom', {'point', 'line'}, 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('impulsive choice evolution according to cost')
    g.facet_grid(chtot.task , chtot.costLL, 'scale', 'independent');  % split in subplots
    g.axe_property('XTick', [0:5],'XTickLabel', {'Calib', 'part_2', 'part_3', 'part_4', 'part_5'});
    % g.set_color_options('map', cfg.fig.colormap([1 2 6], :));
    g.set_names('column', '', 'x', ' ', 'y', 'impulsive choice rate', 'color', 'condition');
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end

%% CALIBRATION - REVERSE CHECK
% =========================================================================
savewhere = strcat(cfg.dir.group, cfg.studyname, '_group_choiceREVERSEcalib_data');  % where to save the final table?

%objective : restrict the choices in calib to choices presented in crea and switch run
%choices matches the indifference points of participants
%warning : at the end not the same amount of choices per participants

if ~exist([savewhere '.mat'], 'file')  % only run, when output does not exist
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceALL_data.mat');
    load(filename);  % load choiceALL_data
    chALL = chtot; %save for after
    choiceALL_data = chtot;
    chtot = [];
    
    indx_50 = find(choiceALL_data.rewardSS >= 50); %target all catch trials
    choiceALL_data(indx_50,:) = []; %delete them
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_calibration_data.mat');
    load(filename); % load calib_data
    
    %check the boundaries for each participant according to delay and task
    max_D_part =  varfun(@max,choiceALL_data,'InputVariables',{'rewardSS'},'GroupingVariables',{'subvisit','costLL','task'});
    min_D_part =  varfun(@min,choiceALL_data,'InputVariables',{'rewardSS'},'GroupingVariables',{'subvisit','costLL','task'});
    indif_D_part =  varfun(@mean,choiceALL_data,'InputVariables',{'empirical_indiff_point'},'GroupingVariables',{'subvisit','costLL','task'});
    indif_D_part.max(:) = NaN ;
    indif_D_part.min(:) = NaN ;
    
    %create index to target precisely the participant-task-delay combo
    indif_D_part.idx(:) = strcat(indif_D_part.subvisit(:),'_',indif_D_part.task(:), '_',num2str(indif_D_part.costLL(:)));
    max_D_part.idx(:) = strcat(max_D_part.subvisit(:),'_',max_D_part.task(:), '_',num2str(max_D_part.costLL(:)));
    min_D_part.idx(:) = strcat(min_D_part.subvisit(:),'_',min_D_part.task(:), '_',num2str(min_D_part.costLL(:)));
    
    subjects = unique(indif_D_part.idx); % list all subjects and session
    
    for i_s = 1 : length(subjects)
        
        max_indx = strcmp(max_D_part.idx,subjects{i_s});
        min_indx = strcmp(min_D_part.idx,subjects{i_s});
        
        max_data = max_D_part.max_rewardSS(max_indx);
        min_data = min_D_part.min_rewardSS(min_indx);
        indif_D_part.max(strcmp(max_D_part.idx, subjects{i_s})) = max_data; %localise le numero sujet
        indif_D_part.min(strcmp(min_D_part.idx, subjects{i_s})) = min_data; %localise le numero sujet
    end
    
    tmp_data = calib_data;
    tmp_data.idx(:) = strcat(tmp_data.subvisit(:),'_',tmp_data.task(:), '_',num2str(tmp_data.costLL(:)));
    
    subjects = unique(tmp_data.idx); % list all subjects and session
    for i_d = 1 : length(subjects)
        val_indx = strcmp(indif_D_part.idx,subjects{i_d});
        if sum(val_indx) > 0
            upper = indif_D_part.max(val_indx); %upper value for the participant-task-delay combo
            lower = indif_D_part.min(val_indx); %lower value for the participant-task-delay combo
            mid = indif_D_part.mean_empirical_indiff_point(val_indx); %indif point for the participant-task-delay combo
            
            tmp_data.empirical_indiff_point(strcmp(tmp_data.idx, subjects{i_d})) = mid; %loacate subject number
            
            indx_upper = tmp_data.rewardSS > upper; % loacate all upper values
            indx_lower = tmp_data.rewardSS < lower;  % locate all lower values
            whom = strcmp(tmp_data.idx, subjects{i_d}); %locate the subject
            new_u_idx = and(whom,indx_upper);  %select both the subject and all upper values
            new_l_idx = and(whom,indx_lower);  %select both the subject and all under values
            
            %
            tmp_data(new_u_idx,:) = []; %erase all the calibration proposal out of the spectrum of the choices
            tmp_data(new_l_idx,:) = []; %erase all the calibration proposal out of the spectrum of the choices
            
        else
        end
        
    end
    
    tmp_data.upper = tmp_data.rewardSS > tmp_data.empirical_indiff_point ;
    tmp_data.equal = tmp_data.rewardSS == tmp_data.empirical_indiff_point;
    chALL.upper = chALL.rewardSS > chALL.empirical_indiff_point;
    chALL.equal = chALL.rewardSS == chALL.empirical_indiff_point;
    
    %clear data to combine
    tmp_data.block_split(:) = 0;
    chALL.session = [];
    chALL.block = [];
    chALL.within_block = [];
    chALL.trial_type = [];
    chALL.good_delay = [];
    chALL.empirical_DV = [];
    tmp_data.part (:) = {'0'};
    tmp_data.idx = [];
    chtot = [tmp_data;chALL];
    chtot.part = str2double(chtot.part);
    chtot_new = chtot(strcmp(chtot.version_test,'new'),:);
    
    for i_l = 1 : height(chtot)
        
        if chtot.block_split(i_l) == 0
            chtot.part(i_l) = 0;
        elseif chtot.block_split(i_l) == 1 | chtot.block_split(i_l) == 2
            chtot.part(i_l) = 1;
        elseif chtot.block_split(i_l) == 3 | chtot.block_split(i_l) == 4
            chtot.part(i_l) = 2;
        else
        end
        
    end
    
    %rename block_split into parts
    chtot.part_n(chtot.part == 0) = cellstr('Calib');
    chtot.part_n(chtot.part == 1) = cellstr('HOC');
    chtot.part_n(chtot.part == 2) = cellstr('Switch');
    
    %reencode diploma
    chtot.diploma_score(:) = 0;
    chtot.diploma_score(strcmp(chtot.last_diploma,'Aucun'),:) = 0;
    chtot.diploma_score(strcmp(chtot.last_diploma,'CAP, BEP'),:) = 1;
    chtot.diploma_score(strcmp(chtot.last_diploma,'BAC'),:) = 2;
    chtot.diploma_score(strcmp(chtot.last_diploma,'BAC + 2'),:) = 3;
    chtot.diploma_score(strcmp(chtot.last_diploma,'BAC + 3'),:) = 4;
    chtot.diploma_score(strcmp(chtot.last_diploma,'BAC + 4'),:) = 5;
    chtot.diploma_score(strcmp(chtot.last_diploma,'BAC + 5'),:) = 6;
    chtot.diploma_score(strcmp(chtot.last_diploma,'BAC + 8'),:) = 7;
    
    %str2double for several
    chtot.chimio = str2double(chtot.chimio);
    chtot.radio = str2double(chtot.radio);
    chtot.FSS_score = str2double(chtot.FSS_score);
    chtot.HAD_Ascore = str2double(chtot.HAD_Ascore);
    chtot.HAD_Dscore = str2double(chtot.HAD_Dscore);
    chtot.STARK_score = str2double(chtot.STARK_score);
    chtot.BIS_score = str2double(chtot.BIS_score);
    chtot.vol_preOP = str2double(chtot.vol_preOP);
    chtot.annoy_calib = str2double(chtot.annoy_calib);
    chtot.annoy_crea = str2double(chtot.annoy_crea);
    chtot.annoy_nswitch = str2double(chtot.annoy_nswitch);
    
    %rename side
    chtot.r_arrow = chtot.choice_side == 37;
    
    % save data for all subjects
    save([savewhere '.mat'],'chtot');  % save as matlab data
    writetable(chtot, [savewhere '.csv']);  % save as csv table
    
    chtot_new = chtot(strcmp(chtot.version_test,'new'),:);
    chtot_old = chtot(strcmp(chtot.version_test,'old'),:);
    save([savewhere '_new.mat'],'chtot_new');  % save as matlab data
    save([savewhere '_old.mat'],'chtot_old');  % save as matlab data
    
    %checking outliers and removing them
    save1 = strcat(cfg.dir.group, cfg.studyname, '_group_choiceREVERSEcalib_data_outliers');  % where to save the final table?
    
    if ~exist([save1 '.txt'], 'file')
        chtot_out = chtot;
        chtot_out(chtot_out.block_split == 0,:) = [];
        
        outliers = varfun(@mean ,chtot_out,'InputVariables',{'choiceSS'},'GroupingVariables',{'subvisit','new_condition'});
        out = outliers.mean_choiceSS > 0.9 | outliers.mean_choiceSS < 0.1; %outliers if more than 90% or less than 10% just after calib
        outliers = outliers(out,:);
        outliers = sortrows(outliers,'mean_choiceSS','ascend');
        
        % write result to text file
        diary([save1 '.txt']);
        disp(outliers);
        diary off;
    end
    
    
    %removing outliers
    for i_o = 1 : length(outliers.subvisit)
        chtot(strcmp(chtot.subvisit,outliers.subvisit{i_o}),:)=[];
    end
    
    
    chtot_new = chtot(strcmp(chtot.version_test,'new'),:);
    chtot_old = chtot(strcmp(chtot.version_test,'old'),:);
    save([save1 '_new.mat'],'chtot_new');  % save as matlab data
    save([save1 '_old.mat'],'chtot_old');  % save as matlab data
    
end

%% CHOICE_SS - PLOT block x patient
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_all_03_subvisit', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 6 ;  % 4 plus legend
    n_y_subplots = 3 ;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data'
    %'_group_choiceALL_data.mat'
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choice&calib_data_outliers_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    
    exclude = varfun(@mean ,chtot,'InputVariables',{'choiceSS'},'GroupingVariables',{'subvisit'});
    exclude = sortrows(exclude,'mean_choiceSS','ascend');
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(exclude);
    diary off;
    
    % gramm - plot
    g = gramm('x', chtot.part_n, 'y', chtot.choiceSS, 'color', chtot.new_condition, 'group', chtot.subvisit);  % define data
    g.stat_summary('geom', {'point', 'line'}, 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.axe_property('YLim', [0 1], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'impulsive choice rate', 'color', 'condition');
    g.set_line_options('styles', {':'})
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    %     % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist


%% CHOICE_SS - PLOT BLOCK CONDITION
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_all_fullexp_04_arrow', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data'
    %'_group_choiceALL_data.mat'
    
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choice&calib_data_outliers_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    
    list_sub = unique(chtot.subvisit);
    ch = [];
    for i_l = 1 : length(list_sub)
        ch_sub = chtot(strcmp(chtot.subvisit, list_sub{i_l}),:);
        if length(unique(ch_sub.block_split)) == 5
            ch = [ch ; ch_sub];
        else
        end
    end
    chtot = ch;
    
    % gramm - plot
    g = gramm('x', chtot.part, 'y', chtot.r_arrow, 'color', chtot.new_condition, 'group', chtot.subvisit);  % define data
    g.stat_summary('type','sem','geom',{'point','line'});% this does plot the mean point
    g.axe_property('YLim',[0 1] ,'XTick', [0:2],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'right arrow pressed', 'color', 'condition');
    g.set_line_options('styles', {':'})
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% CHOICE_SS - PLOT BLOCK CONDITION
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_all_fullexp_05_altern', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data'
    %'_group_choiceALL_data.mat'
    
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choice&calib_data_outliers_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    
    list_sub = unique(chtot.subvisit);
    ch = [];
    for i_l = 1 : length(list_sub)
        ch_sub = chtot(strcmp(chtot.subvisit, list_sub{i_l}),:);
        if length(unique(ch_sub.block_split)) == 5
            ch = [ch ; ch_sub];
        else
        end
    end
    chtot = ch;
    
    % gramm - plot
    g = gramm('x', chtot.part, 'y', chtot.diff_press, 'color', chtot.new_condition, 'group', chtot.subvisit);  % define data
    g.stat_summary('type','sem','geom',{'point','line'});% this does plot the mean point
    g.axe_property('YLim',[0 1] ,'XTick', [0:2],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'diff key pressed', 'color', 'condition');
    g.set_line_options('styles', {':'})
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% CHOICE_SS - PLOT BLOCK PATIENT
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_all_fullexp_06_visit', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data_new'
    %'_group_choiceALL_data.mat'
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choice&calib_data_outliers_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    
    %select only full experiment performers
    list_sub = unique(chtot.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = chtot(strcmp(chtot.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part_n)) == 3
            par = [par ; par_sub];
        else
        end
    end
    chtot = par;
    
    chSS_mean = varfun(@nanmean ,chtot,'InputVariables','choiceSS','GroupingVariables',{'subvisit','visit','type_cond','part','part_n'});
    chSS_mean.choiceSS= chSS_mean.nanmean_choiceSS;
    n_part = unique(chtot.part);
    
    for i_pt = 0 : max(n_part)
        %select the part we want to investigate
        idx_tmp = chSS_mean.part == i_pt;
        ch_tmp = chSS_mean(idx_tmp, :);
        %dispatch according to group
        chSS_grp2 = ch_tmp(ch_tmp.visit==2, :);
        chSS_grp1 = ch_tmp(ch_tmp.visit==1, :);
        %ttest accordingly
        %[H,P,CI,STATS]
        [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.choiceSS, chSS_grp1.choiceSS);
        
        % write result to text file
        diary([savewhere '.txt']);
        disp(strcat('part: ', sprintf('%0.f',i_pt)));
        disp(strcat('mean visit 1 (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.choiceSS))));
        disp(strcat('mean visit 2 (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.choiceSS))));
        
        disp(strcat('V1 vs V2 - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
        
        diary off
        clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
        
    end
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'choiceSS ~ 1 + visit + part_n + part_n:visit ';
    % do the analysis
    chtot = sortrows(chtot,'part_n','ascend'); %put calib part first (as default)
    chtot = sortrows(chtot,'type_cond','ascend'); %put control first (as default)
    res1 = fitglm(chtot, model_formula,'distribution', 'binomial', 'link', 'logit');
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
    diary off;
    
    % gramm - plot
    g = gramm('x', chSS_mean.part_n, 'y', chSS_mean.choiceSS, 'color', chSS_mean.visit);  % define data
    g.stat_summary('type','sem','geom',{'point', 'errorbar'},'dodge',0.4);% this does plot the mean point
    g.set_point_options('base_size',9);
    g.axe_property('YLim', [0.2 1], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'impulsive choice rate', 'color', 'visit');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% CHOICE_SS - PLOT BLOCK CONDITION
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_all_fullexp_07_default', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data'
    %'_group_choiceALL_data.mat'
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choice&calib_data_outliers_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    
    %select only full experiment performers
    list_sub = unique(chtot.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = chtot(strcmp(chtot.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part_n)) == 3
            par = [par ; par_sub];
        else
        end
    end
    chtot = par;
    
    %report mean default choices per subvisit
    save1 = strcat(savewhere, '_list');
    exclude = varfun(@nanmean ,chtot,'InputVariables',{'choice_default'},'GroupingVariables',{'subvisit','new_condition','part','part_n'});
    exclude.choice_default= exclude.nanmean_choice_default;
    exclude = sortrows(exclude,'choice_default','descend');
    % write result to text file
    diary([save1 '.txt']);
    disp(exclude);
    diary off;
    
    
    n_part = unique(exclude.part);
    for i_pt = 0 : max(n_part)
        %select the part we want to investigate
        idx_tmp = exclude.part == i_pt;
        ch_tmp = exclude(idx_tmp, :);
        
        %dispatch according to group
        chSS_grp3 = ch_tmp(strcmp(ch_tmp.new_condition, 'preOP'), :);
        chSS_grp2 = ch_tmp(strcmp(ch_tmp.new_condition, 'postOP'), :);
        chSS_grp1 = ch_tmp(strcmp(ch_tmp.new_condition, 'control'), :);
        
        %ttest accordingly
        %[H,P,CI,STATS]
        [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.choice_default, chSS_grp1.choice_default);
        [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3.choice_default, chSS_grp1.choice_default);
        [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3.choice_default, chSS_grp2.choice_default);
        
        % write result to text file
        diary([savewhere '.txt']);
        disp(strcat('part: ', sprintf('%0.f',i_pt)));
        disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.choice_default))));
        disp(strcat('mean preOP   (n=', sprintf('%0.f',length(unique(chSS_grp3.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp3.choice_default))));
        disp(strcat('mean postOP  (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.choice_default))));
        
        disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
        disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
        disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
        
        diary off
        clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
    end
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'choice_default ~ 1 + visit + part_n + new_condition + part_n:new_condition';        % do the analysis
    chtot = sortrows(chtot,'part_n','ascend'); %put calib part first (as default)
    chtot = sortrows(chtot,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(chtot, model_formula,'distribution', 'binomial', 'link', 'logit');
    % write result to text file
    diary([savewhere '.txt']);
    disp(res1);
    disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
    diary off;
    
    
    % gramm - plot
    g = gramm('x', exclude.part_n, 'y', exclude.choice_default, 'color', exclude.new_condition);  % define data
    g.stat_summary('type','sem','geom',{'point', 'errorbar'},'dodge',0.4);% this does plot the mean point
    g.set_point_options('base_size',9);
    g.axe_property('YLim', [0.2 1], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'proportion default choice', 'color', 'condition');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist


%% CHOICE_SS  - PLOT BLOCK EFFECT
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_all_08_catch_trials', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data'
    %'_group_choiceALL_data.mat'
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceALL_data.mat');  % where is the bic choice data table
    load(filename);  % load choice_data
    
    idx_ct = chtot.trial_type == 99;
    chtot = chtot(idx_ct, :); %keep only the data of the catch trials
    %rename block_split into parts
    chtot.part_n(strcmp(chtot.part,'0'),:) = cellstr('Calib');
    chtot.part_n(strcmp(chtot.part,'1'),:) = cellstr('HOC');
    chtot.part_n(strcmp(chtot.part,'2'),:) = cellstr('Switch');
    
    part_idx = varfun(@nanmean ,chtot,'InputVariables',{'choiceSS'},'GroupingVariables',{'subvisit','new_condition','part','part_n'});
    part_idx.chSScatch = 1 -  part_idx.nanmean_choiceSS;
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(part_idx);
    diary off;
    
    % gramm - plot
    g = gramm('x', part_idx.part_n, 'y', part_idx.chSScatch, 'color', part_idx.new_condition);  % define data
    g.stat_summary('type','sem','geom',{'point', 'errorbar'},'dodge',0.4);% this does plot the mean point
    g.set_point_options('base_size',9);
    g.axe_property('YLim', [0 0.5], 'XTick', [1:3],'XTickLabel', {'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', ' proportion misssed catch trials', 'color', 'condition');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% MAIN EFFECT OF SESSION IN PARAMETERS OF EXPONENTIAL + BIAS WITHIN SUBJECT
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_all_fullexp_09_subj_follow', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 5 plus legend
    n_y_subplots = 3;
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceREVERSEcalib_data_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    chtot.task(strcmp(chtot.task,'TD_ID'),:)={'IvD'};
    chtot.task(strcmp(chtot.task,'TD_DD'),:)={'DvD'};
    
    
    ch_vis = [];
    list = unique(chtot.subject_id);
    
    for i_s = 1 : length(list)
        ch_sub = chtot(strcmp(chtot.subject_id, list{i_s}), :);
        
        savesub = sprintf('%s%s_choice_tasks_all_09_follow_%s', cfg.dir.group_subject, cfg.studyname,list{i_s});   % construct file name
        if ~exist([savesub '.png'], 'file')  % only run, when output does not exist
            % open figure
            f = figure(1); clf;
            set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
                'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
                'PaperUnits','centimeters', 'Visible', 'off');
            
            % gramm - plot
            g = gramm('x', ch_sub.part_n, 'y', ch_sub.choiceSS, 'color', ch_sub.new_condition, 'line', ch_sub.date, 'marker',ch_sub.task);  % define data
            g.stat_summary('geom', {'point', 'line'}, 'setylim', false, 'type', 'sem');  % this does plot the mean point
            g.set_title(strcat(list{i_s},' - visit(s) : ', sprintf('%0.f',length(unique(ch_sub.subvisit)))));
            g.set_point_options('base_size',9);
            g.set_line_options('styles', {':','-','-.','--'})
            g.axe_property('YLim',[0 1] ,'XLim', [0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
            g.set_names('column', '', 'x', ' ', 'y', 'impulsive choice rate', 'color', 'group','line','date(YYYYMMDD)','marker','Task');
            g.set_text_options('base_size',17);
            g.set_color_options('map',cfg.fig.colormap);
            g.draw();
            
            % save figure
            export_fig([savesub '.png'], '-r600', '-nocrop');  % save figure
            close(f);  % close figure
        end
        if length(unique(ch_sub.condition)) > 1
            ch_sub(ch_sub.visit > 2,:) = []; %remove visit 3
            if length(unique(ch_sub(strcmp(ch_sub.condition,'preOP'),:).part_n)) == length(unique(ch_sub(strcmp(ch_sub.condition,'postOP'),:).part_n))
                ch_vis = [ch_vis ; ch_sub];
            else
            end
        else
        end
        
    end
    
    chSS_mean = varfun(@nanmean ,ch_vis,'InputVariables','choiceSS','GroupingVariables',{'subvisit','visit','condition','part','part_n'});
    chSS_mean.choiceSS= chSS_mean.nanmean_choiceSS;
    n_part = unique(chSS_mean.part);
    ch_prepost = [];
    for i_pt = 0 : max(n_part)
        %select the part we want to investigate
        idx_tmp = chSS_mean.part == i_pt;
        ch_tmp = chSS_mean(idx_tmp, :);
        %dispatch according to group
        chSS_grp1 = ch_tmp(strcmp(ch_tmp.condition, 'preOP'), :);
        chSS_grp2 = ch_tmp(strcmp(ch_tmp.condition, 'postOP'), :);
        ch_cmp = chSS_grp1(:,1:5);
        ch_cmp.chSSpre = chSS_grp1.choiceSS;
        ch_cmp.chSSpost = chSS_grp2.choiceSS;
        ch_prepost = [ch_prepost ; ch_cmp];
        %ttest accordingly
        %[H,P,CI,STATS]
        [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest(chSS_grp2.choiceSS, chSS_grp1.choiceSS);
        
        % write result to text file
        diary([savewhere '.txt']);
        disp(strcat('part: ', sprintf('%0.f',i_pt)));
        disp(strcat('mean preOP  (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.choiceSS))));
        disp(strcat('mean postOP (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.choiceSS))));
        
        disp(strcat('preOP vs postOP - repeated ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
        
        diary off
        clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
        
    end
    
    % statistics -generalized linear regression model for repeated measures
    % set up formula
    ch_prepost.condition = [];
    prepost = table([1 2]','VariableNames',{'chSS'});
    rm =fitrm(ch_prepost,'chSSpre-chSSpost~part_n','WithinDesign',prepost);
    ranovatbl = ranova(rm);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITrepM')
    disp(rm.Coefficients);
    disp(rm.Covariance);
    disp(ranovatbl);
    diary off;
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'choiceSS ~ 1 + condition + part_n + part_n:condition ';
    % do the analysis
    ch_vis = sortrows(ch_vis,'part_n','ascend'); %put calib part first (as default)
    ch_vis = sortrows(ch_vis,'type_cond','ascend'); %put control first (as default)
    res1 = fitglm(ch_vis, model_formula,'distribution', 'binomial', 'link', 'logit');
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
    diary off;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    
    % gramm - plot
    g = gramm('x', chSS_mean.part_n, 'y', chSS_mean.choiceSS, 'color', chSS_mean.condition);  % define data
    g.stat_summary('type','sem','geom',{'point', 'errorbar'},'dodge',0.4);% this does plot the mean point
    g.set_point_options('base_size',9);
    g.axe_property('YLim', [0.2 1], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'impulsive choice rate', 'color', 'visit');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end

%% MAIN EFFECT OF SESSION IN PARAMETERS OF EXPONENTIAL + BIAS WITHIN SUBJECT
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_all_fullexp_09_cont_follow', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 5 plus legend
    n_y_subplots = 3;
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceREVERSEcalib_data_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    chtot.task(strcmp(chtot.task,'TD_ID'),:)={'IvD'};
    chtot.task(strcmp(chtot.task,'TD_DD'),:)={'DvD'};
    
    chtot = chtot(strcmp(chtot.condition,'control'),:);
    ch_vis = [];
    list = unique(chtot.subject_id);
    
    for i_s = 1 : length(list)
        ch_sub = chtot(strcmp(chtot.subject_id, list{i_s}), :);
        if length(unique(ch_sub.visit)) > 1
            ch_sub(ch_sub.visit > 2,:) = []; %remove visit 3
            if length(unique(ch_sub(strcmp(ch_sub.condition,'preOP'),:).part_n)) == length(unique(ch_sub(strcmp(ch_sub.condition,'postOP'),:).part_n))
                ch_vis = [ch_vis ; ch_sub];
            else
            end
        else
        end
        
    end
    
    chSS_mean = varfun(@nanmean ,ch_vis,'InputVariables','choiceSS','GroupingVariables',{'subvisit','visit','part','part_n'});
    chSS_mean.choiceSS= chSS_mean.nanmean_choiceSS;
    n_part = unique(chSS_mean.part);
    ch_prepost = [];
    for i_pt = 0 : max(n_part)
        %select the part we want to investigate
        idx_tmp = chSS_mean.part == i_pt;
        ch_tmp = chSS_mean(idx_tmp, :);
        %dispatch according to group
        chSS_grp1 = ch_tmp((ch_tmp.visit ==1), :);
        chSS_grp2 = ch_tmp((ch_tmp.visit ==2), :);
        ch_cmp = chSS_grp1(:,1:4);
        ch_cmp.chSSpre = chSS_grp1.choiceSS;
        ch_cmp.chSSpost = chSS_grp2.choiceSS;
        ch_prepost = [ch_prepost ; ch_cmp];
        %ttest accordingly
        %[H,P,CI,STATS]
        [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest(chSS_grp2.choiceSS, chSS_grp1.choiceSS);
        
        % write result to text file
        diary([savewhere '.txt']);
        disp(strcat('part: ', sprintf('%0.f',i_pt)));
        disp(strcat('mean visit 1  (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.choiceSS))));
        disp(strcat('mean visit 2 (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.choiceSS))));
        
        disp(strcat('test vs retest - repeated ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
        
        diary off
        clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
        
    end
    
    % statistics -generalized linear regression model for repeated measures
    % set up formula
    ch_prepost.visit = [];
    prepost = table([1 2]','VariableNames',{'chSS'});
    rm =fitrm(ch_prepost,'chSSpre-chSSpost~part_n','WithinDesign',prepost);
    ranovatbl = ranova(rm);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITrepM')
    disp(rm.Coefficients);
    disp(rm.Covariance);
    disp(ranovatbl);
    diary off;
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'choiceSS ~ 1 + visit + part_n + part_n:visit ';
    % do the analysis
    ch_vis = sortrows(ch_vis,'part_n','ascend'); %put calib part first (as default)
    ch_vis = sortrows(ch_vis,'type_cond','ascend'); %put control first (as default)
    res1 = fitglm(ch_vis, model_formula,'distribution', 'binomial', 'link', 'logit');
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
    diary off;
    
    chSS_mean.test(:) = {'retest'};
    chSS_mean.test(chSS_mean.visit == 1,:) = {'test'};
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    
    % gramm - plot
    g = gramm('x', chSS_mean.part_n, 'y', chSS_mean.choiceSS, 'color', chSS_mean.test);  % define data
    g.stat_summary('type','sem','geom',{'point', 'errorbar'},'dodge',0.4);% this does plot the mean point
    g.set_point_options('base_size',9);
    g.axe_property('YLim', [0.2 1], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'impulsive choice rate', 'color', 'visit');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end

%% CHOICE_SS - PLOT BLOCK CONDITION
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_all_10_condition_comparison', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data'
    %'_group_choiceALL_data.mat'
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choice&calib_data_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    
    
    chtot.outliers(:) = 0;
    chtot.quitters(:) = 1;
    
    %indexing outliers
    filename = strcat(cfg.dir.group, cfg.studyname,  '_outliers_list.mat');
    load(filename); % this load choice data
    for i_o = 1 : length(outliers.subvisit)
        chtot.outliers(strcmp(chtot.subvisit,outliers.subvisit{i_o}),:)=1;
    end
    
    %indexing quitters
    idx_nq = chtot.part > 1;
    idx_nq = unique(chtot.subvisit(idx_nq,:));
    for i_q = 1 : length(idx_nq)
        chtot.quitters(strcmp(chtot.subvisit,idx_nq{i_q}),:)=0;
    end
    
    
    out_quit =  varfun(@mean,chtot,'InputVariables',{'choiceSS'},'GroupingVariables',{'subvisit','visit','condition','outliers','quitters'});
    save([savewhere 'outquit_new.mat'],'out_quit');
    
    %open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', chtot.part, 'y', chtot.choiceSS, 'color', chtot.new_condition, 'lightness', chtot.quitters, 'marker', chtot.outliers);  % define data
    g.stat_summary('geom', {'point'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.axe_property('XLim', [-0.5 2.5], 'XTick', [0:2],'XTickLabel', {'Calib', 'HOC', 'Switch'},'YLim',[0 1]);
    g.set_color_options('chroma',200,'lightness_range',[60 90],'chroma_range',[60 30]);
    g.set_point_options('base_size',8);
    g.set_text_options('base_size',17);
    g.set_names('column', '', 'x', 'run', 'y', 'impulsive choice rate', 'color', 'group', 'lightness', 'quitters', 'marker', 'outliers');
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end % if files does not exist



%% CHOICE_SS - PLOT BLOCK PATIENT
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_cleaned_fullexp_11_indiv', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data_new'
    %'_group_choiceALL_data.mat'
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceREVERSEcalib_data_outliers_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    chtot.task(strcmp(chtot.task,'TD_ID'),:)={'IvD'};
    chtot.task(strcmp(chtot.task,'TD_DD'),:)={'DvD'};
    
    %keep only first visits
    chtot(chtot.visit >1,:)=[];
    
    %exclusion criteria
    chtot = chtot(strcmp(chtot.exclusion_reason,'RAS'),:);
    
    
    %select only full experiment performers
    list_sub = unique(chtot.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = chtot(strcmp(chtot.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part_n)) == 3
            par = [par ; par_sub];
        else
        end
    end
    chtot = par;
    
    
    task = unique(chtot.task);
    for i_t = 1:length(task)
        ch_t = chtot(strcmp(chtot.task,task(i_t)),:);
        
        diary([savewhere '.txt']);
        disp(strcat('task: ',task{i_t}));
        diary off
        
        n_part = unique(ch_t.part);
        for i_pt = 0 : max(n_part)
            %select the part we want to investigate
            idx_tmp = ch_t.part == i_pt;
            ch_tmp = ch_t(idx_tmp, :);
            %dispatch according to group
            chSS_grp2 = ch_tmp(strcmp(ch_tmp.type_cond, 'patient'), :);
            chSS_grp1 = ch_tmp(strcmp(ch_tmp.type_cond, 'control'), :);
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.choiceSS, chSS_grp1.choiceSS);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.choiceSS))));
            disp(strcat('mean patient (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.choiceSS))));
            
            disp(strcat('patient vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            
            diary off
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
            
            
            
        end
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'choiceSS ~ 1 + part_n + type_cond + part_n:type_cond ';
        % do the analysis
        ch_t = sortrows(ch_t,'part_n','ascend'); %put calib part first (as default)
        ch_t = sortrows(ch_t,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(ch_t, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
    end
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'choiceSS ~ 1 + part_n + type_cond + part_n:type_cond ';
    % do the analysis
    chtot = sortrows(chtot,'part_n','ascend'); %put calib part first (as default)
    chtot = sortrows(chtot,'type_cond','ascend'); %put control first (as default)
    res1 = fitglm(chtot, model_formula);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM MIXED TASKS')
    disp(res1);
    disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
    diary off;
    
    chtot(strcmp(chtot.part_n,'crea'),:)=[];
    
    descrdata.participant = varfun(@nanmean ,chtot,'InputVariables','choiceSS','GroupingVariables',{'subvisit','type_cond','part_n'});
    descrdata.mean = varfun(@nanmean ,descrdata.participant,'InputVariables','nanmean_choiceSS','GroupingVariables',{'type_cond','part_n'});
    descrdata.std = varfun(@nanstd ,descrdata.participant,'InputVariables','nanmean_choiceSS','GroupingVariables',{'type_cond','part_n'});
    descrdata.min = varfun(@min ,descrdata.participant,'InputVariables','nanmean_choiceSS','GroupingVariables',{'type_cond','part_n'});
    descrdata.max = varfun(@max ,descrdata.participant,'InputVariables','nanmean_choiceSS','GroupingVariables',{'type_cond','part_n'});
    
    %calculate delta from calibration for each parameters
    descrdata.participant.delta_chSS(:) = nan;
    list_sub = unique(descrdata.participant.subvisit);
    for i_s = 1 : length(list_sub)
        idx_s =   strcmp(descrdata.participant.subvisit, list_sub{i_s});
        calib = descrdata.participant(idx_s,:).nanmean_choiceSS(strcmp(descrdata.participant(idx_s,:).part_n,'Calib'));
        if length(calib) == 1
            descrdata.participant.delta_chSS(idx_s,:) = (descrdata.participant.nanmean_choiceSS(idx_s,:) - calib)*100 ;
        else
        end
    end
    descrdata.participant.delta_10 = abs(descrdata.participant.delta_chSS) > 10;
    list=descrdata.participant.subvisit(descrdata.participant.delta_10);
    for i_s = 1 : length(list)
        descrdata.participant.delta_group(strcmp(descrdata.participant.subvisit,list(i_s)),:)=1;
    end
    
    descrdata.participant(strcmp(descrdata.participant.part_n,'HOC'),:)=[];
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm - plot
    g = gramm('x', descrdata.participant.part_n, 'y', descrdata.participant.nanmean_choiceSS*100, 'color', ...
        descrdata.participant.type_cond, 'group', descrdata.participant.subvisit);  % define data
    g.stat_summary('type','sem','geom',{'point','line'},'dodge',0.4);% this does plot the mean point
    g.set_point_options('base_size',5);
    g.set_line_options('styles', {':'})
    g.axe_property('XLim',[0.75 2.25],'XTick', [1 2],'XTickLabel', {'Calib', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'impulsive choice rate', 'color', 'group', 'linestyle', '+/- 10% variation');
    g.set_text_options('base_size',17);
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end % if files does not exist

%% CHOICE_SS - PLOT BLOCK PATIENT
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_cleaned_fullexp_12_delta', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data_new'
    %'_group_choiceALL_data.mat'
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceREVERSEcalib_data_outliers_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    chtot.task(strcmp(chtot.task,'TD_ID'),:)={'IvD'};
    chtot.task(strcmp(chtot.task,'TD_DD'),:)={'DvD'};
    
    %keep only first visits
    chtot(chtot.visit >1,:)=[];
    
    %exclusion criteria
    chtot = chtot(strcmp(chtot.exclusion_reason,'RAS'),:);
    
    %select only full experiment performers
    list_sub = unique(chtot.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = chtot(strcmp(chtot.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part_n)) == 3
            par = [par ; par_sub];
        else
        end
    end
    chtot = par;
    
    chtot(strcmp(chtot.part_n,'crea'),:)=[];
    
    descrdata.participant = varfun(@nanmean ,chtot,'InputVariables','choiceSS','GroupingVariables',{'subvisit','type_cond','part_n'});
    descrdata.mean = varfun(@nanmean ,descrdata.participant,'InputVariables','nanmean_choiceSS','GroupingVariables',{'type_cond','part_n'});
    descrdata.std = varfun(@nanstd ,descrdata.participant,'InputVariables','nanmean_choiceSS','GroupingVariables',{'type_cond','part_n'});
    descrdata.min = varfun(@min ,descrdata.participant,'InputVariables','nanmean_choiceSS','GroupingVariables',{'type_cond','part_n'});
    descrdata.max = varfun(@max ,descrdata.participant,'InputVariables','nanmean_choiceSS','GroupingVariables',{'type_cond','part_n'});
    
    %calculate delta from calibration for each parameters
    descrdata.participant.delta_chSS(:) = nan;
    list_sub = unique(descrdata.participant.subvisit);
    for i_s = 1 : length(list_sub)
        idx_s =   strcmp(descrdata.participant.subvisit, list_sub{i_s});
        calib = descrdata.participant(idx_s,:).nanmean_choiceSS(strcmp(descrdata.participant(idx_s,:).part_n,'Calib'));
        if length(calib) == 1
            descrdata.participant.delta_chSS(idx_s,:) = (descrdata.participant.nanmean_choiceSS(idx_s,:) - calib)*100 ;
        else
        end
    end
    
    descrdata.participant(strcmp(descrdata.participant.part_n,'HOC'),:)=[];
    descrdata.participant(strcmp(descrdata.participant.part_n,'Calib'),:)=[];
    
    
    %dispatch according to group
    chSS_grp2 = descrdata.participant(strcmp(descrdata.participant.type_cond, 'patient'), :);
    chSS_grp1 = descrdata.participant(strcmp(descrdata.participant.type_cond, 'control'), :);
    %ttest accordingly
    [h_hnull_1,p_hnull_1] = ttest(chSS_grp1.delta_chSS);
    [h_hnull_2,p_hnull_2]  = ttest(chSS_grp2.delta_chSS);
    %comparison group
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.delta_chSS, chSS_grp1.delta_chSS);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2.delta_chSS, chSS_grp1.delta_chSS,'Tail','right');
    [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2.delta_chSS,chSS_grp1.delta_chSS);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp('ONLY SWITCH DELTA');
    disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1.delta_chSS)),') = ',sprintf('%0.3f',mean(chSS_grp1.delta_chSS))));
    disp(strcat('mean patient (n=', sprintf('%0.f',length(chSS_grp2.delta_chSS)),') = ',sprintf('%0.3f',mean(chSS_grp2.delta_chSS))));
    
    disp(strcat('delta patient against null - ttest : h=',sprintf('%0.0f',h_hnull_2),', p = ' ,sprintf('%0.3f',p_hnull_2)));
    disp(strcat('delta control against null - ttest : h=',sprintf('%0.0f',h_hnull_1),', p = ' ,sprintf('%0.3f',p_hnull_1)));
    
    
    disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
    diary off
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    set(gca,'FontSize',15)
    % gramm - plot
    g = violinplot(descrdata.participant.delta_chSS,descrdata.participant.type_cond,'ShowMean',true);  % define data
    g(1, 1).ViolinColor = [ 0    0.4470    0.7410];
    g(1, 2).ViolinColor = [ 0.8500    0.3250    0.0980];
    g(1, 1).ScatterPlot.MarkerFaceAlpha = 1;
    g(1, 2).ScatterPlot.MarkerFaceAlpha = 1;
    ylabel('impulsive choice rate delta');
    xlabel('(Switch - Calib) run');
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end % if files does not exist


%% CHOICE_SS - PLOT BLOCK PATIENT
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_cleaned_fullexp_13_RT', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data_new'
    %'_group_choiceALL_data.mat'
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceREVERSEcalib_data_outliers_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    chtot.task(strcmp(chtot.task,'TD_ID'),:)={'IvD'};
    chtot.task(strcmp(chtot.task,'TD_DD'),:)={'DvD'};
    
    %keep only first visits
    chtot(chtot.visit >1,:)=[];
    
    %exclusion criteria
    chtot = chtot(strcmp(chtot.exclusion_reason,'RAS'),:);
    
    
    %select only full experiment performers
    list_sub = unique(chtot.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = chtot(strcmp(chtot.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part_n)) == 3
            par = [par ; par_sub];
        else
        end
    end
    chtot = par;
    
    RT_mean = varfun(@nanmean ,chtot,'InputVariables','raw_response_time','GroupingVariables',{'subvisit','type_cond','part','part_n'});
    
    n_part = unique(RT_mean.part);
    for i_pt = 0 : max(n_part)
        %select the part we want to investigate
        idx_tmp = RT_mean.part == i_pt;
        ch_tmp = RT_mean(idx_tmp, :);
        %dispatch according to group
        chSS_grp2 = ch_tmp(strcmp(ch_tmp.type_cond, 'patient'), :);
        chSS_grp1 = ch_tmp(strcmp(ch_tmp.type_cond, 'control'), :);
        %ttest accordingly
        %[H,P,CI,STATS]
        [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.nanmean_raw_response_time, chSS_grp1.nanmean_raw_response_time);
        
        % write result to text file
        diary([savewhere '.txt']);
        disp(strcat('part: ', sprintf('%0.f',i_pt)));
        disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.nanmean_raw_response_time))));
        disp(strcat('mean patient (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.nanmean_raw_response_time))));
        
        disp(strcat('patient vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
        
        diary off
        clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
        
        
    end
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'nanmean_raw_response_time ~ 1 + part_n + type_cond + part_n:type_cond ';
    model_formula2 = 'nanmean_raw_response_time ~ 1 + part + type_cond + part:type_cond ';
    % do the analysis
    RT_mean = sortrows(RT_mean,'part_n','ascend'); %put calib part first (as default)
    RT_mean = sortrows(RT_mean,'type_cond','ascend'); %put control first (as default)
    res1 = fitglm(RT_mean, model_formula);
    res2 = fitglm(RT_mean, model_formula2);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM MIXED TASKS')
    disp('BLOCK AS NOMINAL VARIABLE')
    disp(res1);
    disp('BLOCK AS CONTINUOUS VARIABLE')
    disp(res2);
    disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
    diary off;
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm - plot
    g = gramm('x', RT_mean.part, 'y', RT_mean.nanmean_raw_response_time, 'color', RT_mean.type_cond);  % define data
    g.stat_summary('type','sem','geom',{'line' 'point', 'errorbar'},'dodge',0.4);% this does plot the mean point
    g.set_point_options('base_size',9);
    g.set_line_options('styles', {':'})
    g.axe_property('XLim',[-0.5 2.5] ,'YLim',[1 5] ,'XTick', [0:2],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'response time (in s)', 'color', 'group', 'marker', 'task');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end % if files does not exist

%% CHOICE_SS - PLOT BLOCK PATIENT
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_cleaned_fullexp_14_chSS', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data_new'
    %'_group_choiceALL_data.mat'
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceREVERSEcalib_data_outliers_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    chtot.task(strcmp(chtot.task,'TD_ID'),:)={'IvD'};
    chtot.task(strcmp(chtot.task,'TD_DD'),:)={'DvD'};
    
    %keep only first visits
    chtot(chtot.visit >1,:)=[];
    
    %exclusion criteria
    chtot = chtot(strcmp(chtot.exclusion_reason,'RAS'),:);
    
    
    %select only full experiment performers
    list_sub = unique(chtot.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = chtot(strcmp(chtot.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part_n)) == 3
            par = [par ; par_sub];
        else
        end
    end
    chtot = par;
    
    chSS_mean = varfun(@nanmean ,chtot,'InputVariables','choiceSS','GroupingVariables',{'subvisit','type_cond','part','part_n'});
    
    n_part = unique(chSS_mean.part);
    for i_pt = 0 : max(n_part)
        %select the part we want to investigate
        idx_tmp = chSS_mean.part == i_pt;
        ch_tmp = chSS_mean(idx_tmp, :);
        %dispatch according to group
        chSS_grp2 = ch_tmp(strcmp(ch_tmp.type_cond, 'patient'), :);
        chSS_grp1 = ch_tmp(strcmp(ch_tmp.type_cond, 'control'), :);
        %ttest accordingly
        %[H,P,CI,STATS]
        [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.nanmean_choiceSS, chSS_grp1.nanmean_choiceSS);
        
        % write result to text file
        diary([savewhere '.txt']);
        disp(strcat('part: ', sprintf('%0.f',i_pt)));
        disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.nanmean_choiceSS))));
        disp(strcat('mean patient (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.nanmean_choiceSS))));
        
        disp(strcat('patient vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
        
        diary off
        clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
        
        
    end
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'choiceSS ~ 1 + part_n + type_cond + part_n:type_cond ';
    model_formula2 = 'choiceSS ~ 1 + part + type_cond + part:type_cond ';
    % do the analysis
    chtot = sortrows(chtot,'part_n','ascend'); %put calib part first (as default)
    chtot = sortrows(chtot,'type_cond','ascend'); %put control first (as default)
    res1 = fitglm(chtot, model_formula ,'distribution', 'binomial', 'link', 'logit');
    res2 = fitglm(chtot, model_formula2,'distribution', 'binomial', 'link', 'logit');
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM MIXED TASKS')
    disp('BLOCK AS NOMINAL VARIABLE')
    disp(res1);
    disp('BLOCK AS CONTINUOUS VARIABLE')
    disp(res2);
    disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
    diary off;
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    
    % gramm - plot
    g = gramm('x', chSS_mean.part, 'y', chSS_mean.nanmean_choiceSS*100, 'color', chSS_mean.type_cond);  % define data
    g.stat_summary('type','sem','geom',{'line' 'point', 'errorbar'},'dodge',0.4);% this does plot the mean point
    g.set_point_options('base_size',9);
    g.set_line_options('styles', {':'})
    g.axe_property('XLim',[-0.5 2.5] ,'YLim',[40 65] ,'XTick', [0:2],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'impulsive choice rate', 'color', 'group', 'marker', 'task');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end % if files does not exist

%% CHOICE_SS - PLOT BLOCK PATIENT
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_cleaned_fullexp_15_DD_ID', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data_new'
    %'_group_choiceALL_data.mat'
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceREVERSEcalib_data_outliers_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    chtot.task(strcmp(chtot.task,'TD_ID'),:)={'IvD'};
    chtot.task(strcmp(chtot.task,'TD_DD'),:)={'DvD'};
    
    %exclusion criteria
    chtot = chtot(strcmp(chtot.exclusion_reason,'RAS'),:);
    
    %keep only first visits
    chtot(chtot.visit >1,:)=[];
    
    %select only full experiment performers
    list_sub = unique(chtot.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = chtot(strcmp(chtot.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part_n)) == 3
            par = [par ; par_sub];
        else
        end
    end
    chtot = par;
    
    
    chtot=chtot(strcmp(chtot.part_n,'Switch'),:);
    chSS_mean = varfun(@nanmean ,chtot,'InputVariables','choiceSS','GroupingVariables',{'subvisit','type_cond','task'});
    chSS_mean.choiceSS = chSS_mean.nanmean_choiceSS;
    
    task = unique(chSS_mean.task);
    for i_t = 1:length(task)
        ch_t = chSS_mean(strcmp(chSS_mean.task,task(i_t)),:);
        
        diary([savewhere '.txt']);
        disp('ONLY SWITCH TASK');
        disp(strcat('task: ',task{i_t}));
        diary off
        
        ch_tmp = ch_t;
        %dispatch according to group
        chSS_grp2 = ch_tmp(strcmp(ch_tmp.type_cond, 'patient'), :);
        chSS_grp1 = ch_tmp(strcmp(ch_tmp.type_cond, 'control'), :);
        %ttest accordingly
        %[H,P,CI,STATS]
        [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.choiceSS, chSS_grp1.choiceSS);
        
        % write result to text file
        diary([savewhere '.txt']);
        disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.choiceSS))));
        disp(strcat('mean patient (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.choiceSS))));
        
        disp(strcat('patient vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
        
        diary off
        clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
        
    end
    
    %AVERAGE CHOICES PER PARTICIPANTS
    ch_one = varfun(@mean ,chtot,'InputVariables','choiceSS','GroupingVariables',{'subvisit'});
    ch_combined = [];
    list_sub = unique(ch_one.subvisit);
    for i_ls = 1:length(list_sub)
        ch_s = chtot(strcmp(chtot.subvisit,list_sub(i_ls)),:);
        ch_tmp = [ch_one(strcmp(ch_one.subvisit,list_sub(i_ls)),:), ch_s(1,17:19),ch_s(1,21:60)];
        ch_combined = [ch_combined;ch_tmp];
    end
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula_lesion = 'mean_choiceSS ~ 1 +  LocDG + LocFront + vol_preOP';
    model_formula_psysoc = 'mean_choiceSS ~ 1 + age + sex + diploma_score + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
    model_formula_treat  = 'mean_choiceSS ~ 1 + new_condition + radio + chimio + tt_antiep';
    model_formula_orthoass  = 'mean_choiceSS ~ 1 + flu_anx + flu_P + TMT_A + TMT_B + TMT_flex + span_for + span_back ';
    model_formula_treat2  = 'mean_choiceSS ~ 1 + days_since_OP + radio + chimio + tt_antiep';
    
    % do the analysis
    ch_combined.days_since_OP(isnan(ch_combined.days_since_OP)) = 0;
    par_pat = ch_combined(strcmp(ch_combined.type_cond,'patient'),:);
    par_ctrl = ch_combined(strcmp(ch_combined.type_cond,'control'),:);
    par_pat = sortrows(par_pat,'part_n','ascend'); %put calib part first (as default)
    par_pat = sortrows(par_pat,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(par_pat, model_formula_lesion);
    res2 = fitglm(par_pat, model_formula_psysoc);
    res2ctrl = fitglm(par_ctrl, model_formula_psysoc);
    res3 = fitglm(par_pat, model_formula_treat);
    res3b = fitglm(par_pat, model_formula_treat2);
    res4 = fitglm(par_pat, model_formula_orthoass);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp('LESIONAL FACTOR')
    disp(res1);
    disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
    disp('PSYCHOSOCIOLOGICAL FACTOR')
    disp('PATIENTS')
    disp(res2);
    disp(strcat('R2 =',sprintf('%0.3f',res2.Rsquared.Ordinary)));
    disp('CONTROL')
    disp(res2ctrl);
    disp(strcat('R2 =',sprintf('%0.3f',res2ctrl.Rsquared.Ordinary)));
    disp('TREATMENT FACTOR')
    disp(res3);
    disp(strcat('R2 =',sprintf('%0.3f',res3.Rsquared.Ordinary)));
    disp(res3b);
    disp(strcat('R2 =',sprintf('%0.3f',res3b.Rsquared.Ordinary)));
    disp('ORTHO ASSESSMENT FACTOR')
    disp(res4);
    disp(strcat('R2 =',sprintf('%0.3f',res4.Rsquared.Ordinary)));
    diary off;
    
    
    chSS_mean = varfun(@nanmean ,chtot,'InputVariables','choiceSS','GroupingVariables',{'subvisit','type_cond','part_n','task'});
    
    % define size of figure by number of sub plots
    n_x_subplots = 3;  % 4 plus legend
    n_y_subplots = 2.5;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm - plot
    g = gramm('x', chSS_mean.part_n, 'y', chSS_mean.nanmean_choiceSS*100, 'color', chSS_mean.type_cond, 'marker',chSS_mean.task);  % define data
    g.stat_summary('type','sem','geom',{'point', 'errorbar'},'dodge',0.4);% this does plot the mean point
    g.set_point_options('base_size',9);
    g.axe_property('YLim',[40 65] );
    g.set_names('column', '', 'row', '', 'x', '', 'y', 'impulsive choice rate', 'color', 'group', 'marker', 'task');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.results.stat_summary(1).point_handle.XData = 0.8;
    g.results.stat_summary(1).errorbar_handle.XData = [0.8 0.8 NaN 0.78 0.82 NaN 0.78 0.82 NaN];
    g.results.stat_summary(2).point_handle.XData = 1.2;
    g.results.stat_summary(2).errorbar_handle.XData = [1.2 1.2 NaN 1.18 1.22 NaN 1.18 1.22 NaN];
    g.results.stat_summary(3).point_handle.XData = 0.9;
    g.results.stat_summary(3).errorbar_handle.XData = [0.9 0.9 NaN 0.88 0.92 NaN 0.88 0.92 NaN];
    g.results.stat_summary(4).point_handle.XData = 1.3;
    g.results.stat_summary(4).errorbar_handle.XData = [1.3 1.3 NaN 1.28 1.32 NaN 1.28 1.32 NaN];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist



%% CHOICE_SS - PLOT BLOCK PATIENT
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_old_vs_new_fullexp_16_chSS', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data_new'
    %'_group_choiceALL_data.mat'
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choice&calib_data');  % where is the bic choice data table
    load(filename);  % load choice_data
    %     chtot(strcmp(chtot.type_cond,'control'),:)=[];
    %     chtot = chtot(chtot.visit == 1,:);
    chtot_old = chtot(strcmp(chtot.version_test,'old'),:);
    
    %         filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceREVERSEcalib_data_old');  % where is the bic choice data table
    %         load(filename);  % load choice_data
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceREVERSEcalib_data_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    %         chtot_old.part_n(strcmp(chtot_old.part_n,'Calib'),:)={'HOC'};
    %         chtot_old.part(chtot_old.part == 0,:)= 1;
    chtot_new = [chtot_new(:,1:53) chtot_new(:,57:60)];
    chtot = [chtot_old(:,1:57);chtot_new];
    chtot(strcmp(chtot.type_cond,'control'),:)=[];
    chtot = chtot(chtot.visit == 1,:);
    
    %select only full experiment performers
    list_sub = unique(chtot.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = chtot(strcmp(chtot.subvisit, list_sub{i_l}),:);
        if strcmp(unique(par_sub.version_test),'new')
            if length(unique(par_sub.part_n)) == 3
                par = [par ; par_sub];
            else
            end
        else
            if length(unique(par_sub.part_n)) == 2
                par = [par ; par_sub];
            else
            end
        end
    end
    chtot = par;
    
    
    chtot.task(strcmp(chtot.task,'TD_ID'),:)={'IvD'};
    chtot.task(strcmp(chtot.task,'TD_DD'),:)={'DvD'};
    %rename block_split into parts
    chtot.part_n(strcmp(chtot.part,'0'),:) = cellstr('Calib');
    chtot.part_n(strcmp(chtot.part,'1'),:) = cellstr('HOC');
    chtot.part_n(strcmp(chtot.part,'2'),:) = cellstr('Switch');
    
    chSS_mean = varfun(@nanmean ,chtot,'InputVariables','choiceSS','GroupingVariables',{'version_test','subvisit','part','part_n'});
    
    n_part = unique(chSS_mean.part);
    for i_pt = 1 : max(n_part)
        %select the part we want to investigate
        idx_tmp = chSS_mean.part == i_pt;
        ch_tmp = chSS_mean(idx_tmp, :);
        %dispatch according to group
        chSS_grp2 = ch_tmp(strcmp(ch_tmp.version_test, 'new'), :);
        chSS_grp1 = ch_tmp(strcmp(ch_tmp.version_test, 'old'), :);
        %ttest accordingly
        %[H,P,CI,STATS]
        [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.nanmean_choiceSS, chSS_grp1.nanmean_choiceSS);
        
        % write result to text file
        diary([savewhere '.txt']);
        disp(strcat('part: ', sprintf('%0.f',i_pt)));
        disp(strcat('mean old (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.nanmean_choiceSS))));
        disp(strcat('mean new (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.nanmean_choiceSS))));
        
        disp(strcat('old vs new protocol - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
        
        diary off
        clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
        
        
    end
    
    
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    
    % gramm - plot
    g = gramm('x', chSS_mean.part, 'y', chSS_mean.nanmean_choiceSS*100, 'color', chSS_mean.version_test);  % define data
    g.stat_summary('type','sem','geom',{'line' 'point', 'errorbar'},'dodge',0.4);% this does plot the mean point
    g.set_point_options('base_size',9);
    g.set_line_options('styles', {':'})
    g.axe_property('XLim',[-0.5 2.5] ,'YLim',[25 70] ,'XTick', [0:2],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'impulsive choice rate', 'color', 'protocol', 'marker', 'task');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end % if files does not exist

%% CHOICE_SS - PLOT BLOCK PATIENT
% =========================================================================
savewhere = sprintf('%s%s_choice_tasks_cleaned_frontal_17_chSS', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    %Choose your file :
    %'_group_choice&calib_data_new'
    %'_group_choiceALL_data.mat'
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_choiceREVERSEcalib_data_outliers_new');  % where is the bic choice data table
    load(filename);  % load choice_data
    chtot = chtot_new;
    chtot.task(strcmp(chtot.task,'TD_ID'),:)={'IvD'};
    chtot.task(strcmp(chtot.task,'TD_DD'),:)={'DvD'};
    
    %keep only first visits
    chtot(chtot.visit >1,:)=[];
    
    %exclusion criteria
    chtot = chtot(strcmp(chtot.exclusion_reason,'RAS'),:);
    chtot_frontomesial = chtot(strcmp(chtot.LocEM,'frontomesial'),:);
    chtot_frontolateral= chtot(strcmp(chtot.LocEM,'frontolateral'),:);
    chtot = [chtot_frontomesial;chtot_frontolateral];
    
    %select only full experiment performers
    list_sub = unique(chtot.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = chtot(strcmp(chtot.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part_n)) == 3
            par = [par ; par_sub];
        else
        end
    end
    chtot = par;
    
    chSS_mean = varfun(@nanmean ,chtot,'InputVariables','choiceSS','GroupingVariables',{'subvisit','LocEM','part','part_n'});
    
    n_part = unique(chSS_mean.part);
    for i_pt = 0 : max(n_part)
        %select the part we want to investigate
        idx_tmp = chSS_mean.part == i_pt;
        ch_tmp = chSS_mean(idx_tmp, :);
        %dispatch according to group
        chSS_grp2 = ch_tmp(strcmp(ch_tmp.LocEM, 'frontomesial'), :);
        chSS_grp1 = ch_tmp(strcmp(ch_tmp.LocEM, 'frontolateral'), :);
        %ttest accordingly
        %[H,P,CI,STATS]
        [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.nanmean_choiceSS, chSS_grp1.nanmean_choiceSS);
        
        % write result to text file
        diary([savewhere '.txt']);
        disp(strcat('part: ', sprintf('%0.f',i_pt)));
        disp(strcat('mean frontolateral (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.nanmean_choiceSS))));
        disp(strcat('mean frontomesial (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.nanmean_choiceSS))));
        
        disp(strcat('frontolateral vs frontomesial - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
        
        diary off
        clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
        
        
    end
     chSS_mean2 = varfun(@nanmean ,chtot,'InputVariables','choiceSS','GroupingVariables',{'subvisit','LocEM'});
   
    %dispatch according to group
    chSS_grp2 = chSS_mean2(strcmp(chSS_mean2.LocEM, 'frontomesial'), :);
    chSS_grp1 = chSS_mean2(strcmp(chSS_mean2.LocEM, 'frontolateral'), :);
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.nanmean_choiceSS, chSS_grp1.nanmean_choiceSS);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('run pooled together'));
    disp(strcat('mean frontolateral (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.nanmean_choiceSS))));
    disp(strcat('mean frontomesial (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.nanmean_choiceSS))));
    
    disp(strcat('frontolateral vs frontomesial - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    
    diary off
    clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'choiceSS ~ 1 + part_n + LocEM + part_n:LocEM ';
    model_formula2 = 'choiceSS ~ 1 + part + LocEM + part:LocEM ';
    % do the analysis
    chtot = sortrows(chtot,'part_n','ascend'); %put calib part first (as default)
    chtot = sortrows(chtot,'type_cond','ascend'); %put control first (as default)
    res1 = fitglm(chtot, model_formula ,'distribution', 'binomial', 'link', 'logit');
    res2 = fitglm(chtot, model_formula2,'distribution', 'binomial', 'link', 'logit');
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM MIXED TASKS')
    disp('BLOCK AS NOMINAL VARIABLE')
    disp(res1);
    disp('BLOCK AS CONTINUOUS VARIABLE')
    disp(res2);
    disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
    diary off;
    
    % define size of figure by number of sub plots
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 3;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    
    % gramm - plot
    g = gramm('x', chSS_mean.part, 'y', chSS_mean.nanmean_choiceSS*100, 'color', chSS_mean.LocEM);  % define data
    g.stat_summary('type','sem','geom',{'line' 'point', 'errorbar'},'dodge',0.4);% this does plot the mean point
    g.set_point_options('base_size',9);
    g.set_line_options('styles', {':'})
    g.axe_property('XLim',[-0.5 2.5] ,'YLim',[30 70] ,'XTick', [0:2],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'impulsive choice rate', 'color', 'group', 'marker', 'task');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end % if files does not exist
