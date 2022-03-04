% This script analyses rating data from MRSFATIGUE_02
%
% Author: Antonius Wiehler <antonius.wiehler@gmail.com>
% Original: 2018-05-15, copied from gameday_pilot_05_analyze_choice_01.m
% Modified: 2018-05-15

%% PREPARATION
% =========================================================================
clear all;
close all;
clc;


%% SET CONFIGURATION PARAMETERS
% =========================================================================
% this should everything that is used in multiple scripts.

cfg.studyname = 'GLIOMA_FATIGUE';

cfg.conditions = {'GLIOMAFATIGUE01control', 'GLIOMAFATIGUE01preOP', 'GLIOMAFATIGUE01postOP'}; % name of the condition to analyze


cfg.mainexp.n_sessions = 1;

% directories
cfg.dir.meta         = 'outputs/meta/';
cfg.dir.rating       = 'outputs/ratings/';
cfg.dir.group        = 'analyze/ratings/group/';
cfg.dir.plots        = 'analyze/ratings/within_subject/';
cfg.dir.group_plots  = 'analyze/ratings/group/plots/';

cfg.dir.modelling       = 'analyze/choice/modelling/'; %choice ref


create_missing_directories(cfg.dir);

cfg.dir.export_fig   = 'export_fig/';
addpath(genpath(cfg.dir.export_fig));  % add export_fig for saving plots

% figure
cfg.fig.width              = 20; % cm
cfg.fig.height             = 10;   % cm
cfg.fig.width_per_subplot  = cfg.fig.width / 4.5; % cm
cfg.fig.height_per_subplot = cfg.fig.height / 2; % cm
cfg.fig.height             = 10;   % cm
cfg.fig.fsiz               = 15;    % font size
cfg.fig.offset             = 0.5;  % offset of x axis to better plot two groups

cfg.fig.colormap     = [ 0    0.4470    0.7410;  % blue
    0.8500    0.3250    0.0980;  % red/orange
    0.9290    0.6940    0.1250;  % yellow
    0.4940    0.1840    0.5560;  % purple
    0.4660    0.6740    0.1880;  % green
    0.3010    0.7450    0.9330;  % light blue
    0.6350    0.0780    0.1840];  % dark red


%% COMBINE THE DATA OF ALL SUBJECTS
% =========================================================================
combine_rating_data_group_v01(cfg);



%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_01_ratings_patient', cfg.dir.group_plots, cfg.studyname);   % construct file name (%s = convert string the cfg)

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5 ;  % 4 plus legend
    n_y_subplots = 2.5 ;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_rating_data.mat');  % whre to save the final table?
    load(filename);  % load rating_data
    rating_data(strcmp(rating_data.part_n,'switch_p'),:)= []; %erase all the data from second visit
    
    
    q_name = unique(rating_data.question);
    
    for i_pn = 1 : length(q_name)
        q_data = rating_data(strcmp(rating_data.question, q_name{i_pn}), :);
        n_part = unique(q_data.part_n);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(q_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            test_data = q_data(strcmp(q_data.part_n,n_part{i_pt}), :);
            
            %dispatch according to group
            chSS_grp2 = test_data(strcmp(test_data.type_cond, 'patient'), :);
            chSS_grp1 = test_data(strcmp(test_data.type_cond, 'control'), :);
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.rating, chSS_grp1.rating);
            
            % write result to text file
            diary([savewhere '.txt']);
            sprintf(n_part{i_pt})
            disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.rating))));
            disp(strcat('mean patient (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.rating))));
            
            disp(strcat('patient vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            
            diary off
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
        end
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'rating ~ 1 + visit + part_n + type_cond + part_n:type_cond ';
        % do the analysis
        q_data = sortrows(q_data,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(q_data, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        diary off;
        
    end
    
    % gramm
    g = gramm('x', rating_data.part_n, 'y', rating_data.rating, 'color', rating_data.type_cond);  % define data
    g.stat_summary('geom',{'point', 'line', 'errorbar'}, 'setylim', true, 'type', 'sem');  % this does plot the mean point
    g.facet_grid('', rating_data.question);  % split in subplots
    g.axe_property('YLim', [0 100], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'time', 'y', 'rating', 'color', 'group');
    g.set_color_options('map',cfg.fig.colormap);     g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_02_ratings_condition', cfg.dir.group_plots, cfg.studyname);   % construct file name (%s = convert string the cfg)

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5 ;  % 4 plus legend
    n_y_subplots = 2.5 ;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_rating_data.mat');  % whre to save the final table?
    load(filename);  % load rating_data
    rating_data(strcmp(rating_data.part_n,'switch_p'),:)= []; %erase all the data from second visit
    
    
    q_name = unique(rating_data.question);
    
    for i_pn = 1 : length(q_name)
        q_data = rating_data(strcmp(rating_data.question, q_name{i_pn}), :);
        n_part = unique(q_data.part_n);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(q_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            test_data = q_data(strcmp(q_data.part_n,n_part{i_pt}), :);
            
            %dispatch according to group
            chSS_grp3 = test_data(strcmp(test_data.new_condition, 'preOP'), :);
            chSS_grp2 = test_data(strcmp(test_data.new_condition, 'postOP'), :);
            chSS_grp1 = test_data(strcmp(test_data.new_condition, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.rating, chSS_grp1.rating);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3.rating, chSS_grp1.rating);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3.rating, chSS_grp2.rating);
            
            % write result to text file
            diary([savewhere '.txt']);
            sprintf(n_part{i_pt})
            disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.rating))));
            disp(strcat('mean preOP   (n=', sprintf('%0.f',length(unique(chSS_grp3.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp3.rating))));
            disp(strcat('mean postOP  (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.rating))));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3
        end
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'rating ~ 1 + visit + part_n + new_condition + part_n:new_condition ';
        % do the analysis
        q_data = sortrows(q_data,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(q_data, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        diary off;
        
    end
    
    % gramm
    g = gramm('x', rating_data.part_n, 'y', rating_data.rating, 'color', rating_data.new_condition);  % define data
    g.stat_summary('geom',{'point', 'line', 'errorbar'}, 'setylim', true, 'type', 'sem');  % this does plot the mean point
    g.facet_grid('', rating_data.question);  % split in subplots
    g.axe_property('YLim', [0 100], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'time', 'y', 'rating', 'color', 'group');
    g.set_color_options('map',cfg.fig.colormap);     g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_03_ratings_patient_samesub', cfg.dir.group_plots, cfg.studyname);   % construct file name (%s = convert string the cfg)

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 2 ;  % 4 plus legend
    n_y_subplots = 2 ;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off','DefaultAxesFontSize',12);
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_rating_data.mat');  % whre to save the final table?
    load(filename);  % load rating_data
    rating_data(strcmp(rating_data.part_n,'switch_p'),:)= []; %erase all the data from second visit
    
    %rely on the same subject as cognitive fatigue
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    
    %keep only first visits
    parameters(parameters.visit >1,:)=[];
    
    list_sub = unique(parameters.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = parameters(strcmp(parameters.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part)) == 3
            par = [par ; par_sub];
        else
        end
    end
    parameters = par;
    
    subjects = unique(parameters.subvisit);
    gr=[];
    for i_s = 1 : length(subjects)
        gr_sub = rating_data(strcmp(rating_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    rating_data = gr;
    
    
    q_name = unique(rating_data.question);
    for i_pn = 1 : length(q_name)
        q_data = rating_data(strcmp(rating_data.question, q_name{i_pn}), :);
        n_part = unique(q_data.part_n);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(q_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            test_data = q_data(strcmp(q_data.part_n,n_part{i_pt}), :);
            
            %dispatch according to group
            chSS_grp2 = test_data(strcmp(test_data.type_cond, 'patient'), :);
            chSS_grp1 = test_data(strcmp(test_data.type_cond, 'control'), :);
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.rating, chSS_grp1.rating);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2.rating, chSS_grp1.rating,'Tail','right');
            [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2.rating,chSS_grp1.rating);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1.rating)),') = ',sprintf('%0.3f',mean(chSS_grp1.rating))));
            disp(strcat('mean patient (n=', sprintf('%0.f',length(chSS_grp2.rating)),') = ',sprintf('%0.3f',mean(chSS_grp2.rating))));
            
            disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
            diary off
            
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
        end
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'rating ~ 1 + part_n + type_cond + part_n:type_cond ';
        % do the analysis
        q_data = sortrows(q_data,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(q_data, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        diary off;
        
    end
    
    % define size of figure by number of sub plots
    n_x_subplots = 2;  % 4 plus legend
    n_y_subplots = 2;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    rating_data = rating_data(strcmp(rating_data.question,'Fatigue'),:);
    % gramm
    g = gramm('x', rating_data.part_n, 'y', rating_data.rating, 'color', rating_data.type_cond);  % define data
    g.stat_summary('type','sem','geom',{'point', 'errorbar'},'dodge',0.3);% this does plot the mean point
    g.set_point_options('base_size',8);
    g.axe_property('YLim', [10 70], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'fatigue rating', 'color', 'group');
    g.set_color_options('map',cfg.fig.colormap);     g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_03_ratings_patient_cleaned_samesub', cfg.dir.group_plots, cfg.studyname);   % construct file name (%s = convert string the cfg)

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 2 ;  % 4 plus legend
    n_y_subplots = 2 ;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off','DefaultAxesFontSize',12);
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_rating_data.mat');  % whre to save the final table?
    load(filename);  % load rating_data
    rating_data(strcmp(rating_data.part_n,'switch_p'),:)= []; %erase all the data from second visit
    
    %rely on the same subject as cognitive fatigue
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    
    %exclusion criteria
    parameters = parameters(strcmp(parameters.exclusion_reason,'RAS'),:);
    
    %keep only first visits
    parameters(parameters.visit >1,:)=[];
    
    list_sub = unique(parameters.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = parameters(strcmp(parameters.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part)) == 3
            par = [par ; par_sub];
        else
        end
    end
    parameters = par;
    
    subjects = unique(parameters.subvisit);
    gr=[];
    for i_s = 1 : length(subjects)
        gr_sub = rating_data(strcmp(rating_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    rating_data = gr;
    
    rating_data.vol_preOP = str2double(rating_data.vol_preOP);
        %reencode diploma
    rating_data.diploma_score(:) = 0;
    rating_data.diploma_score(strcmp(rating_data.last_diploma,'Aucun'),:) = 0;
    rating_data.diploma_score(strcmp(rating_data.last_diploma,'CAP, BEP'),:) = 1;
    rating_data.diploma_score(strcmp(rating_data.last_diploma,'BAC'),:) = 2;
    rating_data.diploma_score(strcmp(rating_data.last_diploma,'BAC + 2'),:) = 3;
    rating_data.diploma_score(strcmp(rating_data.last_diploma,'BAC + 3'),:) = 4;
    rating_data.diploma_score(strcmp(rating_data.last_diploma,'BAC + 4'),:) = 5;
    rating_data.diploma_score(strcmp(rating_data.last_diploma,'BAC + 5'),:) = 6;
    rating_data.diploma_score(strcmp(rating_data.last_diploma,'BAC + 8'),:) = 7;

    
    
    q_name = unique(rating_data.question);
    for i_pn = 1 : length(q_name)
        q_data = rating_data(strcmp(rating_data.question, q_name{i_pn}), :);
        n_part = unique(q_data.part_n);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(q_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            test_data = q_data(strcmp(q_data.part_n,n_part{i_pt}), :);
            
            %dispatch according to group
            chSS_grp2 = test_data(strcmp(test_data.type_cond, 'patient'), :);
            chSS_grp1 = test_data(strcmp(test_data.type_cond, 'control'), :);
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.rating, chSS_grp1.rating);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2.rating, chSS_grp1.rating,'Tail','right');
            [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2.rating,chSS_grp1.rating);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1.rating)),') = ',sprintf('%0.3f',mean(chSS_grp1.rating))));
            disp(strcat('mean patient (n=', sprintf('%0.f',length(chSS_grp2.rating)),') = ',sprintf('%0.3f',mean(chSS_grp2.rating))));
            
            disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
            diary off
            
            
                        % statistics -generalized linear regression model
            % set up formula
        model_formula_lesion = 'delta_rating ~ 1 + LocDG + LocFront + vol_preOP';
        model_formula_psysoc = 'delta_rating ~ 1 + age + sex + diploma_score + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
        model_formula_treat  = 'delta_rating ~ 1 + new_condition + radio + chimio + tt_antiep';
        model_formula_orthoass  = 'delta_rating ~ 1 + flu_anx + flu_P + TMT_A + TMT_B + TMT_flex + span_for + span_back ';
            model_formula_orthoass2  = 'delta_rating ~ 1 + flu_anx + flu_P + TMT_A + TMT_flex + span_for + span_back ';
            model_formula_orthoassqual  = 'delta_rating ~ 1 + new_condition + synd_dysexe ';

            % do the analysis
            par_pat = test_data(strcmp(test_data.type_cond,'patient'),:);
            par_ctrl = test_data(strcmp(test_data.type_cond,'control'),:);
            par_pat = sortrows(par_pat,'part_n','ascend'); %put calib part first (as default)
            par_pat = sortrows(par_pat,'new_condition','ascend'); %put control first (as default)
            res1 = fitglm(par_pat, model_formula_lesion);
            res2 = fitglm(par_pat, model_formula_psysoc);
            res2ctrl = fitglm(par_ctrl, model_formula_psysoc);
            res3 = fitglm(par_pat, model_formula_treat);
            res4 = fitglm(par_pat, model_formula_orthoass);
            res4bis = fitglm(par_pat, model_formula_orthoass2);
            res4ter =fitglm(par_pat, model_formula_orthoassqual);
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
            disp('ORTHO ASSESSMENT FACTOR')
            disp(res4);
            disp(res4bis);
            disp(res4ter);
            disp(strcat('R2 =',sprintf('%0.3f',res4.Rsquared.Ordinary)));
            diary off;
            
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp

            
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
            
        end
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'rating ~ 1 + part_n + type_cond + part_n:type_cond ';
        % do the analysis
        q_data = sortrows(q_data,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(q_data, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        diary off;
        
    end
    
    % define size of figure by number of sub plots
    n_x_subplots = 2;  % 4 plus legend
    n_y_subplots = 2;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    rating_data = rating_data(strcmp(rating_data.question,'Fatigue'),:);
        n_subjects = length(unique(rating_data.subvisit));
    y = rating_data.rating;
    sem_subject = @(y)([nanmean(y); (nanmean(y) - nanstd(y)./ sqrt(n_subjects)); (nanmean(y) + nanstd(y)./ sqrt(n_subjects))]);

    % gramm
    g = gramm('x', rating_data.part_n, 'y', rating_data.rating, 'color', rating_data.type_cond);  % define data
    g.stat_summary('type',sem_subject,'geom',{'line' 'point', 'errorbar'},'dodge',0.3);% this does plot the mean point
    g.set_point_options('base_size',8);
    g.set_line_options('styles', {':'})
    g.axe_property('YLim', [10 70], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'fatigue rating', 'color', 'group');
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_04_ratings_condition_samesub', cfg.dir.group_plots, cfg.studyname);   % construct file name (%s = convert string the cfg)

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5 ;  % 4 plus legend
    n_y_subplots = 2.5 ;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_rating_data.mat');  % whre to save the final table?
    load(filename);  % load rating_data
    rating_data(strcmp(rating_data.part_n,'switch_p'),:)= []; %erase all the data from second visit
    
    %rely on the same subject as cognitive fatigue
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    list_sub = unique(parameters.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = parameters(strcmp(parameters.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part)) == 3
            par = [par ; par_sub];
        else
        end
    end
    parameters = par;
    subjects = unique(parameters.subvisit);
    gr=[];
    for i_s = 1 : length(subjects)
        gr_sub = rating_data(strcmp(rating_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    rating_data = gr;
    
    
    
    q_name = unique(rating_data.question);
    for i_pn = 1 : length(q_name)
        q_data = rating_data(strcmp(rating_data.question, q_name{i_pn}), :);
        n_part = unique(q_data.part_n);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(q_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            test_data = q_data(strcmp(q_data.part_n,n_part{i_pt}), :);
            
            %dispatch according to group
            chSS_grp3 = test_data(strcmp(test_data.new_condition, 'preOP'), :);
            chSS_grp2 = test_data(strcmp(test_data.new_condition, 'postOP'), :);
            chSS_grp1 = test_data(strcmp(test_data.new_condition, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.rating, chSS_grp1.rating);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3.rating, chSS_grp1.rating);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3.rating, chSS_grp2.rating);
            
            % write result to text file
            diary([savewhere '.txt']);
            sprintf(n_part{i_pt})
            disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.rating))));
            disp(strcat('mean preOP   (n=', sprintf('%0.f',length(unique(chSS_grp3.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp3.rating))));
            disp(strcat('mean postOP  (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.rating))));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3
        end
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'rating ~ 1 + visit + part_n + new_condition + part_n:new_condition ';
        % do the analysis
        q_data = sortrows(q_data,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(q_data, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
    end
    
    % gramm
    g = gramm('x', rating_data.part_n, 'y', rating_data.rating, 'color', rating_data.new_condition);  % define data
    g.stat_summary('type','sem','geom',{'point', 'errorbar'},'dodge',0.3);% this does plot the mean point
    g.set_point_options('base_size',8);
    g.facet_grid('', rating_data.question);  % split in subplots
    g.axe_property('YLim', [0 100], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'parts', 'y', 'rating', 'color', 'condition');
    g.set_color_options('map',cfg.fig.colormap);     g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_05_ratings_hemis', cfg.dir.group_plots, cfg.studyname);   % construct file name (%s = convert string the cfg)

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5 ;  % 4 plus legend
    n_y_subplots = 2.5 ;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_rating_data.mat');  % whre to save the final table?
    load(filename);  % load rating_data
    rating_data(strcmp(rating_data.part_n,'switch_p'),:)= []; %erase all the data from second visit
    
    
    q_name = unique(rating_data.question);
    
    for i_pn = 1 : length(q_name)
        q_data = rating_data(strcmp(rating_data.question, q_name{i_pn}), :);
        n_part = unique(q_data.part_n);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(q_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            test_data = q_data(strcmp(q_data.part_n,n_part{i_pt}), :);
            
            %dispatch according to group
            chSS_grp3 = test_data(strcmp(test_data.LocDG, 'right'), :);
            chSS_grp2 = test_data(strcmp(test_data.LocDG, 'left'), :);
            chSS_grp1 = test_data(strcmp(test_data.LocDG, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.rating, chSS_grp1.rating);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3.rating, chSS_grp1.rating);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3.rating, chSS_grp2.rating);
            
            % write result to text file
            diary([savewhere '.txt']);
            sprintf(n_part{i_pt})
            disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.rating))));
            disp(strcat('mean right   (n=', sprintf('%0.f',length(unique(chSS_grp3.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp3.rating))));
            disp(strcat('mean left    (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.rating))));
            
            disp(strcat('left vs control  - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('right vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('right vs left    - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3
        end
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'rating ~ 1 + visit + part_n + LocDG + part_n:LocDG ';
        % do the analysis
        q_data = sortrows(q_data,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(q_data, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        diary off;
        
    end
    
    % gramm
    g = gramm('x', rating_data.part_n, 'y', rating_data.rating, 'color', rating_data.LocDG);  % define data
    g.stat_summary('geom',{'point', 'line', 'errorbar'}, 'setylim', true, 'type', 'sem');  % this does plot the mean point
    g.facet_grid('', rating_data.question);  % split in subplots
    g.axe_property('YLim', [0 100], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'time', 'y', 'rating', 'color', 'group');
    g.set_color_options('map',cfg.fig.colormap);     g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_05_ratings_sex', cfg.dir.group_plots, cfg.studyname);   % construct file name (%s = convert string the cfg)

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5 ;  % 4 plus legend
    n_y_subplots = 2.5 ;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_rating_data.mat');  % whre to save the final table?
    load(filename);  % load rating_data
    rating_data(strcmp(rating_data.part_n,'switch_p'),:)= []; %erase all the data from second visit
    
    
    q_name = unique(rating_data.question);
    
    for i_pn = 1 : length(q_name)
        q_data = rating_data(strcmp(rating_data.question, q_name{i_pn}), :);
        n_part = unique(q_data.part_n);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(q_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            test_data = q_data(strcmp(q_data.part_n,n_part{i_pt}), :);
            
            %dispatch according to group
            chSS_grp2 = test_data(strcmp(test_data.sex, 'M'), :);
            chSS_grp1 = test_data(strcmp(test_data.sex, 'F'), :);
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.rating, chSS_grp1.rating);
            
            % write result to text file
            diary([savewhere '.txt']);
            sprintf(n_part{i_pt})
            disp(strcat('mean female (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.rating))));
            disp(strcat('mean male   (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.rating))));
            
            disp(strcat('male vs female - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            
            diary off
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
        end
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'rating ~ 1 + visit + part_n + sex + part_n:sex ';
        % do the analysis
        q_data = sortrows(q_data,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(q_data, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        diary off;
        
    end
    
    % gramm
    g = gramm('x', rating_data.part_n, 'y', rating_data.rating, 'color', rating_data.sex);  % define data
    g.stat_summary('geom',{'point', 'line', 'errorbar'}, 'setylim', true, 'type', 'sem');  % this does plot the mean point
    g.facet_grid('', rating_data.question);  % split in subplots
    g.axe_property('YLim', [0 100], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'time', 'y', 'rating', 'color', 'group');
    g.set_color_options('map',cfg.fig.colormap);     g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_05_ratings_ampm', cfg.dir.group_plots, cfg.studyname);   % construct file name (%s = convert string the cfg)

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5 ;  % 4 plus legend
    n_y_subplots = 2.5 ;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_rating_data.mat');  % whre to save the final table?
    load(filename);  % load rating_data
    rating_data(strcmp(rating_data.part_n,'switch_p'),:)= []; %erase all the data from second visit
    
    
    q_name = unique(rating_data.question);
    
    for i_pn = 1 : length(q_name)
        q_data = rating_data(strcmp(rating_data.question, q_name{i_pn}), :);
        n_part = unique(q_data.part_n);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(q_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            test_data = q_data(strcmp(q_data.part_n,n_part{i_pt}), :);
            
            %dispatch according to group
            chSS_grp2 = test_data(strcmp(test_data.am_pm, 'pm'), :);
            chSS_grp1 = test_data(strcmp(test_data.am_pm, 'am'), :);
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.rating, chSS_grp1.rating);
            
            % write result to text file
            diary([savewhere '.txt']);
            sprintf(n_part{i_pt})
            disp(strcat('mean am (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.rating))));
            disp(strcat('mean pm (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.rating))));
            
            disp(strcat('pm vs am - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            
            diary off
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
        end
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'rating ~ 1 + visit + part_n + am_pm + part_n:am_pm ';
        % do the analysis
        q_data = sortrows(q_data,'am_pm','ascend'); %put control first (as default)
        res1 = fitglm(q_data, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        diary off;
        
    end
    
    % gramm
    g = gramm('x', rating_data.part_n, 'y', rating_data.rating, 'color', rating_data.am_pm);  % define data
    g.stat_summary('geom',{'point', 'line', 'errorbar'}, 'setylim', true, 'type', 'sem');  % this does plot the mean point
    g.facet_grid('', rating_data.question);  % split in subplots
    g.axe_property('YLim', [0 100], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'time', 'y', 'rating', 'color', 'group');
    g.set_color_options('map',cfg.fig.colormap);     g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_03_ratings_annoy_patient_cleaned_samesub', cfg.dir.group_plots, cfg.studyname);   % construct file name (%s = convert string the cfg)

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5 ;  % 4 plus legend
    n_y_subplots = 2.5 ;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off','DefaultAxesFontSize',12);
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_rating_data.mat');  % whre to save the final table?
    load(filename);  % load rating_data
    rating_data(strcmp(rating_data.part_n,'switch_p'),:)= []; %erase all the data from second visit
    
    %rely on the same subject as cognitive fatigue
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    
    %exclusion criteria
    parameters = parameters(strcmp(parameters.exclusion_reason,'RAS'),:);
    
    %keep only first visits
    parameters(parameters.visit >1,:)=[];
    
    list_sub = unique(parameters.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = parameters(strcmp(parameters.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part)) == 3
            par = [par ; par_sub];
        else
        end
    end
    parameters = par;
    
    subjects = unique(parameters.subvisit);
    gr=[];
    for i_s = 1 : length(subjects)
        gr_sub = rating_data(strcmp(rating_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    rating_data = gr;
    rating_data.annoy(:) = 0;
    rating_data.annoy_calib = str2double(rating_data.annoy_calib);
    rating_data.annoy_crea = str2double(rating_data.annoy_crea);
    rating_data.annoy_nswitch = str2double(rating_data.annoy_nswitch);
    
    rating_data.annoy(strcmp(rating_data.part_n,'calib'),:)=rating_data.annoy_calib(strcmp(rating_data.part_n,'calib'),:);
    rating_data.annoy(strcmp(rating_data.part_n,'crea'),:)=rating_data.annoy_crea(strcmp(rating_data.part_n,'crea'),:);
    rating_data.annoy(strcmp(rating_data.part_n,'nswitch'),:)=rating_data.annoy_nswitch(strcmp(rating_data.part_n,'nswitch'),:);
    
    n_part = unique(rating_data.part_n);
    
    for i_pt = 1 : length(n_part)
        %select the part we want to investigate
        test_data = rating_data(strcmp(rating_data.part_n,n_part{i_pt}), :);
        
        %dispatch according to group
        chSS_grp2 = test_data(strcmp(test_data.type_cond, 'patient'), :);
        chSS_grp1 = test_data(strcmp(test_data.type_cond, 'control'), :);
        %ttest accordingly
        %[H,P,CI,STATS]
        [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.annoy, chSS_grp1.annoy);
        [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2.annoy, chSS_grp1.annoy,'Tail','right');
        [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2.annoy,chSS_grp1.annoy);
        
        % write result to text file
        diary([savewhere '.txt']);
        disp(strcat('part: ', sprintf('%0.f',i_pt)));
        disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp1.annoy))));
        disp(strcat('mean patient (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp2.annoy))));
        
        disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
        disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
        disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
        diary off
        
        clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
    end
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'annoy ~ 1 + part_n + type_cond + part_n:type_cond ';
    % do the analysis
    rating_data = sortrows(rating_data,'type_cond','ascend'); %put control first (as default)
    res1 = fitglm(rating_data, model_formula);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    diary off;
    
    
    % define size of figure by number of sub plots
    n_x_subplots = 2;  % 4 plus legend
    n_y_subplots = 2;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    rating_data = rating_data(strcmp(rating_data.question,'Fatigue'),:);
    % gramm
    g = gramm('x', rating_data.part_n, 'y', rating_data.annoy, 'color', rating_data.type_cond);  % define data
    g.stat_summary('type','sem','geom',{'point', 'errorbar'},'dodge',0.3);% this does plot the mean point
    g.set_point_options('base_size',8);
    g.axe_property('YLim', [0.5 3.5], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'row', '', 'x', 'run', 'y', 'boredom rating', 'color', 'group');
    g.set_color_options('map',cfg.fig.colormap);     g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist
