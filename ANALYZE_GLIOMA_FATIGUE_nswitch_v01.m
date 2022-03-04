% This script analyses choice data from game day pilots 05
%
% Author: Antonius Wiehler <antonius.wiehler@gmail.com>
% Original: 2018-04-17, copied from TIMECONTROL_analyze_nswitch_v02.m


%% PREPARATION
% =========================================================================
clear all;
close all;
clc;


%% SET CONFIGURATION PARAMETERS
% =========================================================================
% this should be everything that is used in multiple scripts.
cfg.studyname = 'GLIOMA_FATIGUE';

cfg.conditions = {'GLIOMAFATIGUE01control', 'GLIOMAFATIGUE01preOP', 'GLIOMAFATIGUE01postOP'}; % name of the condition to analyze
cfg.task      = 'nswitch';  % task name

cfg.mainexp.n_blocks   = 10;  % from fatstim 05 main testing script
cfg.mainexp.n_sessions = 1;

% directories
cfg.dir.meta         = 'outputs/meta/';
cfg.dir.nswitch      = 'outputs/nswitch/';

cfg.dir.group        = 'analyze/nswitch/group/';
cfg.dir.group_plots  = 'analyze/nswitch/group/plots/';

cfg.dir.modelling       = 'analyze/choice/modelling/'; %choice ref


create_missing_directories(cfg.dir);

% add toolboxes
[~, hostname]=system('hostname');
hostname = deblank(hostname);

if strcmp(hostname, 'UMR-PESSI-WP001')
    cfg.dir.export_fig   = 'C:/Users/student/Documents/Matlab_Toolbox/export_fig-master';
    cfg.dir.gramm = 'C:/Users/student/Documents/Matlab_Toolbox/gramm-master';
    addpath(genpath(cfg.dir.export_fig));  % add export_fig for saving plots
    addpath(genpath(cfg.dir.gramm)); % add gramm for easier plots
end

% figure
cfg.fig.width              = 20; % cm
cfg.fig.height             = 10;   % cm
cfg.fig.width_per_subplot  = cfg.fig.width / 4.5; % cm
cfg.fig.height_per_subplot = cfg.fig.height / 2; % cm
cfg.fig.fsiz               = 7;    % font size
cfg.fig.offset             = 0.05;  % offset of x axis to better plot two groups



cfg.fig.colormap     = [ 0    0.4470    0.7410;  % blue
    0.8500    0.3250    0.0980;  % red/orange
    0.9290    0.6940    0.1250;  % yellow
    0.4940    0.1840    0.5560;  % purple
    0.4660    0.6740    0.1880;  % green
    0.3010    0.7450    0.9330;  % light blue
    0.6350    0.0780    0.1840];  % dark red



%% COMBINE THE DATA OF ALL SUBJECTS PER TASK
% =========================================================================
combine_nswitch_data_group_v02(cfg, cfg.conditions);



%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_01_performance_patient', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 1.5;  % 4 plus legend
    n_y_subplots = 1;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
    %dispatch according to group
    chSS_grp2 = test_data(strcmp(test_data.type_cond, 'patient'), :);
    chSS_grp1 = test_data(strcmp(test_data.type_cond, 'control'), :);
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.performance, chSS_grp1.performance);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.performance))));
    disp(strcat('mean patient (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.performance))));
    
    disp(strcat('patient vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    
    diary off
    clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'performance ~ 1 + visit + trial*type_cond';
    % do the analysis
    test_data = sortrows(test_data,'type_cond','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    diary off;
    
    % gramm
    g = gramm('x', test_data.type_cond, 'y', test_data.performance);  % define data
    g.stat_summary('geom',{'point', 'errorbar'}, 'setylim', true, 'type', 'sem');  % this does plot the mean point
    g.set_names('column', '', 'row', '', 'x', '', 'y', 'performance', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);     g.set_text_options('base_size',17);g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist


%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_01_response_time_patient', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 1.5;  % 4 plus legend
    n_y_subplots = 1;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
    %removing first trial as it is not time limited
    test_data(test_data.trial == 1,:) = [];
    
    %dispatch according to group
    chSS_grp2 = test_data(strcmp(test_data.type_cond, 'patient'), :);
    chSS_grp1 = test_data(strcmp(test_data.type_cond, 'control'), :);
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp1.response_time))));
    disp(strcat('mean patient (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp2.response_time))));
    
    disp(strcat('patient vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    
    diary off
    clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'response_time ~ 1 + visit + trial*type_cond';
    % do the analysis
    test_data = sortrows(test_data,'type_cond','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    diary off;
    
    % gramm
    g = gramm('x', test_data.type_cond, 'y', test_data.response_time);  % define data
    g.stat_summary('geom',{'point', 'errorbar'}, 'setylim', true, 'type', 'sem');  % this does plot the mean point
    g.set_names('column', '', 'row', '', 'x', '', 'y', 'response time (s)', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);     g.set_text_options('base_size',17);g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_02_performance_condition', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 1.5;  % 4 plus legend
    n_y_subplots = 1;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
    %dispatch according to group
    chSS_grp3 = test_data(strcmp(test_data.new_condition, 'preOP'), :);
    chSS_grp2 = test_data(strcmp(test_data.new_condition, 'postOP'), :);
    chSS_grp1 = test_data(strcmp(test_data.new_condition, 'control'), :);
    
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.performance, chSS_grp1.performance);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3.performance, chSS_grp1.performance);
    [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3.performance, chSS_grp2.performance);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.performance))));
    disp(strcat('mean preOP   (n=', sprintf('%0.f',length(unique(chSS_grp3.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp3.performance))));
    disp(strcat('mean postOP  (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.performance))));
    
    disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
    
    diary off
    clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3
    
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'performance ~ 1 + visit + trial*new_condition';
    % do the analysis
    test_data = sortrows(test_data,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    diary off;
    
    
    % gramm
    g = gramm('x', test_data.condition, 'y', test_data.performance);  % define data
    g.stat_summary('geom',{'point', 'errorbar'}, 'setylim', true, 'type', 'sem');  % this does plot the mean point
    g.set_names('column', '', 'row', '', 'x', '', 'y', 'performance', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);     g.set_text_options('base_size',17);g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_02_response_time_condition', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 1.5;  % 4 plus legend
    n_y_subplots = 1;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
    %removing first trial as it is not time limited
    test_data(test_data.trial == 1,:) = [];
    
    %dispatch according to group
    chSS_grp3 = test_data(strcmp(test_data.new_condition, 'preOP'), :);
    chSS_grp2 = test_data(strcmp(test_data.new_condition, 'postOP'), :);
    chSS_grp1 = test_data(strcmp(test_data.new_condition, 'control'), :);
    
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3.response_time, chSS_grp1.response_time);
    [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3.response_time, chSS_grp2.response_time);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp1.response_time))));
    disp(strcat('mean preOP   (n=', sprintf('%0.f',length(unique(chSS_grp3.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp3.response_time))));
    disp(strcat('mean postOP  (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp2.response_time))));
    
    disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
    
    diary off
    clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3
    
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'response_time ~ 1 + visit + trial*new_condition';
    % do the analysis
    test_data = sortrows(test_data,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    diary off;
    
    
    % gramm
    g = gramm('x', test_data.condition, 'y', test_data.response_time);  % define data
    g.stat_summary('geom',{'point', 'errorbar'}, 'setylim', true, 'type', 'sem');  % this does plot the mean point
    g.set_names('column', '', 'row', '', 'x', '', 'y', 'response time (s)', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);     g.set_text_options('base_size',17);g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_02_performance_trial_condition', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 1;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
    %dispatch according to group
    chSS_grp3 = test_data(strcmp(test_data.new_condition, 'preOP'), :);
    chSS_grp2 = test_data(strcmp(test_data.new_condition, 'postOP'), :);
    chSS_grp1 = test_data(strcmp(test_data.new_condition, 'control'), :);
    
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.performance, chSS_grp1.performance);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3.performance, chSS_grp1.performance);
    [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3.performance, chSS_grp2.performance);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.performance))));
    disp(strcat('mean preOP   (n=', sprintf('%0.f',length(unique(chSS_grp3.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp3.performance))));
    disp(strcat('mean postOP  (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.performance))));
    
    disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
    
    diary off
    clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3
    
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'performance ~ 1 + visit + trial*new_condition';
    % do the analysis
    test_data = sortrows(test_data,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    diary off;
    
    
    % gramm
    g = gramm('x', test_data.trial, 'y', test_data.performance, 'color',test_data.new_condition);  % define data
    g.stat_summary('bin_in', 10, 'geom',{'point', 'errorbar','line'}, 'setylim', true, 'type', 'sem');  % this does plot the mean point
    g.set_names('column', '', 'row', '', 'x', '', 'y', 'performance', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);     g.set_text_options('base_size',17);g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_02_response_time_trial_condition', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 1;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
    %removing first trial as it is not time limited
    test_data(test_data.trial == 1,:) = [];
    
    
    %dispatch according to group
    chSS_grp3 = test_data(strcmp(test_data.new_condition, 'preOP'), :);
    chSS_grp2 = test_data(strcmp(test_data.new_condition, 'postOP'), :);
    chSS_grp1 = test_data(strcmp(test_data.new_condition, 'control'), :);
    
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3.response_time, chSS_grp1.response_time);
    [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3.response_time, chSS_grp2.response_time);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp1.response_time))));
    disp(strcat('mean preOP   (n=', sprintf('%0.f',length(unique(chSS_grp3.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp3.response_time))));
    disp(strcat('mean postOP  (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp2.response_time))));
    
    disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
    
    diary off
    clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3
    
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'response_time ~ 1 + visit + trial*new_condition';
    % do the analysis
    test_data = sortrows(test_data,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    diary off;
    
    % gramm
    g = gramm('x', test_data.trial, 'y', test_data.response_time, 'color',test_data.new_condition);  % define data
    g.stat_summary('bin_in', 10, 'geom',{'point', 'errorbar','line'}, 'setylim', true, 'type', 'sem');  % this does plot the mean point
    g.set_names('column', '', 'row', '', 'x', '', 'y', 'response time (s)', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);     g.set_text_options('base_size',17);g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_03_performance_condition_samesub', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 1.5;  % 4 plus legend
    n_y_subplots = 1;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
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
        gr_sub = test_data(strcmp(test_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    test_data = gr;
    
    %dispatch according to group
    chSS_grp3 = test_data(strcmp(test_data.new_condition, 'preOP'), :);
    chSS_grp2 = test_data(strcmp(test_data.new_condition, 'postOP'), :);
    chSS_grp1 = test_data(strcmp(test_data.new_condition, 'control'), :);
    
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.performance, chSS_grp1.performance);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3.performance, chSS_grp1.performance);
    [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3.performance, chSS_grp2.performance);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.performance))));
    disp(strcat('mean preOP   (n=', sprintf('%0.f',length(unique(chSS_grp3.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp3.performance))));
    disp(strcat('mean postOP  (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.performance))));
    
    disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
    
    diary off
    clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3
    
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'performance ~ 1 + visit + trial*new_condition';
    % do the analysis
    test_data = sortrows(test_data,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    diary off;
    
    
    % gramm
    g = gramm('x', test_data.condition, 'y', test_data.performance);  % define data
    g.stat_summary('geom',{'point', 'errorbar'}, 'setylim', true, 'type', 'sem');  % this does plot the mean point
    g.set_names('column', '', 'row', '', 'x', '', 'y', 'performance', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);     g.set_text_options('base_size',17);g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_03_response_time_condition_samesub', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 1.5;  % 4 plus legend
    n_y_subplots = 1;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
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
        gr_sub = test_data(strcmp(test_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    test_data = gr;
    
    %removing first trial as it is not time limited
    test_data(test_data.trial == 1,:) = [];
    
    %dispatch according to group
    chSS_grp3 = test_data(strcmp(test_data.new_condition, 'preOP'), :);
    chSS_grp2 = test_data(strcmp(test_data.new_condition, 'postOP'), :);
    chSS_grp1 = test_data(strcmp(test_data.new_condition, 'control'), :);
    
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3.response_time, chSS_grp1.response_time);
    [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3.response_time, chSS_grp2.response_time);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp1.response_time))));
    disp(strcat('mean preOP   (n=', sprintf('%0.f',length(unique(chSS_grp3.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp3.response_time))));
    disp(strcat('mean postOP  (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp2.response_time))));
    
    disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
    
    diary off
    clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3
    
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'response_time ~ 1 + visit + trial*new_condition';
    % do the analysis
    test_data = sortrows(test_data,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    diary off;
    
    
    % gramm
    g = gramm('x', test_data.condition, 'y', test_data.response_time);  % define data
    g.stat_summary('geom',{'point', 'errorbar'}, 'setylim', true, 'type', 'sem');  % this does plot the mean point
    g.set_names('column', '', 'row', '', 'x', '', 'y', 'response time (s)', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);     g.set_text_options('base_size',17);g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_03_performance_trial_condition_samesub', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 1;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
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
        gr_sub = test_data(strcmp(test_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    test_data = gr;
    
    %dispatch according to group
    chSS_grp3 = test_data(strcmp(test_data.new_condition, 'preOP'), :);
    chSS_grp2 = test_data(strcmp(test_data.new_condition, 'postOP'), :);
    chSS_grp1 = test_data(strcmp(test_data.new_condition, 'control'), :);
    
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.performance, chSS_grp1.performance);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3.performance, chSS_grp1.performance);
    [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3.performance, chSS_grp2.performance);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.performance))));
    disp(strcat('mean preOP   (n=', sprintf('%0.f',length(unique(chSS_grp3.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp3.performance))));
    disp(strcat('mean postOP  (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.performance))));
    
    disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
    
    diary off
    clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3
    
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'performance ~ 1 + visit + trial*new_condition';
    % do the analysis
    test_data = sortrows(test_data,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    diary off;
    
    n_subjects = length(unique(test_data.subvisit));
    y = test_data.performance;
    sem_subject = @(y)([nanmean(y); (nanmean(y) - nanstd(y)./ sqrt(n_subjects)); (nanmean(y) + nanstd(y)./ sqrt(n_subjects))]);
    
    % gramm
    g = gramm('x', test_data.trial, 'y', test_data.performance, 'color',test_data.new_condition);  % define data
    g.stat_summary('bin_in', 10, 'geom',{'point', 'area','line'}, 'setylim', true, 'type', sem_subject);  % this does plot the mean point
    g.set_names('column', '', 'row', '', 'x', '', 'y', 'performance', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);     g.set_text_options('base_size',17);g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_03_response_time_trial_condition_samesub', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 4;  % 4 plus legend
    n_y_subplots = 1;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
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
        gr_sub = test_data(strcmp(test_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    test_data = gr;
    
    %removing first trial as it is not time limited
    test_data(test_data.trial == 1,:) = [];
    
    %dispatch according to group
    chSS_grp3 = test_data(strcmp(test_data.new_condition, 'preOP'), :);
    chSS_grp2 = test_data(strcmp(test_data.new_condition, 'postOP'), :);
    chSS_grp1 = test_data(strcmp(test_data.new_condition, 'control'), :);
    
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3.response_time, chSS_grp1.response_time);
    [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3.response_time, chSS_grp2.response_time);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp1.response_time))));
    disp(strcat('mean preOP   (n=', sprintf('%0.f',length(unique(chSS_grp3.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp3.response_time))));
    disp(strcat('mean postOP  (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp2.response_time))));
    
    disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
    
    diary off
    clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3
    
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'response_time ~ 1 + visit + trial*new_condition';
    % do the analysis
    test_data = sortrows(test_data,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    diary off;
    
    n_subjects = length(unique(test_data.subvisit));
    y = test_data.response_time;
    sem_subject = @(y)([nanmean(y); (nanmean(y) - nanstd(y)./ sqrt(n_subjects)); (nanmean(y) + nanstd(y)./ sqrt(n_subjects))]);
    
    % gramm
    g = gramm('x', test_data.trial, 'y', test_data.response_time, 'color',test_data.new_condition);  % define data
    g.stat_summary('bin_in', 10, 'geom',{'point', 'area','line'}, 'setylim', true, 'type', sem_subject);  % this does plot the mean point
    g.set_names('column', '', 'row', '', 'x', '', 'y', 'response time (s)', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);
    g.set_text_options('base_size',17);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end % if files does not exist


%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_04_response_time_trial_patient_samesub', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
    %rely on the same subject as cognitive fatigue
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    %keep only first visits & included patients
    parameters(parameters.visit >1,:)=[];
    parameters = parameters(strcmp(parameters.exclusion_reason,'RAS'),:);
    
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
        gr_sub = test_data(strcmp(test_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    test_data = gr;
    
    %removing first trial as it is not time limited
    test_data(test_data.trial == 1,:) = [];
    
    %dispatch according to group
    chSS_grp2 = test_data(strcmp(test_data.type_cond, 'patient'), :);
    chSS_grp1 = test_data(strcmp(test_data.type_cond, 'control'), :);
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time,'Tail','right');
    [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2.response_time,chSS_grp1.response_time);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp1.response_time))));
    disp(strcat('mean patient (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp2.response_time))));
    
    disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
    diary off
    
    clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'response_time ~ 1 + trial*type_cond';
    model_formula2 = 'response_time ~ 1 + trial*type_cond + trial*age';
    % do the analysis
    test_data = sortrows(test_data,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    res2= fitglm(test_data, model_formula2);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    disp(res2);
    diary off;
    
    n_x_subplots = 3;  % 2 plus legend
    n_y_subplots = 2;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    n_subjects = length(unique(test_data.subvisit));
    y = test_data.response_time;
    sem_subject = @(y)([nanmean(y); (nanmean(y) - nanstd(y)./ sqrt(n_subjects)); (nanmean(y) + nanstd(y)./ sqrt(n_subjects))]);
    
    % gramm
    g = gramm('x', test_data.block, 'y', test_data.response_time, 'color',test_data.type_cond);  % define data
    g.stat_summary('geom',{'point', 'area','line'},'setylim', false, 'type',sem_subject);  % this does plot the mean point
    g.axe_property('YLim',[1 2]);
    g.set_names('column', '', 'row', '', 'x', 'block', 'y', 'response time (s)', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);
    g.set_text_options('base_size',17);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end % if files does not exist


%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_04_response_time_trial_patient_samesub_switchcost', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
    %rely on the same subject as cognitive fatigue
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    %keep only first visits & included patients
    parameters(parameters.visit >1,:)=[];
    parameters = parameters(strcmp(parameters.exclusion_reason,'RAS'),:);
    
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
        gr_sub = test_data(strcmp(test_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    test_data = gr;
    
    %removing first trial as it is not time limited
    test_data(test_data.trial == 1,:) = [];
    
    %dispatch according to group
    chSS_grp2 = test_data(strcmp(test_data.type_cond, 'patient'), :);
    chSS_grp1 = test_data(strcmp(test_data.type_cond, 'control'), :);
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time,'Tail','right');
    [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2.response_time,chSS_grp1.response_time);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp1.response_time))));
    disp(strcat('mean patient (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp2.response_time))));
    
    disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
    diary off
    
    clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'response_time ~ 1 + trial*type_cond + type_cond*rule_switch';
    model_formula2 = 'response_time ~ 1 + trial*type_cond + type_cond*rule_switch + trial*age';
    % do the analysis
    test_data = sortrows(test_data,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    res2= fitglm(test_data, model_formula2);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    disp(res2);
    diary off;
    
    descrdata = varfun(@nanmean ,test_data,'InputVariables','response_time','GroupingVariables',{'subvisit','type_cond','rule_switch','block'});
    descrdata1 = descrdata(descrdata.rule_switch == 0,:);
    descrdata1.response_time_A = descrdata1.nanmean_response_time;
    descrdata1.response_time_B = descrdata.nanmean_response_time(descrdata.rule_switch == 1,:);
    descrdata1.cost_switch = descrdata1.response_time_B - descrdata1.response_time_A ;
    
    descrdata1.nanmean_response_time =[];
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'cost_switch ~ 1 + block*type_cond ';
    % do the analysis
    descrdata1 = sortrows(descrdata1,'type_cond','ascend'); %put control first (as default)
    res2 = fitglm(descrdata1, model_formula);
    descrdata2 = descrdata1;
    descrdata2(descrdata2.block < 5,:) =[];
    res3 = fitglm(descrdata2, model_formula);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM - ALL BLOCKS')
    disp(res2);
    disp('BLOCK 5 to 25 only');
    disp(res3);
    diary off;
    
    n_x_subplots = 3;  % 2 plus legend
    n_y_subplots = 2;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    
    n_subjects = length(unique(descrdata1.subvisit));
    y = descrdata1.cost_switch;
    sem_subject = @(y)([nanmean(y); (nanmean(y) - nanstd(y)./ sqrt(n_subjects)); (nanmean(y) + nanstd(y)./ sqrt(n_subjects))]);
    % gramm
    g = gramm('x', descrdata1.block, 'y', descrdata1.cost_switch, 'color',descrdata1.type_cond);  % define data
    g.stat_summary('geom', {'point', 'area','line'}, 'setylim', false, 'type', sem_subject);  % this does plot the mean point
    g.axe_property('YLim',[0 1]);
    g.set_names('column', '', 'row', '', 'x', 'block', 'y', 'switch cost (s)', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);
    g.set_text_options('base_size',17);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_04_accuracy_trial_patient_samesub', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
    %rely on the same subject as cognitive fatigue
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    %keep only first visits & included patients
    parameters(parameters.visit >1,:)=[];
    parameters = parameters(strcmp(parameters.exclusion_reason,'RAS'),:);
    
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
        gr_sub = test_data(strcmp(test_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    test_data = gr;
    
    %removing first trial as it is not time limited
    test_data(test_data.trial == 1,:) = [];
    
    %dispatch according to group
    chSS_grp2 = test_data(strcmp(test_data.type_cond, 'patient'), :);
    chSS_grp1 = test_data(strcmp(test_data.type_cond, 'control'), :);
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.performance, chSS_grp1.performance);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2.performance, chSS_grp1.performance,'Tail','right');
    [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2.performance,chSS_grp1.performance);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp1.performance))));
    disp(strcat('mean patient (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp2.performance))));
    
    disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
    diary off
    
    clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'performance ~ 1 + trial*type_cond';
    model_formula2 = 'response_time ~ 1 + trial*type_cond + trial*age';
    % do the analysis
    test_data = sortrows(test_data,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    res2= fitglm(test_data, model_formula2);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    disp(res2);
    diary off;
    
    n_x_subplots = 3;  % 2 plus legend
    n_y_subplots = 2;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    
    n_subjects = length(unique(test_data.subvisit));
    y = test_data.performance;
    sem_subject = @(y)([nanmean(y); (nanmean(y) - nanstd(y)./ sqrt(n_subjects)); (nanmean(y) + nanstd(y)./ sqrt(n_subjects))]);
    
    % gramm
    g = gramm('x', test_data.block, 'y', test_data.performance, 'color',test_data.type_cond);  % define data
    g.stat_summary('geom',{'point', 'area','line'} , 'setylim', false, 'type', sem_subject);  % this does plot the mean point
    g.axe_property('YLim',[0.9 1]);
    g.set_names('column', '', 'row', '', 'x', 'block', 'y', 'accuracy', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);     g.set_text_options('base_size',17);g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_05_med_response_time_trial_patient_samesub', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
    %rely on the same subject as cognitive fatigue
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    %keep only first visits & included patients
    parameters(parameters.visit >1,:)=[];
    parameters = parameters(strcmp(parameters.exclusion_reason,'RAS'),:);
    
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
        gr_sub = test_data(strcmp(test_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    test_data = gr;
    
    %removing first trial as it is not time limited
    test_data(test_data.trial == 1,:) = [];
    
    RT_mean = varfun(@nanmean ,test_data,'InputVariables','response_time','GroupingVariables',{'subvisit','type_cond','block'});
    RT_med = varfun(@nanmedian ,test_data,'InputVariables','response_time','GroupingVariables',{'subvisit','type_cond','block'});
    RT_med.response_time = RT_med.nanmedian_response_time;
    
    %dispatch according to group
    chSS_grp2 = RT_med(strcmp(RT_med.type_cond, 'patient'), :);
    chSS_grp1 = RT_med(strcmp(RT_med.type_cond, 'control'), :);
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time,'Tail','right');
    [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2.response_time,chSS_grp1.response_time);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp1.response_time))));
    disp(strcat('mean patient (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp2.response_time))));
    
    disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
    diary off
    
    clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'response_time ~ 1 + trial*type_cond';
    model_formula2 = 'response_time ~ 1 + trial*type_cond + trial*age';
    % do the analysis
    test_data = sortrows(test_data,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(test_data, model_formula);
    res2= fitglm(test_data, model_formula2);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    disp(res2);
    diary off;
    
    n_x_subplots = 3;  % 2 plus legend
    n_y_subplots = 2;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', RT_med.block, 'y', RT_med.response_time, 'color',RT_med.type_cond);  % define data
    g.stat_summary('geom',{'point', 'area','line'},'setylim', false);  % this does plot the mean point
    g.axe_property('YLim',[0.8 1.8]);
    g.set_names('column', '', 'row', '', 'x', 'block', 'y', 'median response time (s)', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);
    g.set_text_options('base_size',17);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_06_var_response_time_trial_patient_samesub', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
    %rely on the same subject as cognitive fatigue
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    %keep only first visits & included patients
    parameters(parameters.visit >1,:)=[];
    parameters = parameters(strcmp(parameters.exclusion_reason,'RAS'),:);
    
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
        gr_sub = test_data(strcmp(test_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    test_data = gr;
    
    %removing first trial as it is not time limited
    test_data(test_data.trial == 1,:) = [];
    
    RT_mean = varfun(@nanmean ,test_data,'InputVariables','response_time','GroupingVariables',{'subvisit','type_cond','block'});
    RT_var = varfun(@nanvar ,test_data,'InputVariables','response_time','GroupingVariables',{'subvisit','type_cond','block'});
    RT_var.response_time = RT_var.nanvar_response_time;
    
    %dispatch according to group
    chSS_grp2 = RT_var(strcmp(RT_var.type_cond, 'patient'), :);
    chSS_grp1 = RT_var(strcmp(RT_var.type_cond, 'control'), :);
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time,'Tail','right');
    [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2.response_time,chSS_grp1.response_time);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean control (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp1.response_time))));
    disp(strcat('mean patient (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp2.response_time))));
    
    disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
    diary off
    
    clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'response_time ~ 1 + block*type_cond';
    model_formula2 = 'response_time ~ 1 + trial*type_cond + trial*age';
    % do the analysis
    test_data = sortrows(test_data,'new_condition','ascend'); %put control first (as default)
    res1 = fitglm(RT_var, model_formula);
    res2= fitglm(test_data, model_formula2);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    disp(res2);
    diary off;
    
    n_x_subplots = 3;  % 2 plus legend
    n_y_subplots = 2;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', RT_var.block, 'y', RT_var.response_time, 'color',RT_var.type_cond);  % define data
    g.stat_summary('geom',{'point', 'area','line'},'setylim', false);  % this does plot the mean point
    g.axe_property('YLim',[0 1.8]);
    g.set_names('column', '', 'row', '', 'x', 'block', 'y', 'response time variance (s²)', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);
    g.set_text_options('base_size',17);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end % if files does not exist

%% PLOT CONDITION EFFECT IN PERFORMANCE
% =========================================================================
savewhere = sprintf('%s%s_07_response_time_trial_patient_samesub', cfg.dir.group_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    filename = strcat(cfg.dir.group, cfg.studyname, '_group_nswitch_data');
    load(filename);  % load rating_data
    
    
    % CONDITION
    test_data = nswitch_data;
    
    %rely on the same subject as cognitive fatigue
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    %keep only first visits & included patients
    parameters(parameters.visit >1,:)=[];
    parameters = parameters(strcmp(parameters.exclusion_reason,'RAS'),:);
    
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
        gr_sub = test_data(strcmp(test_data.subvisit,subjects{i_s}),:);
        if isempty(gr_sub)
            disp(subjects{i_s})
        else
        end
        gr = [gr;gr_sub];
    end
    
    test_data = gr;
    %removing first trial as it is not time limited
    test_data(test_data.trial == 1,:) = [];
    
    %median-splitting  participant based on age
    RT_mean = varfun(@nanmean ,test_data,'InputVariables','response_time','GroupingVariables',{'subvisit','age','type_cond','block'});
    age_RT = varfun(@nanmean ,test_data,'InputVariables','age','GroupingVariables',{'subvisit','age'});
    age_med  = nanmedian(age_RT.nanmean_age);
    RT_mean.age_med(:) = {'young'};
    RT_mean.age_med(RT_mean.age > age_med,:) = {'old'} ;
    RT_mean.response_time = RT_mean.nanmean_response_time;


    
    %dispatch according to group
    chSS_grp2 = RT_mean(strcmp(RT_mean.age_med, 'old'), :);
    chSS_grp1 = RT_mean(strcmp(RT_mean.age_med, 'young'), :);
    %ttest accordingly
    %[H,P,CI,STATS]
    [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time);
    [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2.response_time, chSS_grp1.response_time,'Tail','right');
    [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2.response_time,chSS_grp1.response_time);
    
    % write result to text file
    diary([savewhere '.txt']);
    disp(strcat('mean young (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp1.response_time))));
    disp(strcat('mean old (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',nanmean(chSS_grp2.response_time))));
    
    disp(strcat('young vs old - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
    disp(strcat('young vs old - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
    disp(strcat('young vs old - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
    diary off
    
    clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
    
    % statistics -generalized linear regression model
    % set up formula
    model_formula = 'response_time ~ 1 + block*type_cond';
    model_formula2 = 'response_time ~ 1 + block*type_cond + block*age_med';
    % do the analysis
    RT_mean = sortrows(RT_mean,'type_cond','ascend'); %put control first (as default)
    res1 = fitglm(RT_mean, model_formula);
    res2= fitglm(RT_mean, model_formula2);
    % write result to text file
    diary([savewhere '.txt']);
    disp('FITGLM')
    disp(res1);
    disp(res2);
    diary off;
    
    n_x_subplots = 3;  % 2 plus legend
    n_y_subplots = 2;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
        
    % gramm
    g = gramm('x', RT_mean.block, 'y', RT_mean.response_time, 'color',RT_mean.age_med);  % define data
    g.stat_summary('geom',{'point', 'area','line'},'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.axe_property('YLim',[0.8 2]);
    g.set_names('column', '', 'row', '', 'x', 'block', 'y', 'response time (s)', 'color','group');
    g.set_color_options('map',cfg.fig.colormap);
    g.set_text_options('base_size',17);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end % if files does not exist
