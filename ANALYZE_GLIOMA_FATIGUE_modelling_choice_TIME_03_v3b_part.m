% This script models choice data from GLIOMA_FATIGUE
%
% Author: Antonius Wiehler <antonius.wiehler@gmail.com>
% Original: 2018-09-14
% Modified: 2018-09-17



%% PREPARATION
% =========================================================================
clear all;
close all;
clc;


%% SET CONFIGURATION PARAMETERS
% =========================================================================
% this should be everything that is used in multiple scripts.
cfg.studyname = 'GLIOMA_FATIGUE';

cfg.conditions = {'control', 'preOP', 'postOP'};

task_list(1, 1) = {'TD_ID'};  % time discounting imm-delayed
task_list(1, 2) = {'TD_DD'};  % time discounting delayed-delayed

task_list = reshape(task_list', [], 1);  % reshape into one line of tasks, to undo this use reshape(task_list, 2, [])'

cfg.mainexp.n_blocks   = 23;  % from main testing script
cfg.mainexp.n_sessions = 1;



% directories

cfg.dir.modelling       = 'analyze/choice/modelling/';
cfg.dir.modelling_plots = 'analyze/choice/modelling/plots/part/';
cfg.dir.modelling_subj  = 'analyze/choice/modelling/subject/part/';
cfg.dir.group           = 'analyze/choice/group/';
cfg.dir.group_plots     = 'analyze/choice/group/plots/';
create_missing_directories(cfg.dir);

cfg.dir.choice_models = 'choice_models/';
addpath(genpath(cfg.dir.choice_models));  % add export_fig for saving plots

% add toolboxes
[~, hostname]=system('hostname');
hostname = deblank(hostname);

if strcmp(hostname, 'UMR-PESSI-WP001')
    cfg.dir.export_fig   = 'C:/Users/student/Documents/Matlab_Toolbox/export_fig-master';
    cfg.dir.gramm = 'C:/Users/student/Documents/Matlab_Toolbox/gramm-master';
    addpath(genpath(cfg.dir.export_fig));  % add export_fig for saving plots
    addpath(genpath(cfg.dir.gramm)); % add gramm for easier plots
    
    % setup VBA toolbox
    if ~contains(path, 'VBA-toolbox')  % if VBA is not in path, run the setup
        cd('C:/Users/student/Documents/Matlab_Toolbox/VBA-toolbox');
        VBA_setup;
        cd('C:\Users\student\Documents\GitHub\GLIOMA_FATIGUE');
        
    end
    
    
    
end

% figure
cfg.fig.width              = 20; % cm
cfg.fig.height             = 10;   % cm
cfg.fig.width_per_subplot  = cfg.fig.width / 4.5; % cm
cfg.fig.height_per_subplot = cfg.fig.height / 2; % cm
cfg.fig.height             = 10;   % cm
cfg.fig.fsiz               = 7;    % font size
cfg.fig.offset             = 0.5;  % offset of x axis to better plot two groups

cfg.fig.colormap     = [ 0    0.4470    0.7410;  % blue
    0.8500    0.3250    0.0980;  % red/orange
    0.9290    0.6940    0.1250;  % yellow
    0.4940    0.1840    0.5560;  % purple
    0.4660    0.6740    0.1880;  % green
    0.3010    0.7450    0.9330;  % light blue
    0.6350    0.0780    0.1840];  % dark red


%% MODEL CALIBRATION AND THE 4 PARTS OF THE EXPERIMENT SEPARATELY
% =========================================================================

% load calibaration data
filename = strcat(cfg.dir.group, cfg.studyname, '_group_choice&calib_data.mat');
tmp = load(filename); % this load calib data
choiceALL_data = tmp.chtot;

choiceALL_data = choiceALL_data(strcmp(choiceALL_data.domain, 'time'), :);

choice_data_full = choiceALL_data;  % save as backup for later

conditions = cfg.conditions;


% load model set up
model_space = choice_model_setup_TIME_v03;  % load models
model_space = model_space([1 4 5]); % models preselections

% run VBA for every condition
for i_c = 1 : length(conditions)
    
    % reduce data to current condition :
    choiceALL_data = choice_data_full(strcmp(choice_data_full.new_condition, conditions{i_c}), :);
    %choiceALL_data = choice_data_full ;
    
    % run VBA for every model
    for i_m = 1 : length(model_space)
        
        close all
        clear in dim g_fname group_options VBA_priors options_group output_subject posterior_subject subjects
        
        new_filename = sprintf('%s%s_VBA_TIME_modelling_fatigue_model_allvis_part_%02i_%s.mat', cfg.dir.modelling, conditions{i_c}, i_m, model_space{i_m}.model_name);
        
        
        subjects = unique(choiceALL_data.subvisit);
        n_subjects = length(subjects);
        
        
        
        % CHECK IF FILE EXISTS
        if ~exist(new_filename, 'file')
            fprintf('Running %s...\n', new_filename);  % print current filename to console
            
            for i_s = 1 : n_subjects
                
                % prepare subject data
                clear subject_data
                subject_data = choiceALL_data(strcmp(choiceALL_data.subvisit, subjects{i_s}), :);
                
                
                for i_h = 1 : 3
                    clear subject_data_tmp
                    if i_h == 1
                        idx_tmp = subject_data.part == 0;
                    elseif i_h == 2
                        idx_tmp = subject_data.part == 1;
                    elseif i_h == 3
                        idx_tmp = subject_data.part == 2;
                    end
                    
                    subject_data_tmp = subject_data(idx_tmp, :);
                    
                    if ~isempty(subject_data_tmp)
                        
                        % reformat data. VBA indices will be set in the choice model setup
                        ntrials(i_s, i_h) = length(subject_data_tmp.choiceSS);
                        u(i_s, i_h) = {[subject_data_tmp.costSS'; subject_data_tmp.rewardSS'; subject_data_tmp.costLL'; subject_data_tmp.rewardLL']};
                        y(i_s, i_h) = {subject_data_tmp.choiceSS'};
                        
                        
                        g_fname = model_space{i_m}.g_fname;  % load function handle from model structure
                        
                        % load parts of dim that depend on model
                        dim = model_space{i_m}.dim;
                        
                        % add parts of dim that depent on data
                        dim.n_t = ntrials(i_s, i_h);  % number of trials
                        
                        % load parts of in that depend on model
                        in = model_space{i_m}.in;
                        
                        % add parts of in that depent on data
                        in.ind.t    = [1; 3]; % index of temporal horizon (in u)
                        in.ind.R    = [2; 4]; % index of expected reward (in u)
                        
                        
                        % Build options for model inversion
                        options.dim        = dim;  % copy dim structure (s. above)
                        options.inG        = in;  %  copy indices (s. above), passed to model function
                        options.priors     = model_space{i_m}.VBA_priors;  % copy prior structure (s. above)
                        
                        % copy priors fom calibration
                        %                     sub_indx = strcmp(calib.subjects, subjects{i_s});
                        %                     options.priors.muPhi = calib.posterior_subject{sub_indx}.muPhi;
                        %                     options.priors.SigmaPhi = calib.posterior_subject{sub_indx}.SigmaPhi;
                        
                        
                        options.binomial   = 1;  % use this if response is 0/1 and not continous
                        options.DisplayWin = 0;  % show display during fitting?
                        
                        % Call inversion routine
                        
                        [posterior_subject{i_s, i_h}, output_subject{i_s, i_h}] = VBA_NLStateSpaceModel(y{i_s, i_h}, u{i_s, i_h}, [], g_fname, dim, options); % PROBLEM HERE ->double           %VBA_ReDisplay(posterior, output);  % to get the graphical output again
                        close all
                        % SAVE RESULT
                        save(new_filename, 'posterior_subject', 'output_subject', 'subjects');
                    else
                    end  % if sub data is not empty
                    
                end  % loop halfs
            end  % for loop subjects
            
        end  % if file does not exist
    end  % loop models
end %loop conditions


% CALIBRATION: MODEL SELECTION
% =========================================================================

savewhere = sprintf('%sVBA_TIME_parameters_fatigue_model_selection_allvis_part.mat', cfg.dir.modelling);

if ~exist(savewhere, 'file')
    
    clear L posterior out options
    
    conditions = cfg.conditions;
    
    
    for i_c = 1 : length(conditions)
        
        search_filename = sprintf('%s%s_VBA_TIME_modelling_fatigue_model_allvis_part_*.mat', cfg.dir.modelling, conditions{i_c});
        files = dir(search_filename);
        files = files (1:3);
        
        % load files
        for i_f = 1 : length(files)
            
            fit = load(fullfile(cfg.dir.modelling, files(i_f).name));
            
            for i_s = 1 : length(fit.subjects)
                if ~isempty(fit.output_subject{i_s})
                    L(i_f, i_s) = fit.output_subject{i_s}.F;
                else
                end
            end  % loop subjects
            
        end  % loop files
        
        
        % DO MODEL SELECTION
        options.DisplayWin = 1;  % show display during fitting?
        [posterior{i_c}, output{i_c}] = VBA_groupBMC(L, options);
        winning_model(i_c) = find(output{i_c}.ep == (max(output{i_c}.ep)));
        
        %         [posterior{1}, output{1}] = VBA_groupBMC(L, options);
        %         winning_model(1) = find(output{1}.ep == (max(output{1}.ep)));
        %
    end  % loop conditions
    
    save(savewhere, 'posterior', 'output', 'winning_model', 'L', 'options', 'conditions');
    
end  % if file does not exist



%% COMBINE PARAMETERS OF ALL MODELS INTO ONE TABLE - MAIN EXP
% =========================================================================

% load model set up
model_space = choice_model_setup_TIME_v03;  % load models
model_space = model_space([1 4 5]); % models preselections

big_filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part', cfg.dir.modelling);

if ~exist(big_filename, 'file')
    
    all_parameters = [];
    
    % run data combination for every condition
    for i_c = 1 : length(conditions)
        
        
        % run data combination for every model
        for i_m = 1 : length(model_space)
            
            
            % choice data of main exp to extract condition
            filename = strcat(cfg.dir.group, cfg.studyname, '_group_choice&calib_data.mat');
            load(filename); % this load choice data
            chtot = chtot(strcmp(chtot.domain, 'time'), :);
            chtot = chtot(strcmp(chtot.new_condition, conditions{i_c}), :);
            
            
            model_filename = sprintf('%s%s_VBA_TIME_modelling_fatigue_model_allvis_part_%02i_%s.mat', cfg.dir.modelling,conditions{i_c}, i_m, model_space{i_m}.model_name);
            modelling = load(model_filename); % this load mainexp choice parameter data
            
            parameters = [];  % to add data later
            
            for i_s = 1 : length(modelling.subjects)
                
                sub_indx_choice = strcmp(chtot.subvisit, modelling.subjects(i_s));
                
                
                % meta data
                choice_tmp = chtot(sub_indx_choice, :);
                
                
                for i_h = 1 : 3
                    clear tmp_session
                    if ~isempty(modelling.output_subject{i_s, i_h})
                        
                        for i_p = 1 : length(model_space{i_m}.parameter_names)  % loop parameters
                            
                            tmp_session = table(choice_tmp.condition(1), choice_tmp.domain(1), choice_tmp.effective_dimensions(1), ...
                                choice_tmp.subvisit(1),choice_tmp.subject_id(1), choice_tmp.date(1), choice_tmp.visit(1),choice_tmp.type_cond(1),choice_tmp.new_condition(1),choice_tmp.version_test(1),...
                                choice_tmp.sex(1), choice_tmp.am_pm(1), choice_tmp.last_diploma(1), choice_tmp.CSP(1), choice_tmp.chimio(1), choice_tmp.radio(1),choice_tmp.tt_antiep(1), choice_tmp.LocR(1),choice_tmp.LocDG(1),choice_tmp.LocFront(1),choice_tmp.LocEM(1),choice_tmp.vol_preOP(1),choice_tmp.days_since_OP(1),choice_tmp.age(1),choice_tmp.exclusion_reason(1),...
                                choice_tmp.FSS_score(1),choice_tmp.HAD_Ascore(1),choice_tmp.HAD_Dscore(1),choice_tmp.STARK_score(1),choice_tmp.BIS_score(1), ...
                                choice_tmp.annoy_calib(1),choice_tmp.annoy_crea(1),choice_tmp.annoy_nswitch(1),...
                                choice_tmp.flu_anx(1),choice_tmp.flu_P(1),choice_tmp.TMT_A(1),choice_tmp.TMT_B(1),choice_tmp.TMT_flex(1),choice_tmp.span_for(1),choice_tmp.span_back(1),choice_tmp.synd_dysexe(1), ...
                                {num2str(i_h)}, {model_space{i_m}.model_name}, modelling.output_subject{i_s, i_h}.F, ...
                                modelling.output_subject{i_s, i_h}.fit.bacc, model_space{i_m}.parameter_names(i_p), modelling.posterior_subject{i_s, i_h}.muPhi(i_p));
                            
                            tmp_session.Properties.VariableNames = {'condition', 'domain', 'effective_dimensions', 'subvisit','subject_id', 'date', 'visit', 'type_cond', 'new_condition','version_test', ...
                                'sex', 'am_pm', 'last_diploma', 'CSP', 'chimio', 'radio', 'tt_antiep','LocR', 'LocDG', 'LocFront' , 'LocEM', 'vol_preOP', 'days_since_OP', 'age','exclusion_reason',...
                                'FSS_score', 'HAD_Ascore', 'HAD_Dscore', 'STARK_score' , 'BIS_score' ,...
                                'annoy_calib', 'annoy_crea', 'annoy_nswitch',...
                                'flu_anx','flu_P','TMT_A','TMT_B','TMT_flex','span_for','span_back','synd_dysexe',...
                                'part','model', 'F', 'bacc', 'parameter_name', 'parameter_value'};
                            
                            parameters = [parameters; tmp_session];  % for model table
                            all_parameters = [all_parameters; tmp_session];  % for all models table
                        end  % loop parameters
                    end  % if data is not empty
                end  % loop halfs
            end  % loop subjects
            
            
        end  % loop models
        
    end % loop conditions
    all_parameters.part = str2double(all_parameters.part);
    
    all_parameters.part_n(all_parameters.part == 1) = cellstr('calib');
    all_parameters.part_n(all_parameters.part == 2) = cellstr('crea');
    all_parameters.part_n(all_parameters.part == 3) = cellstr('nswitch');
    
    %calculate delta from calibration for each parameters
    all_parameters.delta_par(:) = nan;
    
    
    list_model = unique(all_parameters.model);
    for i_m = 1 : length(list_model)
        idx_m = strcmp(all_parameters.model, list_model{i_m});
        
        list_name = unique(all_parameters.parameter_name);
        for i_pn = 1 : length(list_name)
            idx_pn = strcmp(all_parameters.parameter_name, list_name{i_pn});
            
            list_sub = unique(all_parameters.subvisit);
            for i_s = 1 : length(list_sub)
                idx_s =   strcmp(all_parameters.subvisit, list_sub{i_s});
                big_idx = idx_m & idx_pn & idx_s;
                calib = all_parameters(big_idx,:).parameter_value(strcmp(all_parameters(big_idx,:).part_n,'calib'));
                if length(calib) == 1
                    all_parameters.delta_par(big_idx,:) = all_parameters.parameter_value(big_idx,:) - calib ;
                else
                end
            end %loop subvisit
        end   %loop par_name
    end %loop model
    
    all_parameters.real_c = strcat(all_parameters.new_condition,'_',all_parameters.condition);
    
    %reencode diploma
    all_parameters.diploma_score(:) = 0;
    all_parameters.diploma_score(strcmp(all_parameters.last_diploma,'Aucun'),:) = 0;
    all_parameters.diploma_score(strcmp(all_parameters.last_diploma,'CAP, BEP'),:) = 1;
    all_parameters.diploma_score(strcmp(all_parameters.last_diploma,'BAC'),:) = 2;
    all_parameters.diploma_score(strcmp(all_parameters.last_diploma,'BAC + 2'),:) = 3;
    all_parameters.diploma_score(strcmp(all_parameters.last_diploma,'BAC + 3'),:) = 4;
    all_parameters.diploma_score(strcmp(all_parameters.last_diploma,'BAC + 4'),:) = 5;
    all_parameters.diploma_score(strcmp(all_parameters.last_diploma,'BAC + 5'),:) = 6;
    all_parameters.diploma_score(strcmp(all_parameters.last_diploma,'BAC + 8'),:) = 7;
    
    %encode old/new protocol version
    allparameters_new = all_parameters(strcmp(all_parameters.version_test,'new'),:);
    allparameters_old = all_parameters(strcmp(all_parameters.version_test,'old'),:);
    
    save([big_filename '_new.mat'],'allparameters_new');  % save as matlab data
    save([big_filename '_old.mat'],'allparameters_old');  % save as matlab data
    save([big_filename '.mat'],'all_parameters');  % save as matlab data
    
    %checking outliers and removing them
    save1 = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers', cfg.dir.modelling);
    %removing outliers
    filename = strcat(cfg.dir.group, cfg.studyname,  '_outliers_list.mat');
    load(filename); % this load choice data
    for i_o = 1 : length(outliers.subvisit)
        all_parameters(strcmp(all_parameters.subvisit,outliers.subvisit{i_o}),:)=[];
    end
    
    
    allparameters_new = all_parameters(strcmp(all_parameters.version_test,'new'),:);
    allparameters_old = all_parameters(strcmp(all_parameters.version_test,'old'),:);
    save([save1 '_new.mat'],'allparameters_new');  % save as matlab data
    save([save1 '_old.mat'],'allparameters_old');  % save as matlab data
    save([save1 '.mat'],'all_parameters');  % save as matlab data
    
end  % big filename


close all


%% EXPONENTIAL BIAS - PATIENT vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_01_exp_bias_allvis_part_patient', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'patient'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean patient (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('patient vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            
            diary off
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
            
        end %loop part
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + visit + part_n + type_cond + part_n:type_cond ';
        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.type_cond);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.set_title('exponential + bias - all visit')
    g.axe_property('XLim', [0.5 3.5], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist


%% EXPONENTIAL BIAS - PATIENT vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_01_exp_bias_allvis_part_patient_onevis', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'patient'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2, chSS_grp1,'Tail','right');
            [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2,chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean patient (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
            diary off
            
            
            save_1 = sprintf(strcat(savewhere, '_deltapar_', par_name{i_pn} , '_',num2str(i_pt)));   % construct file name
            if ~exist([save_1 '.png'], 'file')  % only run, when output does not exist
                
                f = figure(1); clf;
                set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
                    'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
                    'PaperUnits','centimeters', 'Visible', 'off');
                
                g = gramm('x', par_tmp.delta_par, 'color', par_tmp.type_cond);  % define data
                g.stat_bin('nbins',6,'geom','bar'); % this does plot the mean point
                g.set_title(strcat(par_name{i_pn},' - delta - ',num2str(i_pt)))
                if i_pn == 1
                    g.axe_property('XLim', [-10 10]);
                elseif i_pn == 2
                    g.axe_property('XLim', [-4 4]);
                elseif i_pn == 3
                    g.axe_property('XLim', [-6 6]);
                else
                end
                g.set_names('column', '', 'x', 'delta distribution', 'y', '', 'color', 'group');
                g.set_color_options('map',cfg.fig.colormap);
                g.draw();
                
                % save figure
                export_fig([save_1 '.png'], '-r600', '-nocrop');  % save figure
                close(f);  % close figure
                
                save_2 = sprintf(strcat(savewhere, '_parameter_', par_name{i_pn} , '_',num2str(i_pt)));   % construct file name
                
                f = figure(1); clf;
                set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
                    'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
                    'PaperUnits','centimeters', 'Visible', 'off');
                
                g = gramm('x', par_tmp.parameter_value, 'color', par_tmp.type_cond);  % define data
                g.stat_bin('nbins',6,'geom','bar'); % this does plot the mean point
                g.set_title(strcat(par_name{i_pn},' - pardis - ',num2str(i_pt)))
                g.set_names('column', '', 'x', 'parameter distribution', 'y', '', 'color', 'group');
                g.set_text_options('base_size',17);
                g.set_color_options('map',cfg.fig.colormap);
                g.draw();
                
                % save figure
                export_fig([save_2 '.png'], '-r600', '-nocrop');  % save figure
                close(f);  % close figure
                
                
            end
            
            % statistics -generalized linear regression model
            % set up formula
            model_formula = 'delta_par ~ 1 + new_condition + LocDG + LocFront + sex + age + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
            % do the analysis
            par_pat = par_tmp(strcmp(par_tmp.type_cond,'patient'),:);
            par_pat = sortrows(par_pat,'part_n','ascend'); %put calib part first (as default)
            par_pat = sortrows(par_pat,'type_cond','ascend'); %put control first (as default)
            res1 = fitglm(par_pat, model_formula);
            % write result to text file
            diary([savewhere '.txt']);
            disp('FITGLM')
            disp(res1);
            disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
            diary off;
            
            
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
            
        end %loop part
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + part_n + type_cond + part_n:type_cond ';
        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.type_cond);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.set_title('exponential + bias - all visit')
    g.axe_property('XLim', [0.5 3.5], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PATIENT vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_01_exp_bias_allvis_part_patient_onevis_ns', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'patient'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2, chSS_grp1,'Tail','right');
            [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2,chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean patient (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
            diary off
            
            
            % statistics -generalized linear regression model
            % set up formula
            model_formula_lesion = 'parameter_value ~ 1 + LocDG + LocFront + vol_preOP';
            model_formula_psysoc = 'parameter_value ~ 1 + age + sex + diploma_score + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
            model_formula_treat  = 'parameter_value ~ 1 + new_condition + radio + chimio + tt_antiep';
            model_formula_orthoass  = 'parameter_value ~ 1 + flu_anx + flu_P + TMT_A + TMT_B + TMT_flex + span_for + span_back ';
            
            % do the analysis
            par_pat = par_tmp(strcmp(par_tmp.type_cond,'patient'),:);
            par_ctrl = par_tmp(strcmp(par_tmp.type_cond,'control'),:);
            par_pat = sortrows(par_pat,'part_n','ascend'); %put calib part first (as default)
            par_pat = sortrows(par_pat,'new_condition','ascend'); %put control first (as default)
            res1 = fitglm(par_pat, model_formula_lesion);
            res2 = fitglm(par_pat, model_formula_psysoc);
            res2ctrl = fitglm(par_ctrl, model_formula_psysoc);
            res3 = fitglm(par_pat, model_formula_treat);
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
            disp('ORTHO ASSESSMENT FACTOR')
            disp(res4);
            disp(strcat('R2 =',sprintf('%0.3f',res4.Rsquared.Ordinary)));
            diary off;
            
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
            
        end %loop part
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + part_n*type_cond ';
        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
    end %loop parameter name
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    parameters = parameters(strcmp(parameters.part_n,'nswitch'),:);
    % gramm
    g = gramm('x', parameters.part_n, 'y', parameters.parameter_value, 'color', parameters.type_cond);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [4 7];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-7 -5];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PATIENT vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_01_exp_bias_allvis_cleaned_patient_onevis', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    
    par_corr = parameters(strcmp(parameters.parameter_name,'bias'),:);
    subv = unique(par_corr.subvisit);
    corr = [];
    for i_s = 1:length(subv)
        idx_calib = strcmp(par_corr.subvisit,subv(i_s)) & strcmp(par_corr.part_n,'calib');
        idx_delta = strcmp(par_corr.subvisit,subv(i_s)) & strcmp(par_corr.part_n,'nswitch');
        calib = par_corr.parameter_value(idx_calib,:);
        delta = par_corr.delta_par(idx_delta,:);
        type_cond = par_corr.type_cond(idx_delta,:);
        corr_sub = [calib,delta,type_cond];
        corr = [corr;corr_sub];
    end
    
    [R,P,RL,RU]= corrcoef(cell2mat(corr(:,1)),cell2mat(corr(:,2)));
    
    % gramm
    g = gramm('x', corr(:,1), 'y', corr(:,2), 'color', corr(:,3));  % define data
    g.geom_point();
    g.set_names('column', '', 'x', 'CALIB', 'y', 'delta switch', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + part_n + type_cond + part_n:type_cond ';
        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'type_cond','ascend'); %put control first (as default)
        resa = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM parameter value')
        disp(resa);
        disp(strcat('R2 =',sprintf('%0.3f',resa.Rsquared.Ordinary)));
        diary off;
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'delta_par ~ 1 + part_n + type_cond + part_n:type_cond ';
        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM delta parameter')
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.delta_par(strcmp(par_tmp.type_cond, 'patient'), :);
            chSS_grp1 = par_tmp.delta_par(strcmp(par_tmp.type_cond, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2, chSS_grp1,'Tail','right');
            [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2,chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean patient (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
            diary off
            
            
            save_1 = sprintf(strcat(savewhere, '_deltapar_', par_name{i_pn} , '_',num2str(i_pt)));   % construct file name
            if ~exist([save_1 '.png'], 'file')  % only run, when output does not exist
                
                f = figure(1); clf;
                set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
                    'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
                    'PaperUnits','centimeters', 'Visible', 'off');
                
                g = gramm('x', par_tmp.delta_par, 'color', par_tmp.type_cond);  % define data
                g.stat_bin('nbins',6,'geom','bar'); % this does plot the mean point
                g.set_title(strcat(par_name{i_pn},' - delta - ',num2str(i_pt)))
                if i_pn == 1
                    g.axe_property('XLim', [-10 10]);
                elseif i_pn == 2
                    g.axe_property('XLim', [-4 4]);
                elseif i_pn == 3
                    g.axe_property('XLim', [-6 6]);
                else
                end
                g.set_names('column', '', 'x', 'delta distribution', 'y', '', 'color', 'group');
                g.set_text_options('base_size',17);
                g.set_color_options('map',cfg.fig.colormap);
                g.draw();
                
                % save figure
                export_fig([save_1 '.png'], '-r600', '-nocrop');  % save figure
                close(f);  % close figure
                
                save_2 = sprintf(strcat(savewhere, '_parameter_', par_name{i_pn} , '_',num2str(i_pt)));   % construct file name
                
                f = figure(1); clf;
                set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
                    'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
                    'PaperUnits','centimeters', 'Visible', 'off');
                
                g = gramm('x', par_tmp.parameter_value, 'color', par_tmp.type_cond);  % define data
                g.stat_bin('nbins',6,'geom','bar'); % this does plot the mean point
                g.set_title(strcat(par_name{i_pn},' - pardis - ',num2str(i_pt)))
                g.set_names('column', '', 'x', 'parameter distribution', 'y', '', 'color', 'group');
                g.set_text_options('base_size',17);
                g.set_color_options('map',cfg.fig.colormap);
                g.draw();
                
                % save figure
                export_fig([save_2 '.png'], '-r600', '-nocrop');  % save figure
                close(f);  % close figure
                
                
            end
            
            
            % statistics -generalized linear regression model
            % set up formula
            model_formula_lesion = 'delta_par ~ 1 + LocDG + LocFront + vol_preOP';
            model_formula_psysoc = 'delta_par ~ 1 + age + sex + diploma_score + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
            model_formula_treat  = 'delta_par ~ 1 + new_condition + radio + chimio + tt_antiep';
            model_formula_orthoass  = 'delta_par ~ 1 + flu_anx + flu_P + TMT_A + TMT_B + TMT_flex + span_for + span_back ';
            
            % do the analysis
            par_pat = par_tmp(strcmp(par_tmp.type_cond,'patient'),:);
            par_ctrl = par_tmp(strcmp(par_tmp.type_cond,'control'),:);
            par_pat = sortrows(par_pat,'part_n','ascend'); %put calib part first (as default)
            par_pat = sortrows(par_pat,'new_condition','ascend'); %put control first (as default)
            res1 = fitglm(par_pat, model_formula_lesion);
            res2 = fitglm(par_pat, model_formula_psysoc);
            res2ctrl = fitglm(par_ctrl, model_formula_psysoc);
            res3 = fitglm(par_pat, model_formula_treat);
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
            disp('ORTHO ASSESSMENT FACTOR')
            disp(res4);
            disp(strcat('R2 =',sprintf('%0.3f',res4.Rsquared.Ordinary)));
            diary off;
            
            
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
            
        end %loop part
        
        
        clear par
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.type_cond);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.set_title('exponential + bias - all visit')
    g.axe_property('XLim', [0.5 3.5], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PATIENT vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_01_exp_bias_allvis_cleaned_patient_onevis_ns', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'patient'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2, chSS_grp1,'Tail','right');
            [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2,chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean patient (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
            diary off
            
            % statistics -generalized linear regression model
            % set up formula
            model_formula_lesion = 'parameter_value ~ 1 + LocDG + LocFront + vol_preOP';
            model_formula_psysoc = 'parameter_value ~ 1 + age + sex + diploma_score + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
            model_formula_treat  = 'parameter_value ~ 1 + new_condition + radio + chimio + tt_antiep';
            model_formula_orthoass  = 'parameter_value ~ 1 + flu_anx + flu_P + TMT_A + TMT_B + TMT_flex + span_for + span_back ';
            model_formula_orthoass2  = 'parameter_value ~ 1 + flu_anx + flu_P + TMT_A + TMT_flex + span_for + span_back ';
            model_formula_orthoassqual  = 'parameter_value ~ 1 + new_condition + synd_dysexe ';
            
            % do the analysis
            par_pat = par_tmp(strcmp(par_tmp.type_cond,'patient'),:);
            par_ctrl = par_tmp(strcmp(par_tmp.type_cond,'control'),:);
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
            
        end %loop part
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + part_n*type_cond ';
        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
    end %loop parameter name
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 1.5;
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    parameters = parameters(strcmp(parameters.part_n,'nswitch'),:);
    % gramm
    g = gramm('x', parameters.part_n, 'y', parameters.parameter_value, 'color', parameters.type_cond);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_point_options('base_size',9);
    g.axe_property('XTickLabel', {'Switch'});
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.set_names('column', '', 'x', '', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [3 7];
    g.facet_axes_handles(1,2).YLim = [-0.5 1.5];
    g.facet_axes_handles(1,3).YLim = [-8 -4];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PATIENT vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_01_exp_bias_allvis_cleaned_hemis_onevis_ns', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'patient'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2, chSS_grp1,'Tail','right');
            [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2,chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean patient (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
            diary off
            
            % statistics -generalized linear regression model
            % set up formula
            model_formula_lesion = 'parameter_value ~ 1 + LocDG + LocFront + vol_preOP';
            model_formula_psysoc = 'parameter_value ~ 1 + age + sex + diploma_score + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
            model_formula_treat  = 'parameter_value ~ 1 + new_condition + radio + chimio + tt_antiep';
            model_formula_orthoass  = 'parameter_value ~ 1 + flu_anx + flu_P + TMT_A + TMT_B + TMT_flex + span_for + span_back ';
            
            % do the analysis
            par_pat = par_tmp(strcmp(par_tmp.type_cond,'patient'),:);
            par_ctrl = par_tmp(strcmp(par_tmp.type_cond,'control'),:);
            par_pat = sortrows(par_pat,'part_n','ascend'); %put calib part first (as default)
            par_pat = sortrows(par_pat,'new_condition','ascend'); %put control first (as default)
            res1 = fitglm(par_pat, model_formula_lesion);
            res2 = fitglm(par_pat, model_formula_psysoc);
            res2ctrl = fitglm(par_ctrl, model_formula_psysoc);
            res3 = fitglm(par_pat, model_formula_treat);
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
            disp('ORTHO ASSESSMENT FACTOR')
            disp(res4);
            disp(strcat('R2 =',sprintf('%0.3f',res4.Rsquared.Ordinary)));
            diary off;
            
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
            
        end %loop part
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + part_n*type_cond ';
        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    parameters = parameters(strcmp(parameters.part_n,'nswitch'),:);
    % gramm
    g = gramm('x', parameters.part_n, 'y', parameters.parameter_value, 'color', parameters.LocDG);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_point_options('base_size',9);
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [4 7];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-7 -5];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist



%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_02_exp_bias_allvis_part_fullexp_condition', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp3 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'preOP'), :);
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'postOP'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3, chSS_grp1);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3, chSS_grp2);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean preOP   (n=', sprintf('%0.f',length(chSS_grp3)),') = ',sprintf('%0.3f',mean(chSS_grp3))));
            disp(strcat('mean postOP  (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + visit + part_n + new_condition + part_n:new_condition';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula_2 = 'parameter_value ~ 1 + age + visit + part_n + new_condition + part_n:new_condition';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res2 = fitglm(par, model_formula_2);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res2);
        disp(strcat('R2 =',sprintf('%0.3f',res2.Rsquared.Ordinary)));
        diary off;
        
        
        %one by one analysis
        %indexing every condition of analyses
        %calib vs creativity
        idx_c_cr = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'crea'});
        %calib vs nswitch
        idx_c_ns = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'nswitch'});
        
        part_comp = [{idx_c_cr} ; {idx_c_ns}]; %stock all index
        
        %control vs postop
        idx_c_po = strcmp(par.new_condition, {'control'}) | strcmp(par.new_condition, {'postOP'});
        %control vs preop
        idx_c_pr = strcmp(par.new_condition, {'control'}) | strcmp(par.new_condition, {'preOP'});
        %postop vs preop
        idx_po_pr = strcmp(par.new_condition, {'postOP'}) | strcmp(par.new_condition, {'preOP'});
        
        cond_comp = [{idx_c_po};{idx_c_pr};{idx_po_pr}]; %stock all index
        
        %multiple comparison
        for i_loop_cond = 1: length(cond_comp)
            if i_loop_cond == 1
                diary([savewhere '.txt']);
                disp('---CONTROL vs POSTOP---')
                diary off;
            elseif i_loop_cond == 2
                diary([savewhere '.txt']);
                disp('---CONTROL vs PREOP---')
                diary off;
            elseif i_loop_cond == 3
                diary([savewhere '.txt']);
                disp('---PREOP vs POSTOP---')
                diary off;
            else
            end
            
            for i_loop_part = 1 : length(part_comp)
                par_tmp = par(cond_comp{i_loop_cond} & part_comp{i_loop_part},:); %indexing condition comparison
                if i_loop_part == 1
                    diary([savewhere '.txt']);
                    disp('---CALIB vs CREA---')
                    diary off;
                elseif i_loop_part == 2
                    diary([savewhere '.txt']);
                    disp('---CALIB vs NSWITCH---')
                    diary off;
                else
                end
                % statistics -generalized linear regression model
                % set up formula
                model_formula = 'parameter_value ~ 1 + visit + part_n + new_condition + part_n:new_condition';
                % do the analysis
                par_tmp = sortrows(par_tmp,'part_n','ascend'); %put calib part first (as default)
                par_tmp = sortrows(par_tmp,'new_condition','ascend'); %put control first (as default)
                res1 = fitglm(par_tmp, model_formula);
                % write result to text file
                diary([savewhere '.txt']);
                disp(res1);
                disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
                diary off;
            end
        end
        clear par
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.new_condition);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 ,'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_02_exp_bias_allvis_part_fullexp_condition_w_out', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_new.mat', cfg.dir.modelling);
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp3 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'preOP'), :);
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'postOP'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3, chSS_grp1);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3, chSS_grp2);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean preOP   (n=', sprintf('%0.f',length(chSS_grp3)),') = ',sprintf('%0.3f',mean(chSS_grp3))));
            disp(strcat('mean postOP  (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + visit + part_n + new_condition + part_n:new_condition';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        %one by one analysis
        %indexing every condition of analyses
        %calib vs creativity
        idx_c_cr = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'crea'});
        %calib vs nswitch
        idx_c_ns = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'nswitch'});
        
        part_comp = [{idx_c_cr} ; {idx_c_ns}]; %stock all index
        
        %control vs postop
        idx_c_po = strcmp(par.new_condition, {'control'}) | strcmp(par.new_condition, {'postOP'});
        %control vs preop
        idx_c_pr = strcmp(par.new_condition, {'control'}) | strcmp(par.new_condition, {'preOP'});
        %postop vs preop
        idx_po_pr = strcmp(par.new_condition, {'postOP'}) | strcmp(par.new_condition, {'preOP'});
        
        cond_comp = [{idx_c_po};{idx_c_pr};{idx_po_pr}]; %stock all index
        
        %multiple comparison
        for i_loop_cond = 1: length(cond_comp)
            if i_loop_cond == 1
                diary([savewhere '.txt']);
                disp('---CONTROL vs POSTOP---')
                diary off;
            elseif i_loop_cond == 2
                diary([savewhere '.txt']);
                disp('---CONTROL vs PREOP---')
                diary off;
            elseif i_loop_cond == 3
                diary([savewhere '.txt']);
                disp('---PREOP vs POSTOP---')
                diary off;
            else
            end
            
            for i_loop_part = 1 : length(part_comp)
                par_tmp = par(cond_comp{i_loop_cond} & part_comp{i_loop_part},:); %indexing condition comparison
                if i_loop_part == 1
                    diary([savewhere '.txt']);
                    disp('---CALIB vs CREA---')
                    diary off;
                elseif i_loop_part == 2
                    diary([savewhere '.txt']);
                    disp('---CALIB vs NSWITCH---')
                    diary off;
                else
                end
                % statistics -generalized linear regression model
                % set up formula
                model_formula = 'parameter_value ~ 1 + visit + part_n + new_condition + part_n:new_condition';
                % do the analysis
                par_tmp = sortrows(par_tmp,'part_n','ascend'); %put calib part first (as default)
                par_tmp = sortrows(par_tmp,'new_condition','ascend'); %put control first (as default)
                res1 = fitglm(par_tmp, model_formula);
                % write result to text file
                diary([savewhere '.txt']);
                disp(res1);
                disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
                diary off;
            end
        end
        clear par
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.new_condition);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 ,'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist


%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_02_exp_bias_allvis_part_fullexp_condition_onevis', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp3 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'preOP'), :);
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'postOP'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3, chSS_grp1);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3, chSS_grp2);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean preOP   (n=', sprintf('%0.f',length(chSS_grp3)),') = ',sprintf('%0.3f',mean(chSS_grp3))));
            disp(strcat('mean postOP  (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + part_n + new_condition + part_n:new_condition';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        %one by one analysis
        %indexing every condition of analyses
        %calib vs creativity
        idx_c_cr = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'crea'});
        %calib vs nswitch
        idx_c_ns = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'nswitch'});
        
        part_comp = [{idx_c_cr} ; {idx_c_ns}]; %stock all index
        
        %control vs postop
        idx_c_po = strcmp(par.new_condition, {'control'}) | strcmp(par.new_condition, {'postOP'});
        %control vs preop
        idx_c_pr = strcmp(par.new_condition, {'control'}) | strcmp(par.new_condition, {'preOP'});
        %postop vs preop
        idx_po_pr = strcmp(par.new_condition, {'postOP'}) | strcmp(par.new_condition, {'preOP'});
        
        cond_comp = [{idx_c_po};{idx_c_pr};{idx_po_pr}]; %stock all index
        
        %multiple comparison
        for i_loop_cond = 1: length(cond_comp)
            if i_loop_cond == 1
                diary([savewhere '.txt']);
                disp('---CONTROL vs POSTOP---')
                diary off;
            elseif i_loop_cond == 2
                diary([savewhere '.txt']);
                disp('---CONTROL vs PREOP---')
                diary off;
            elseif i_loop_cond == 3
                diary([savewhere '.txt']);
                disp('---PREOP vs POSTOP---')
                diary off;
            else
            end
            
            for i_loop_part = 1 : length(part_comp)
                par_tmp = par(cond_comp{i_loop_cond} & part_comp{i_loop_part},:); %indexing condition comparison
                if i_loop_part == 1
                    diary([savewhere '.txt']);
                    disp('---CALIB vs CREA---')
                    diary off;
                elseif i_loop_part == 2
                    diary([savewhere '.txt']);
                    disp('---CALIB vs NSWITCH---')
                    diary off;
                else
                end
                % statistics -generalized linear regression model
                % set up formula
                model_formula = 'parameter_value ~ 1 + part_n + new_condition + part_n:new_condition';
                % do the analysis
                par_tmp = sortrows(par_tmp,'part_n','ascend'); %put calib part first (as default)
                par_tmp = sortrows(par_tmp,'new_condition','ascend'); %put control first (as default)
                res1 = fitglm(par_tmp, model_formula);
                % write result to text file
                diary([savewhere '.txt']);
                disp(res1);
                disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
                diary off;
            end
        end
        clear par
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.new_condition);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 ,'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_02_exp_bias_allvis_part_fullexp_condition_wo_preop', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    
    %keep only first visits
    parameters(strcmp(parameters.new_condition,'preOP'),:)=[];
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'postOP'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean postOP  (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            
            diary off
            
            % statistics -generalized linear regression model
            % set up formula
            model_formula = 'delta_par ~ 1 + new_condition + LocDG + LocFront + sex + age + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
            % do the analysis
            par_tmp = sortrows(par_tmp,'part_n','ascend'); %put calib part first (as default)
            par_tmp = sortrows(par_tmp,'type_cond','ascend'); %put control first (as default)
            res1 = fitglm(par_tmp, model_formula);
            % write result to text file
            diary([savewhere '.txt']);
            disp('FITGLM')
            disp(res1);
            disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
            diary off;
            
            clear  chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + visit + part_n*new_condition';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        %one by one analysis
        %indexing every condition of analyses
        %calib vs creativity
        idx_c_cr = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'crea'});
        %calib vs nswitch
        idx_c_ns = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'nswitch'});
        
        part_comp = [{idx_c_cr} ; {idx_c_ns}]; %stock all index
        
        %control vs postop
        idx_c_po = strcmp(par.new_condition, {'control'}) | strcmp(par.new_condition, {'postOP'});
        
        cond_comp = [{idx_c_po}]; %stock all index
        
        %multiple comparison
        for i_loop_cond = 1: length(cond_comp)
            if i_loop_cond == 1
                diary([savewhere '.txt']);
                disp('---CONTROL vs POSTOP---')
                diary off;
            else
            end
            
            for i_loop_part = 1 : length(part_comp)
                par_tmp = par(cond_comp{i_loop_cond} & part_comp{i_loop_part},:); %indexing condition comparison
                if i_loop_part == 1
                    diary([savewhere '.txt']);
                    disp('---CALIB vs CREA---')
                    diary off;
                elseif i_loop_part == 2
                    diary([savewhere '.txt']);
                    disp('---CALIB vs NSWITCH---')
                    diary off;
                else
                end
                % statistics -generalized linear regression model
                % set up formula
                model_formula = 'parameter_value ~ 1 + visit + part_n*new_condition';
                % do the analysis
                par_tmp = sortrows(par_tmp,'part_n','ascend'); %put calib part first (as default)
                par_tmp = sortrows(par_tmp,'new_condition','ascend'); %put control first (as default)
                res1 = fitglm(par_tmp, model_formula);
                % write result to text file
                diary([savewhere '.txt']);
                disp(res1);
                disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
                diary off;
            end
        end
        clear par
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.new_condition);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 ,'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_03_exp_bias_allvis_part_shortexp_condition', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    
    idx_vis = parameters.part > 2;
    parameters(idx_vis,:) = []; %keep only the data of the calib & crea
    
    list_sub = unique(parameters.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = parameters(strcmp(parameters.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part)) == 2
            par = [par ; par_sub];
        else
        end
    end
    parameters = par;
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp3 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'preOP'), :);
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'postOP'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3, chSS_grp1);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3, chSS_grp2);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean preOP   (n=', sprintf('%0.f',length(chSS_grp3)),') = ',sprintf('%0.3f',mean(chSS_grp3))));
            disp(strcat('mean postOP  (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + visit + part_n + new_condition + part_n:new_condition';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
    end %loop parameter name
    
    %open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.new_condition);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - short experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [0.5 3.5], 'XTick', [1:3],'XTickLabel', {'CALIB', 'CREA', ''});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_03_exp_bias_allvis_part_shortexp_condition_w_out', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    
    idx_vis = parameters.part > 2;
    parameters(idx_vis,:) = []; %keep only the data of the calib & crea part
    
    list_sub = unique(parameters.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = parameters(strcmp(parameters.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part)) == 2
            par = [par ; par_sub];
        else
        end
    end
    parameters = par;
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp3 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'preOP'), :);
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'postOP'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3, chSS_grp1);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3, chSS_grp2);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean preOP   (n=', sprintf('%0.f',length(chSS_grp3)),') = ',sprintf('%0.3f',mean(chSS_grp3))));
            disp(strcat('mean postOP  (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + visit + part_n + new_condition + part_n:new_condition';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
    end %loop parameter name
    
    %open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.new_condition);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - short experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [0.5 3.5], 'XTick', [1:3],'XTickLabel', {'CALIB', 'CREA', ''});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_04_onlyexp_allvis_part_fullexp_condition', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential'), :);
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp3 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'preOP'), :);
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'postOP'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3, chSS_grp1);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3, chSS_grp2);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean preOP   (n=', sprintf('%0.f',length(chSS_grp3)),') = ',sprintf('%0.3f',mean(chSS_grp3))));
            disp(strcat('mean postOP  (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + visit + part_n + new_condition + part_n:new_condition';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.new_condition);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 ,'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist


%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL - BACC
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_05_exp_bias_allvis_part_fullexp_balanceaccur', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 3;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    
    sum_exclude = varfun(@mean ,parameters,'InputVariables',{'bacc'},'GroupingVariables',{'new_condition','subvisit'});
    sum_exclude = sortrows(sum_exclude,'mean_bacc','ascend');
    sum_exclude = sortrows(sum_exclude,'new_condition','ascend');
    part_exclude = varfun(@mean ,parameters,'InputVariables',{'bacc'},'GroupingVariables',{'part','subvisit'});
    part_exclude = sortrows(part_exclude,'mean_bacc','ascend');
    part_exclude = sortrows(part_exclude,'part','ascend');
    % write result to text file
    diary([savewhere '.txt']);
    disp(sum_exclude);
    disp(part_exclude);
    diary off;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.bacc, 'color', parameters.new_condition);  % define data
    g.stat_summary('geom', {'point', 'errorbar', 'line'}, 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias');
    g.axe_property('YLim',[0 1] ,'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_text_options('base_size',17);
    g.set_names('column', '', 'x', 'block', 'y', 'balance accuracy', 'color', 'group');
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL - BACC
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_05_exp_bias_allvis_part_shortexp_balanceaccur', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    
    idx_vis = parameters.part > 2;
    parameters(idx_vis,:) = []; %keep only the data of the first visit
    
    list_sub = unique(parameters.subvisit);
    par = [];
    for i_l = 1 : length(list_sub)
        par_sub = parameters(strcmp(parameters.subvisit, list_sub{i_l}),:);
        if length(unique(par_sub.part)) == 2
            par = [par ; par_sub];
        else
        end
    end
    parameters = par;
    sum_exclude = varfun(@mean ,parameters,'InputVariables',{'bacc'},'GroupingVariables',{'new_condition','subvisit'});
    sum_exclude = sortrows(sum_exclude,'mean_bacc','ascend');
    sum_exclude = sortrows(sum_exclude,'new_condition','ascend');
    part_exclude = varfun(@mean ,parameters,'InputVariables',{'bacc'},'GroupingVariables',{'part','subvisit'});
    part_exclude = sortrows(part_exclude,'mean_bacc','ascend');
    part_exclude = sortrows(part_exclude,'part','ascend');
    % write result to text file
    diary([savewhere '.txt']);
    disp(sum_exclude);
    disp(part_exclude);
    diary off;
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.bacc, 'color', parameters.new_condition);  % define data
    g.stat_summary('geom', {'point', 'errorbar', 'line'}, 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit');
    g.axe_property('YLim',[0.5 1] ,'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_text_options('base_size',17);
    g.set_names('column', '', 'x', 'block', 'y', 'balance accuracy', 'color', 'group');
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end  % if file does not exist


%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_06_exp_bias_allvis_part_fullexp_delta', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 2 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp3 = par_tmp.delta_par(strcmp(par_tmp.new_condition, 'preOP'), :);
            chSS_grp2 = par_tmp.delta_par(strcmp(par_tmp.new_condition, 'postOP'), :);
            chSS_grp1 = par_tmp.delta_par(strcmp(par_tmp.new_condition, 'control'), :);
            
            %ttest H0
            [~,p_tt_grp3,~,tt_grp3] = ttest(par_tmp.delta_par(strcmp(par_tmp.new_condition, 'preOP'), :));
            [~,p_tt_grp2,~,tt_grp2] = ttest(par_tmp.delta_par(strcmp(par_tmp.new_condition, 'postOP'), :));
            [~,p_tt_grp1,~,tt_grp1] = ttest(par_tmp.delta_par(strcmp(par_tmp.new_condition, 'control'), :));
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3, chSS_grp1);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3, chSS_grp2);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1)), ', sd = ',sprintf('%0.3f',std(chSS_grp1))));
            disp(strcat('mean preOP   (n=', sprintf('%0.f',length(chSS_grp3)),') = ',sprintf('%0.3f',mean(chSS_grp3)), ', sd = ',sprintf('%0.3f',std(chSS_grp3))));
            disp(strcat('mean postOP  (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2)), ', sd = ',sprintf('%0.3f',std(chSS_grp2))));
            
            disp(strcat('preOP    - ttest : t(',sprintf('%0.0f',tt_grp3.df),') = ',sprintf('%0.3f',tt_grp3.tstat), ', p = ' ,sprintf('%0.3f',p_tt_grp3)));
            disp(strcat('postOP   - ttest : t(',sprintf('%0.0f',tt_grp2.df),') = ',sprintf('%0.3f',tt_grp2.tstat), ', p = ' ,sprintf('%0.3f',p_tt_grp2)));
            disp(strcat('control  - ttest : t(',sprintf('%0.0f',tt_grp1.df),') = ',sprintf('%0.3f',tt_grp1.tstat), ', p = ' ,sprintf('%0.3f',p_tt_grp1)));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3
            
        end %loop part
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'delta_par ~ 1 + visit + part_n + new_condition + part_n:new_condition';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        %one by one analysis
        %indexing every condition of analyses
        %calib vs creativity
        idx_c_cr = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'crea'});
        %calib vs nswitch
        idx_c_ns = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'nswitch'});
        
        part_comp = [{idx_c_cr} ; {idx_c_ns}]; %stock all index
        
        %control vs postop
        idx_c_po = strcmp(par.new_condition, {'control'}) | strcmp(par.new_condition, {'postOP'});
        %control vs preop
        idx_c_pr = strcmp(par.new_condition, {'control'}) | strcmp(par.new_condition, {'preOP'});
        %postop vs preop
        idx_po_pr = strcmp(par.new_condition, {'postOP'}) | strcmp(par.new_condition, {'preOP'});
        
        cond_comp = [{idx_c_po};{idx_c_pr};{idx_po_pr}]; %stock all index
        
        %multiple comparison
        for i_loop_cond = 1: length(cond_comp)
            if i_loop_cond == 1
                diary([savewhere '.txt']);
                disp('---CONTROL vs POSTOP---')
                diary off;
            elseif i_loop_cond == 2
                diary([savewhere '.txt']);
                disp('---CONTROL vs PREOP---')
                diary off;
            elseif i_loop_cond == 3
                diary([savewhere '.txt']);
                disp('---PREOP vs POSTOP---')
                diary off;
            else
            end
            
            for i_loop_part = 1 : length(part_comp)
                par_tmp = par(cond_comp{i_loop_cond} & part_comp{i_loop_part},:); %indexing condition comparison
                if i_loop_part == 1
                    diary([savewhere '.txt']);
                    disp('---CALIB vs CREA---')
                    diary off;
                elseif i_loop_part == 2
                    diary([savewhere '.txt']);
                    disp('---CALIB vs NSWITCH---')
                    diary off;
                else
                end
                % statistics -generalized linear regression model
                % set up formula
                model_formula = 'delta_par ~ 1 + visit + part_n + new_condition + part_n:new_condition';
                % do the analysis
                par_tmp = sortrows(par_tmp,'part_n','ascend'); %put calib part first (as default)
                par_tmp = sortrows(par_tmp,'new_condition','ascend'); %put control first (as default)
                res1 = fitglm(par_tmp, model_formula);
                % write result to text file
                diary([savewhere '.txt']);
                disp(res1);
                disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
                diary off;
            end
        end
        clear par
    end %loop parameter name
    
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.delta_par, 'color', parameters.new_condition);  % define data
    g.geom_point('dodge', 0.6,'alpha',0.4);
    g.stat_summary('geom', {'bar', 'black_errorbar'},'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [1.5 3.5],  'XTick', [2 3], 'XTickLabel', {'CREA', 'NSWITCH'});
    g.set_names('column', '', 'x', 'block', 'y', 'delta_par calib', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    for i_bar = 1 : 9
        set(g.results.stat_summary(i_bar).bar_handle,'FaceAlpha',.7)
    end
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_06_exp_bias_allvis_part_fullexp_delta_onevis', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 2 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp3 = par_tmp.delta_par(strcmp(par_tmp.new_condition, 'preOP'), :);
            chSS_grp2 = par_tmp.delta_par(strcmp(par_tmp.new_condition, 'postOP'), :);
            chSS_grp1 = par_tmp.delta_par(strcmp(par_tmp.new_condition, 'control'), :);
            
            %ttest H0
            [~,p_tt_grp3,~,tt_grp3] = ttest(par_tmp.delta_par(strcmp(par_tmp.new_condition, 'preOP'), :));
            [~,p_tt_grp2,~,tt_grp2] = ttest(par_tmp.delta_par(strcmp(par_tmp.new_condition, 'postOP'), :));
            [~,p_tt_grp1,~,tt_grp1] = ttest(par_tmp.delta_par(strcmp(par_tmp.new_condition, 'control'), :));
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3, chSS_grp1);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3, chSS_grp2);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1)), ', sd = ',sprintf('%0.3f',std(chSS_grp1))));
            disp(strcat('mean preOP   (n=', sprintf('%0.f',length(chSS_grp3)),') = ',sprintf('%0.3f',mean(chSS_grp3)), ', sd = ',sprintf('%0.3f',std(chSS_grp3))));
            disp(strcat('mean postOP  (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2)), ', sd = ',sprintf('%0.3f',std(chSS_grp2))));
            
            disp(strcat('preOP    - ttest : t(',sprintf('%0.0f',tt_grp3.df),') = ',sprintf('%0.3f',tt_grp3.tstat), ', p = ' ,sprintf('%0.3f',p_tt_grp3)));
            disp(strcat('postOP   - ttest : t(',sprintf('%0.0f',tt_grp2.df),') = ',sprintf('%0.3f',tt_grp2.tstat), ', p = ' ,sprintf('%0.3f',p_tt_grp2)));
            disp(strcat('control  - ttest : t(',sprintf('%0.0f',tt_grp1.df),') = ',sprintf('%0.3f',tt_grp1.tstat), ', p = ' ,sprintf('%0.3f',p_tt_grp1)));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3
            
        end %loop part
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'delta_par ~ 1 + visit + part_n + new_condition + part_n:new_condition';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        %one by one analysis
        %indexing every condition of analyses
        %calib vs creativity
        idx_c_cr = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'crea'});
        %calib vs nswitch
        idx_c_ns = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'nswitch'});
        
        part_comp = [{idx_c_cr} ; {idx_c_ns}]; %stock all index
        
        %control vs postop
        idx_c_po = strcmp(par.new_condition, {'control'}) | strcmp(par.new_condition, {'postOP'});
        %control vs preop
        idx_c_pr = strcmp(par.new_condition, {'control'}) | strcmp(par.new_condition, {'preOP'});
        %postop vs preop
        idx_po_pr = strcmp(par.new_condition, {'postOP'}) | strcmp(par.new_condition, {'preOP'});
        
        cond_comp = [{idx_c_po};{idx_c_pr};{idx_po_pr}]; %stock all index
        
        %multiple comparison
        for i_loop_cond = 1: length(cond_comp)
            if i_loop_cond == 1
                diary([savewhere '.txt']);
                disp('---CONTROL vs POSTOP---')
                diary off;
            elseif i_loop_cond == 2
                diary([savewhere '.txt']);
                disp('---CONTROL vs PREOP---')
                diary off;
            elseif i_loop_cond == 3
                diary([savewhere '.txt']);
                disp('---PREOP vs POSTOP---')
                diary off;
            else
            end
            
            for i_loop_part = 1 : length(part_comp)
                par_tmp = par(cond_comp{i_loop_cond} & part_comp{i_loop_part},:); %indexing condition comparison
                if i_loop_part == 1
                    diary([savewhere '.txt']);
                    disp('---CALIB vs CREA---')
                    diary off;
                elseif i_loop_part == 2
                    diary([savewhere '.txt']);
                    disp('---CALIB vs NSWITCH---')
                    diary off;
                else
                end
                % statistics -generalized linear regression model
                % set up formula
                model_formula = 'delta_par ~ 1 + visit + part_n + new_condition + part_n:new_condition';
                % do the analysis
                par_tmp = sortrows(par_tmp,'part_n','ascend'); %put calib part first (as default)
                par_tmp = sortrows(par_tmp,'new_condition','ascend'); %put control first (as default)
                res1 = fitglm(par_tmp, model_formula);
                % write result to text file
                diary([savewhere '.txt']);
                disp(res1);
                disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
                diary off;
            end
        end
        clear par
    end %loop parameter name
    
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.delta_par, 'color', parameters.new_condition);  % define data
    g.geom_point('dodge', 0.6,'alpha',0.4);
    g.stat_summary('geom', {'bar', 'black_errorbar'},'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [1.5 3.5],  'XTick', [2 3], 'XTickLabel', {'CREA', 'NSWITCH'});
    g.set_names('column', '', 'x', 'block', 'y', 'delta calib', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    for i_bar = 1 : 9
        set(g.results.stat_summary(i_bar).bar_handle,'FaceAlpha',.7)
    end
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_06_exp_bias_allvis_part_fullexp_delta_patient_onevis', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 2 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.delta_par(strcmp(par_tmp.type_cond, 'patient'), :);
            chSS_grp1 = par_tmp.delta_par(strcmp(par_tmp.type_cond, 'control'), :);
            
            %ttest H0
            [~,p_tt_grp2,~,tt_grp2] = ttest(par_tmp.delta_par(strcmp(par_tmp.type_cond, 'patient'), :));
            [~,p_tt_grp1,~,tt_grp1] = ttest(par_tmp.delta_par(strcmp(par_tmp.type_cond, 'control'), :));
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1)), ', sd = ',sprintf('%0.3f',std(chSS_grp1))));
            disp(strcat('mean patient (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2)), ', sd = ',sprintf('%0.3f',std(chSS_grp2))));
            
            disp(strcat('patient  - ttest : t(',sprintf('%0.0f',tt_grp2.df),') = ',sprintf('%0.3f',tt_grp2.tstat), ', p = ' ,sprintf('%0.3f',p_tt_grp2)));
            disp(strcat('control  - ttest : t(',sprintf('%0.0f',tt_grp1.df),') = ',sprintf('%0.3f',tt_grp1.tstat), ', p = ' ,sprintf('%0.3f',p_tt_grp1)));
            
            disp(strcat('patient vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            diary off
            clear  chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1
            
        end %loop part
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'delta_par ~ 1 + visit + part_n + type_cond + part_n:type_cond';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        %one by one analysis
        %indexing every condition of analyses
        %calib vs creativity
        idx_c_cr = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'crea'});
        %calib vs nswitch
        idx_c_ns = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'nswitch'});
        
        part_comp = [{idx_c_cr} ; {idx_c_ns}]; %stock all index
        
        %control vs patient
        idx_c_po = strcmp(par.type_cond, {'control'}) | strcmp(par.type_cond, {'patient'});
        
        cond_comp = [{idx_c_po}]; %stock all index
        
        %multiple comparison
        for i_loop_cond = 1: length(cond_comp)
            if i_loop_cond == 1
                diary([savewhere '.txt']);
                disp('---CONTROL vs PATIENT---')
                diary off;
            elseif i_loop_cond == 2
                diary([savewhere '.txt']);
                disp('---CONTROL vs PREOP---')
                diary off;
            elseif i_loop_cond == 3
                diary([savewhere '.txt']);
                disp('---PREOP vs POSTOP---')
                diary off;
            else
            end
            
            for i_loop_part = 1 : length(part_comp)
                par_tmp = par(cond_comp{i_loop_cond} & part_comp{i_loop_part},:); %indexing condition comparison
                if i_loop_part == 1
                    diary([savewhere '.txt']);
                    disp('---CALIB vs CREA---')
                    diary off;
                elseif i_loop_part == 2
                    diary([savewhere '.txt']);
                    disp('---CALIB vs NSWITCH---')
                    diary off;
                else
                end
                % statistics -generalized linear regression model
                % set up formula
                model_formula = 'delta_par ~ 1 + visit + part_n + type_cond + part_n:type_cond';
                % do the analysis
                par_tmp = sortrows(par_tmp,'part_n','ascend'); %put calib part first (as default)
                par_tmp = sortrows(par_tmp,'type_cond','ascend'); %put control first (as default)
                res1 = fitglm(par_tmp, model_formula);
                % write result to text file
                diary([savewhere '.txt']);
                disp(res1);
                disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
                diary off;
            end
        end
        clear par
    end %loop parameter name
    
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.delta_par, 'color', parameters.type_cond);  % define data
    g.geom_point('dodge', 0.6,'alpha',0.4);
    g.stat_summary('geom', {'bar', 'black_errorbar'},'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [1.5 3.5],  'XTick', [2 3], 'XTickLabel', {'CREA', 'NSWITCH'});
    g.set_names('column', '', 'x', 'block', 'y', 'delta calib', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    for i_bar = 1 : 6
        set(g.results.stat_summary(i_bar).bar_handle,'FaceAlpha',.7)
    end
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.delta_par, 'color', parameters.type_cond);  % define data
    %g.geom_point('dodge', 0.6,'alpha',0.4);
    g.stat_boxplot();  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [1.5 3.5],  'XTick', [2 3], 'XTickLabel', {'CREA', 'NSWITCH'});
    g.set_names('column', '', 'x', 'block', 'y', 'delta calib', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere 'boxplot.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_07_exp_bias_allvis_part_fullexp_condition_subvis', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 10;  % 4 plus legend
    n_y_subplots = 5;
    
    
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
    
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.new_condition, 'group', parameters.subvisit);  % define data
    g.stat_summary('geom', {'point', 'line'}, 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_line_options('styles', {':'})
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_08_exp_bias_allvis_part_fullexp_delta_dis_', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    clear par
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        
        n_part = unique(par.part_n);
        for i_pt = 2 : length(n_part)
            
            par_tmp = par(strcmp(par.part_n, n_part{i_pt}), :);
            
            save_1 = sprintf(strcat(savewhere, par_name{i_pn} , '_',n_part{i_pt}));   % construct file name
            if ~exist([save_1 '.png'], 'file')  % only run, when output does not exist
                
                f = figure(1); clf;
                set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
                    'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
                    'PaperUnits','centimeters', 'Visible', 'off');
                
                g = gramm('x', par_tmp.delta_par, 'color', par_tmp.new_condition);  % define data
                g.stat_bin('nbins',6,'geom','bar'); % this does plot the mean point
                g.set_title(strcat(par_name{i_pn},' - delta - ',n_part{i_pt}))
                g.set_names('column', '', 'x', 'delta distribution', 'y', '', 'color', 'group');
                g.set_text_options('base_size',17);
                g.set_color_options('map',cfg.fig.colormap);
                g.draw();
                
                outliers = varfun(@mean ,par_tmp,'InputVariables',{'delta_par'},'GroupingVariables',{'subvisit','new_condition'});
                outliers = sortrows(outliers,'mean_delta_par','ascend');
                
                % write result to text file
                diary([save_1 '.txt']);
                sprintf(strcat(par_name{i_pn},' - delta - ',n_part{i_pt}))
                disp(outliers)
                diary off
                
                % save figure
                export_fig([save_1 '.png'], '-r600', '-nocrop');  % save figure
                close(f);  % close figure
            end %saving
            clear par_tmp
        end %loop part
        clear par
    end %loop parameter name
    
end  % if file does not exist


%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_09_exp_bias_allvis_part_comparison_condition', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    parameters.outliers(:) = 0;
    parameters.quitters(:) = 1;
    
    %indexing outliers
    filename = strcat(cfg.dir.group, cfg.studyname,  '_outliers_list.mat');
    load(filename); % this load choice data
    for i_o = 1 : length(outliers.subvisit)
        parameters.outliers(strcmp(parameters.subvisit,outliers.subvisit{i_o}),:)=1;
    end
    
    %indexing quitters
    idx_nq = parameters.part > 2;
    idx_nq = unique(parameters.subvisit(idx_nq,:));
    for i_q = 1 : length(idx_nq)
        parameters.quitters(strcmp(parameters.subvisit,idx_nq{i_q}),:)=0;
    end
    
    %decribe population
    liste = unique(parameters.subvisit);
    descr = [];
    for i_li = 1 : length(liste)
        subj_descr = parameters(strcmp(parameters.subvisit,liste(i_li)),:);
        descr = [descr; subj_descr(1,:)];
    end
    
    descr_stat.age.mean = varfun(@mean ,descr,'InputVariables','age','GroupingVariables',{'new_condition','outliers','quitters','sex'});
    descr_stat.age.std = varfun(@std ,descr,'InputVariables','age','GroupingVariables',{'new_condition','outliers','quitters','sex'});
    
    descr_stat.hemis.mean = varfun(@mean ,descr,'InputVariables','age','GroupingVariables',{'new_condition','outliers','quitters','LocFront','LocDG'});
    descr_stat.hemis.std = varfun(@std ,descr,'InputVariables','age','GroupingVariables',{'new_condition','outliers','LocFront','LocDG'});
    
    descr_stat.FSS.mean = varfun(@mean ,descr,'InputVariables','FSS_score','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.FSS.std = varfun(@std ,descr,'InputVariables','FSS_score','GroupingVariables',{'new_condition','outliers','quitters'});
    
    descr_stat.HAD_A.mean = varfun(@mean ,descr,'InputVariables','HAD_Ascore','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.HAD_A.std = varfun(@std ,descr,'InputVariables','HAD_Ascore','GroupingVariables',{'new_condition','outliers','quitters'});
    
    descr_stat.HAD_D.mean = varfun(@mean ,descr,'InputVariables','HAD_Dscore','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.HAD_D.std = varfun(@std ,descr,'InputVariables','HAD_Dscore','GroupingVariables',{'new_condition','outliers','quitters'});
    
    descr_stat.STARK.mean = varfun(@mean ,descr,'InputVariables','STARK_score','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.STARK.std = varfun(@std ,descr,'InputVariables','STARK_score','GroupingVariables',{'new_condition','outliers','quitters'});
    
    descr_stat.BIS.mean = varfun(@nanmean ,descr,'InputVariables','BIS_score','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.BIS.std = varfun(@nanstd ,descr,'InputVariables','BIS_score','GroupingVariables',{'new_condition','outliers','quitters'});
    
    
    
    descr = sortrows(descr,'subject_id','ascend');
    descr = descr(:, [1 5:7 9 22 36 37 38] );
    
    list_id = unique(descr.subject_id);
    prepost = [];
    for i_li = 1 : length(list_id)
        subj_descr = descr(strcmp(descr.subject_id,list_id(i_li)),:);
        if height(subj_descr) > 1 && sum(strcmp(subj_descr.condition, 'preOP')) > 0
            prepost = [prepost; subj_descr];
        else
        end
    end
    
    descr_stat.prepost = prepost;
    descr_stat.descr = descr;
    
    save(strcat(cfg.dir.group,'descriptive_statistics.mat'), 'descr_stat');
    
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp3 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'preOP'), :);
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'postOP'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.new_condition, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3, chSS_grp1);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3, chSS_grp2);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean preOP   (n=', sprintf('%0.f',length(chSS_grp3)),') = ',sprintf('%0.3f',mean(chSS_grp3))));
            disp(strcat('mean postOP  (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + visit + outliers + quitters + part_n + new_condition + part_n:new_condition';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
        
    end %loop parameter name
    
    %check if parameters can predict behaviour
    par = parameters(strcmp(parameters.part_n, 'calib'),:);
    par_beta = par(strcmp(par.parameter_name, 'beta'),:);
    par_bias = par(strcmp(par.parameter_name, 'bias'),:);
    par_logk = par(strcmp(par.parameter_name, 'log_k'),:);
    par = par_beta;
    par.beta = par.parameter_value;
    par.bias = par_bias.parameter_value;
    par.logk = par_logk.parameter_value;
    % statistics -generalized linear regression model
    % set up formula
    model_formula_1a = 'outliers ~ 1 + new_condition + visit + beta + bias + logk';
    model_formula_1b = 'outliers ~ 1 + new_condition + visit + sex + age + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
    model_formula_2a = 'quitters ~ 1 + new_condition + visit + beta + bias + logk';
    model_formula_2b = 'quitters ~ 1 + new_condition + visit + sex + age + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
    par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
    par = sortrows(par,'new_condition','ascend'); %put control first (as default)
    res1a = fitglm(par, model_formula_1a);
    res1b = fitglm(par, model_formula_1b);
    par_p = [par(strcmp(par.new_condition,'preOP'),:);par(strcmp(par.new_condition,'postOP'),:)];
    res2a = fitglm(par_p,model_formula_2a);
    res2b = fitglm(par_p,model_formula_2b);
    % write result to text file
    diary([savewhere '.txt']);
    disp('---why outliers ?---');
    disp(res1a);
    disp(strcat('R2 =',sprintf('%0.3f',res1a.Rsquared.Ordinary)));
    disp(res1b);
    disp(strcat('R2 =',sprintf('%0.3f',res1b.Rsquared.Ordinary)));
    disp('---why quitters ?---');
    disp(res2a);
    disp(strcat('R2 =',sprintf('%0.3f',res2a.Rsquared.Ordinary)));
    disp(res2b);
    disp(strcat('R2 =',sprintf('%0.3f',res2b.Rsquared.Ordinary)));
    diary off;
    
    %open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.new_condition, 'lightness', parameters.quitters, 'marker', parameters.outliers);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - short experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [0.5 3.5], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_color_options('chroma',200,'lightness_range',[60 90],'chroma_range',[60 30]);
    g.set_point_options('base_size',6);
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group', 'lightness', 'quitters', 'marker', 'outliers');
    g.set_text_options('base_size',17);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 12];
    g.facet_axes_handles(1,2).YLim = [-3 2];
    g.facet_axes_handles(1,3).YLim = [-12 -1];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end  % if file does not exist


%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_09_exp_bias_allvis_part_comparison_patient', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    
    %keep only first visits
    parameters(parameters.visit >1,:)=[];
    
    parameters.outliers(:) = 0;
    parameters.quitters(:) = 1;
    
    %indexing outliers
    filename = strcat(cfg.dir.group, cfg.studyname,  '_outliers_list.mat');
    load(filename); % this load choice data
    for i_o = 1 : length(outliers.subvisit)
        parameters.outliers(strcmp(parameters.subvisit,outliers.subvisit{i_o}),:)=1;
    end
    
    %indexing quitters
    idx_nq = parameters.part > 2;
    idx_nq = unique(parameters.subvisit(idx_nq,:));
    for i_q = 1 : length(idx_nq)
        parameters.quitters(strcmp(parameters.subvisit,idx_nq{i_q}),:)=0;
    end
    
    %decribe population
    liste = unique(parameters.subvisit);
    descr = [];
    for i_li = 1 : length(liste)
        subj_descr = parameters(strcmp(parameters.subvisit,liste(i_li)),:);
        descr = [descr; subj_descr(1,:)];
    end
    
    descr_stat.age.mean = varfun(@mean ,descr,'InputVariables','age','GroupingVariables',{'new_condition','outliers','quitters','sex'});
    descr_stat.age.std = varfun(@std ,descr,'InputVariables','age','GroupingVariables',{'new_condition','outliers','quitters','sex'});
    
    descr_stat.hemis.mean = varfun(@mean ,descr,'InputVariables','age','GroupingVariables',{'new_condition','outliers','quitters','LocFront','LocDG'});
    descr_stat.hemis.std = varfun(@std ,descr,'InputVariables','age','GroupingVariables',{'new_condition','outliers','LocFront','LocDG'});
    
    descr_stat.FSS.mean = varfun(@mean ,descr,'InputVariables','FSS_score','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.FSS.std = varfun(@std ,descr,'InputVariables','FSS_score','GroupingVariables',{'new_condition','outliers','quitters'});
    
    descr_stat.HAD_A.mean = varfun(@mean ,descr,'InputVariables','HAD_Ascore','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.HAD_A.std = varfun(@std ,descr,'InputVariables','HAD_Ascore','GroupingVariables',{'new_condition','outliers','quitters'});
    
    descr_stat.HAD_D.mean = varfun(@mean ,descr,'InputVariables','HAD_Dscore','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.HAD_D.std = varfun(@std ,descr,'InputVariables','HAD_Dscore','GroupingVariables',{'new_condition','outliers','quitters'});
    
    descr_stat.STARK.mean = varfun(@mean ,descr,'InputVariables','STARK_score','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.STARK.std = varfun(@std ,descr,'InputVariables','STARK_score','GroupingVariables',{'new_condition','outliers','quitters'});
    
    descr_stat.BIS.mean = varfun(@nanmean ,descr,'InputVariables','BIS_score','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.BIS.std = varfun(@nanstd ,descr,'InputVariables','BIS_score','GroupingVariables',{'new_condition','outliers','quitters'});
    
    
    
    descr = sortrows(descr,'subject_id','ascend');
    descr = descr(:, [1 5:7 9 22 36 37 38] );
    
    descr_stat.descr = descr;
    
    save(strcat(cfg.dir.group,'descriptive_statistics_patient.mat'), 'descr_stat');
    
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'patient'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2, chSS_grp1,'Tail','right');
            [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2,chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean patient (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + outliers + quitters + part_n*type_cond';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
        
    end %loop parameter name
    
    %focus on patients
    par_patient = parameters(strcmp(parameters.type_cond,'patient'),:);
    list_patient = unique(par_patient.subvisit);
    par_u_tot = [];
    for i_l = 1 : length(list_patient)
        par_u = par_patient(strcmp(par_patient.subvisit,list_patient{i_l}),:);
        par_u = par_u(1,:);
        par_u_tot = [par_u_tot;par_u];
    end
    par_patient = par_u_tot;
    % statistics -generalized linear regression model
    % set up formula
    model_formula_1_out  = 'outliers ~ 1 + new_condition + sex + age + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
    model_formula_2_quit = 'quitters ~ 1 + new_condition + sex + age + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
    par_patient = sortrows(par_patient,'part_n','ascend'); %put calib part first (as default)
    par_patient = sortrows(par_patient,'new_condition','ascend'); %put preOP first (as default)
    res1out = fitglm(par_patient, model_formula_1_out);
    res1quit = fitglm(par_patient, model_formula_2_quit);
    % write result to text file
    diary([savewhere '.txt']);
    disp('---why outliers ?---');
    disp(res1out);
    disp(strcat('R2 =',sprintf('%0.3f',res1out.Rsquared.Ordinary)));
    disp('---why quitters ?---');
    disp(res1quit);
    disp(strcat('R2 =',sprintf('%0.3f',res1quit.Rsquared.Ordinary)));
    diary off;
    
    %open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.type_cond, 'lightness', parameters.quitters, 'marker', parameters.outliers);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - short experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [0.5 3.5], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_color_options('chroma',200,'lightness_range',[60 90],'chroma_range',[60 30]);
    g.set_point_options('base_size',6);
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group', 'lightness', 'quitters', 'marker', 'outliers');
    g.set_text_options('base_size',17);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 12];
    g.facet_axes_handles(1,2).YLim = [-3 2];
    g.facet_axes_handles(1,3).YLim = [-12 -1];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_09_exp_bias_allvis_cleaned_comparison_patient', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    
    %exclusion criteria
    parameters = parameters(strcmp(parameters.exclusion_reason,'RAS'),:);
    
    %keep only first visits
    parameters(parameters.visit >1,:)=[];
    
    parameters.outliers(:) = 0;
    parameters.quitters(:) = 1;
    
    %indexing outliers
    filename = strcat(cfg.dir.group, cfg.studyname,  '_outliers_list.mat');
    load(filename); % this load choice data
    for i_o = 1 : length(outliers.subvisit)
        parameters.outliers(strcmp(parameters.subvisit,outliers.subvisit{i_o}),:)=1;
    end
    
    %indexing quitters
    idx_nq = parameters.part > 2;
    idx_nq = unique(parameters.subvisit(idx_nq,:));
    for i_q = 1 : length(idx_nq)
        parameters.quitters(strcmp(parameters.subvisit,idx_nq{i_q}),:)=0;
    end
    
    %decribe population
    liste = unique(parameters.subvisit);
    descr = [];
    for i_li = 1 : length(liste)
        subj_descr = parameters(strcmp(parameters.subvisit,liste(i_li)),:);
        descr = [descr; subj_descr(1,:)];
    end
    
    descr_stat.age.mean = varfun(@mean ,descr,'InputVariables','age','GroupingVariables',{'new_condition','outliers','quitters','sex'});
    descr_stat.age.std = varfun(@std ,descr,'InputVariables','age','GroupingVariables',{'new_condition','outliers','quitters','sex'});
    
    descr_stat.hemis.mean = varfun(@mean ,descr,'InputVariables','age','GroupingVariables',{'new_condition','outliers','quitters','LocFront','LocDG'});
    descr_stat.hemis.std = varfun(@std ,descr,'InputVariables','age','GroupingVariables',{'new_condition','outliers','LocFront','LocDG'});
    
    descr_stat.FSS.mean = varfun(@mean ,descr,'InputVariables','FSS_score','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.FSS.std = varfun(@std ,descr,'InputVariables','FSS_score','GroupingVariables',{'new_condition','outliers','quitters'});
    
    descr_stat.HAD_A.mean = varfun(@mean ,descr,'InputVariables','HAD_Ascore','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.HAD_A.std = varfun(@std ,descr,'InputVariables','HAD_Ascore','GroupingVariables',{'new_condition','outliers','quitters'});
    
    descr_stat.HAD_D.mean = varfun(@mean ,descr,'InputVariables','HAD_Dscore','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.HAD_D.std = varfun(@std ,descr,'InputVariables','HAD_Dscore','GroupingVariables',{'new_condition','outliers','quitters'});
    
    descr_stat.STARK.mean = varfun(@mean ,descr,'InputVariables','STARK_score','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.STARK.std = varfun(@std ,descr,'InputVariables','STARK_score','GroupingVariables',{'new_condition','outliers','quitters'});
    
    descr_stat.BIS.mean = varfun(@nanmean ,descr,'InputVariables','BIS_score','GroupingVariables',{'new_condition','outliers','quitters'});
    descr_stat.BIS.std = varfun(@nanstd ,descr,'InputVariables','BIS_score','GroupingVariables',{'new_condition','outliers','quitters'});
    
    descr_tt = descr;
    descr_tt(descr_tt.outliers>0,:)=[];
    descr_tt(descr_tt.quitters>0,:)=[];
    
    for i_psy = 26 : 30
        chSS_grp2 = table2array(descr_tt(strcmp(descr_tt.type_cond, 'patient'), i_psy));
        chSS_grp1 = table2array(descr_tt(strcmp(descr_tt.type_cond, 'control'), i_psy));
        
        %ttest accordingly
        %[H,P,CI,STATS]
        [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
        [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2,chSS_grp1);
        
        % write result to text file
        diary([savewhere '.txt']);
        disp(strcat('psy dimension: ', descr_tt.Properties.VariableNames{i_psy}));
        disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
        disp(strcat('mean patient (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
        
        disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
        disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
        diary off
        
    end
    
    
    descr = sortrows(descr,'subject_id','ascend');
    descr = descr(:, [1 5:7 9 22 36 37 38] );
    
    descr_stat.descr = descr;
    
    save(strcat(cfg.dir.group,'descriptive_statistics_cleaned_patient.mat'), 'descr_stat');
    
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'patient'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.type_cond, 'control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2, chSS_grp1,'Tail','right');
            [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2,chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean patient (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('patient vs control - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('patient vs control - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('patient vs control - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + outliers + quitters + part_n*type_cond';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
        
    end %loop parameter name
    
    %focus on patients
    par_patient = parameters(strcmp(parameters.type_cond,'patient'),:);
    list_patient = unique(par_patient.subvisit);
    par_u_tot = [];
    for i_l = 1 : length(list_patient)
        par_u = par_patient(strcmp(par_patient.subvisit,list_patient{i_l}),:);
        par_u = par_u(1,:);
        par_u_tot = [par_u_tot;par_u];
    end
    par_patient = par_u_tot;
    % statistics -generalized linear regression model
    % set up formula
    model_formula_1_out  = 'outliers ~ 1 + new_condition + sex + age + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
    model_formula_2_quit = 'quitters ~ 1 + new_condition + sex + age + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
    par_patient = sortrows(par_patient,'part_n','ascend'); %put calib part first (as default)
    par_patient = sortrows(par_patient,'new_condition','ascend'); %put preOP first (as default)
    res1out = fitglm(par_patient, model_formula_1_out);
    res1quit = fitglm(par_patient, model_formula_2_quit);
    % write result to text file
    diary([savewhere '.txt']);
    disp('---why outliers ?---');
    disp(res1out);
    disp(strcat('R2 =',sprintf('%0.3f',res1out.Rsquared.Ordinary)));
    disp('---why quitters ?---');
    disp(res1quit);
    disp(strcat('R2 =',sprintf('%0.3f',res1quit.Rsquared.Ordinary)));
    diary off;
    
    %open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.type_cond, 'lightness', parameters.quitters, 'marker', parameters.outliers);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - short experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [0.5 3.5], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_color_options('chroma',200,'lightness_range',[60 90],'chroma_range',[60 30]);
    g.set_point_options('base_size',6);
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group', 'lightness', 'quitters', 'marker', 'outliers');
    g.set_text_options('base_size',17);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 12];
    g.facet_axes_handles(1,2).YLim = [-3 2];
    g.facet_axes_handles(1,3).YLim = [-12 -1];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end  % if file does not exist


%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_10_exp_bias_allvis_part_postOP_old_new', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    parameters = parameters(strcmp(parameters.new_condition, 'postOP'), :);
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name(i_pn)), :);
        n_part = unique(par.part_n);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.version_test, 'new'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.version_test, 'old'), :);
            
            %ttest accordingly
            [h, p_interaction_lin_of_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean old (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean new (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('ttest - pvalue old vs new = ',sprintf('%0.3f',mean(p_interaction_lin_of_exp_1))));
            
            diary off
            
            save_1 = sprintf(strcat(savewhere, par_name{i_pn} , '_',n_part{i_pt}));   % construct file name
            if ~exist([save_1 '.png'], 'file')  % only run, when output does not exist
                
                f = figure(1); clf;
                set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
                    'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
                    'PaperUnits','centimeters', 'Visible', 'off');
                
                g = gramm('x', par_tmp.parameter_value, 'color', par_tmp.version_test);  % define data
                g.stat_bin('nbins',6,'geom','bar'); % this does plot the mean point
                g.set_title(strcat(par_name{i_pn},' - delta - ',n_part{i_pt}))
                g.set_names('column', '', 'x', 'delta distribution', 'y', '', 'color', 'version test');
                g.set_text_options('base_size',17);
                g.set_color_options('map',cfg.fig.colormap);
                g.draw();
                
                % save figure
                export_fig([save_1 '.png'], '-r600', '-nocrop');  % save figure
                close(f);  % close figure
            end
            
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        clear par
        
        
    end %loop parameter name
    
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.new_condition, 'lightness', parameters.version_test);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - postOP only');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.set_color_options('chroma',200,'lightness_range',[40 80],'chroma_range',[90 30]);
    g.axe_property('XLim',[0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
end  % if file does not exist


%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_11_exp_bias_allvis_part_fullexp_visit', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    parameters(strcmp(parameters.new_condition, 'preOP'), :) = [];
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.parameter_value(par_tmp.visit==2, :);
            chSS_grp1 = par_tmp.parameter_value(par_tmp.visit==1, :);
            
            %ttest accordingly
            [h, p_interaction_lin_of_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean visit 1 (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean visit 2 (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('ttest - pvalue visit 1 vs 2 = ',sprintf('%0.3f',mean(p_interaction_lin_of_exp_1))));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_fit_lm = 'parameter_value ~ 1 + visit + part_n + new_condition + visit:new_condition';
        reslm = fitglm(par,model_fit_lm);
        % write result to text file
        diary([savewhere '.txt']);
        disp(reslm);
        diary off;
        
        clear par
        
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.new_condition, 'lightness', parameters.visit);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 ,'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.set_point_options('base_size',6);
    g.set_color_options('chroma',300,'lightness_range',[60 90],'chroma_range',[90 30]);
    g.axe_property('XLim', [0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'lightness', 'visit', 'color', 'group');
    g.set_text_options('base_size',17);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist



%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_12_exp_bias_allvis_part_fullexp_real_c', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    parameters.real_c = strcat(parameters.new_condition,'_',parameters.condition);
    parameters(strcmp(parameters.real_c, 'postOP_preOP'),:)=[];
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            
            %dispatch according to group
            chSS_grp3 = par_tmp.parameter_value(strcmp(par_tmp.real_c, 'preOP_preOP'), :);
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.real_c, 'postOP_postOP'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.real_c, 'control_control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3, chSS_grp1);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3, chSS_grp2);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control_control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean preOP_preOP     (n=', sprintf('%0.f',length(chSS_grp3)),') = ',sprintf('%0.3f',mean(chSS_grp3))));
            disp(strcat('mean postOP_postOP   (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + visit + part_n + real_c + part_n:real_c';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        %one by one analysis
        %indexing every condition of analyses
        %calib vs creativity
        idx_c_cr = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'crea'});
        %calib vs nswitch
        idx_c_ns = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'nswitch'});
        
        part_comp = [{idx_c_cr} ; {idx_c_ns}]; %stock all index
        
        %control vs postop
        idx_c_po = strcmp(par.real_c, {'control_control'}) | strcmp(par.real_c, {'postOP_postOP'});
        %control vs preop
        idx_c_pr = strcmp(par.real_c, {'control_control'}) | strcmp(par.real_c, {'preOP_preOP'});
        %postop vs preop
        idx_po_pr = strcmp(par.real_c, {'postOP_postOP'}) | strcmp(par.real_c, {'preOP_preOP'});
        
        cond_comp = [{idx_c_po};{idx_c_pr};{idx_po_pr}]; %stock all index
        
        %multiple comparison
        for i_loop_cond = 1: length(cond_comp)
            if i_loop_cond == 1
                diary([savewhere '.txt']);
                disp('---CONTROL vs POSTOP---')
                diary off;
            elseif i_loop_cond == 2
                diary([savewhere '.txt']);
                disp('---CONTROL vs PREOP---')
                diary off;
            elseif i_loop_cond == 3
                diary([savewhere '.txt']);
                disp('---PREOP vs POSTOP---')
                diary off;
            else
            end
            
            for i_loop_part = 1 : length(part_comp)
                par_tmp = par(cond_comp{i_loop_cond} & part_comp{i_loop_part},:); %indexing condition comparison
                if i_loop_part == 1
                    diary([savewhere '.txt']);
                    disp('---CALIB vs CREA---')
                    diary off;
                elseif i_loop_part == 2
                    diary([savewhere '.txt']);
                    disp('---CALIB vs NSWITCH---')
                    diary off;
                else
                end
                % statistics -generalized linear regression model
                % set up formula
                model_formula = 'parameter_value ~ 1 + visit + part_n + real_c + part_n:real_c';
                % do the analysis
                par_tmp = sortrows(par_tmp,'part_n','ascend'); %put calib part first (as default)
                par_tmp = sortrows(par_tmp,'new_condition','ascend'); %put control first (as default)
                res1 = fitglm(par_tmp, model_formula);
                % write result to text file
                diary([savewhere '.txt']);
                disp(res1);
                disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
                diary off;
            end
        end
        clear par
    end %loop parameter name
    
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.real_c, 'marker', parameters.visit);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 ,'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.set_point_options('base_size',6);
    g.axe_property('XLim', [0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group','marker','visit');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist



%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_12_exp_bias_allvis_part_fullexp_real_c_postonly', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    parameters.real_c = strcat(parameters.new_condition,'_',parameters.condition);
    parameters = parameters(strcmp(parameters.real_c, 'postOP_postOP'),:);
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            
            %dispatch according to group
            chSS_grp2 = par_tmp.parameter_value(par_tmp.visit==2, :);
            chSS_grp1 = par_tmp.parameter_value(par_tmp.visit==1, :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean visit 1 (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean visit 2 (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('visit 1 vs 2 - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + visit + part_n  + visit:part_n';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
        
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.real_c, 'lightness', parameters.visit);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 ,'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.set_point_options('base_size',6);
    g.set_color_options('chroma',300,'lightness_range',[60 90],'chroma_range',[90 30]);
    g.axe_property('XLim', [0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'lightness', 'condition');
    g.set_text_options('base_size',17);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_12_exp_bias_allvis_part_fullexp_real_c_controlonly', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    parameters.real_c = strcat(parameters.new_condition,'_',parameters.condition);
    parameters = parameters(strcmp(parameters.real_c, 'control_control'),:);
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            
            %dispatch according to group
            chSS_grp2 = par_tmp.parameter_value(par_tmp.visit==2, :);
            chSS_grp1 = par_tmp.parameter_value(par_tmp.visit==1, :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean visit 1 (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean visit 2 (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('visit 1 vs 2 - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + visit + part_n  + visit:part_n';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
        
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.real_c, 'lightness', parameters.visit);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 ,'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.set_point_options('base_size',6);
    g.set_color_options('chroma',300,'lightness_range',[60 90],'chroma_range',[90 30]);
    g.axe_property('XLim', [0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'lightness', 'condition');
    g.set_text_options('base_size',17);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_12_exp_bias_allvis_part_fullexp_real_c_onevis', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    parameters.real_c = strcat(parameters.new_condition,'_',parameters.condition);
    parameters(strcmp(parameters.real_c, 'postOP_preOP'),:)=[];
    
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
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp3 = par_tmp.parameter_value(strcmp(par_tmp.real_c, 'preOP_preOP'), :);
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.real_c, 'postOP_postOP'), :);
            chSS_grp1 = par_tmp.parameter_value(strcmp(par_tmp.real_c, 'control_control'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp3, chSS_grp1);
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest2(chSS_grp3, chSS_grp2);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean control_control (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean preOP_preOP     (n=', sprintf('%0.f',length(chSS_grp3)),') = ',sprintf('%0.3f',mean(chSS_grp3))));
            disp(strcat('mean postOP_postOP   (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('postOP vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('preOP  vs control - independant ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('preOP  vs postOP  - independant ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + visit + part_n + real_c + part_n:real_c';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'new_condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        %one by one analysis
        %indexing every condition of analyses
        %calib vs creativity
        idx_c_cr = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'crea'});
        %calib vs nswitch
        idx_c_ns = strcmp(par.part_n, {'calib'}) | strcmp(par.part_n, {'nswitch'});
        
        part_comp = [{idx_c_cr} ; {idx_c_ns}]; %stock all index
        
        %control vs postop
        idx_c_po = strcmp(par.real_c, {'control_control'}) | strcmp(par.real_c, {'postOP_postOP'});
        %control vs preop
        idx_c_pr = strcmp(par.real_c, {'control_control'}) | strcmp(par.real_c, {'preOP_preOP'});
        %postop vs preop
        idx_po_pr = strcmp(par.real_c, {'postOP_postOP'}) | strcmp(par.real_c, {'preOP_preOP'});
        
        cond_comp = [{idx_c_po};{idx_c_pr};{idx_po_pr}]; %stock all index
        
        %multiple comparison
        for i_loop_cond = 1: length(cond_comp)
            if i_loop_cond == 1
                diary([savewhere '.txt']);
                disp('---CONTROL vs POSTOP---')
                diary off;
            elseif i_loop_cond == 2
                diary([savewhere '.txt']);
                disp('---CONTROL vs PREOP---')
                diary off;
            elseif i_loop_cond == 3
                diary([savewhere '.txt']);
                disp('---PREOP vs POSTOP---')
                diary off;
            else
            end
            
            for i_loop_part = 1 : length(part_comp)
                par_tmp = par(cond_comp{i_loop_cond} & part_comp{i_loop_part},:); %indexing condition comparison
                if i_loop_part == 1
                    diary([savewhere '.txt']);
                    disp('---CALIB vs CREA---')
                    diary off;
                elseif i_loop_part == 2
                    diary([savewhere '.txt']);
                    disp('---CALIB vs NSWITCH---')
                    diary off;
                else
                end
                % statistics -generalized linear regression model
                % set up formula
                model_formula = 'parameter_value ~ 1 + visit + part_n + real_c + part_n:real_c';
                % do the analysis
                par_tmp = sortrows(par_tmp,'part_n','ascend'); %put calib part first (as default)
                par_tmp = sortrows(par_tmp,'new_condition','ascend'); %put control first (as default)
                res1 = fitglm(par_tmp, model_formula);
                % write result to text file
                diary([savewhere '.txt']);
                disp(res1);
                disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
                diary off;
            end
        end
        clear par
    end %loop parameter name
    
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.real_c);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 ,'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_color_options('map',cfg.fig.colormap);
    g.set_text_options('base_size',17);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_13_exp_bias_allvis_part_pre_post', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    
    
    list_sub = unique(parameters.subject_id);
    par = [];
    for i_ls = 1 : length(list_sub)
        par_sub = parameters(strcmp(parameters.subject_id, list_sub{i_ls}),:);
        if length(unique(par_sub.condition)) > 1
            par = [par ; par_sub];
        else
        end
    end
    parameters = par;
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp3 = par_tmp.parameter_value(strcmp(par_tmp.condition, 'preOP'), :);
            chSS_grp2 = par_tmp.parameter_value(strcmp(par_tmp.condition, 'postOP'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_3,~,stats_exp_3] = ttest(chSS_grp3, chSS_grp2);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean preOP   (n=', sprintf('%0.f',length(chSS_grp3)),') = ',sprintf('%0.3f',mean(chSS_grp3))));
            disp(strcat('mean postOP  (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('preOP  vs postOP  - paired ttest : t(',sprintf('%0.0f',stats_exp_3.df),') = ',sprintf('%0.3f',stats_exp_3.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_3)));
            
            diary off
            clear chSS_grp3 chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 p_interaction_lin_of_exp_2 p_interaction_lin_of_exp_3 par_tmp
            
        end %loop part
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + visit + part_n + condition + part_n:condition';        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'condition','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.condition);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 ,'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_title('exponential + bias - all visit - full experiment');
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.axe_property('XLim', [0.5 3.5],'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_color_options('map',cfg.fig.colormap);
    g.set_text_options('base_size',17);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_14_exp_bias_allvis_part_fullexp_condition_w_out', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_new.mat', cfg.dir.modelling);
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
    
    
    delta_pos = parameters.delta_par > 0;
    parameters.delta_pos(:) = {'delta_neg'};
    parameters.delta_pos(delta_pos,:) = {'delta_pos'};
    parameters(strcmp(parameters.part_n,'calib'),:) = [];
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    
    % gramm
    g = gramm('x', parameters.delta_pos, 'color', parameters.new_condition);  % define data
    g.stat_bin();
    g.facet_grid(parameters.parameter_name, parameters.part_n, 'scale', 'fixed');  % split in subplots
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist



%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_14_exp_bias_subjects', cfg.dir.modelling_subj, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 6 ;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    parameters.part = parameters.part - 1 ; % calib as zero
    list_sub = unique(parameters.subject_id);
    par_name = unique(parameters.parameter_name);
    
    for i_s = 1 : length(list_sub)
        par_sub = parameters(strcmp(parameters.subject_id, list_sub{i_s}), :);
        
        savesub = sprintf('%s%s_choice_14_subj_%s', cfg.dir.modelling_subj, cfg.studyname, list_sub{i_s});   % construct file name
        if ~exist([savesub '.png'], 'file')  % only run, when output does not exist
            % open figure
            f = figure(1); clf;
            set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
                'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
                'PaperUnits','centimeters', 'Visible', 'off');
            
            % gramm
            g = gramm('x', par_sub.part, 'y', par_sub.parameter_value, 'color', par_sub.condition, 'group',par_sub.subvisit);  % define data
            g.stat_summary('geom', {'point','line'}, 'setylim', false, 'type', 'sem');  % this does plot the mean point
            g.set_title(strcat(list_sub{i_s},' - visit(s) : ', sprintf('%0.f',length(unique(par_sub.subvisit)))));
            g.facet_grid('', par_sub.parameter_name, 'scale', 'independent');  % split in subplots
            g.axe_property('XLim', [-0.5 2.5],'XTick', [0:2],'XTickLabel', {'Calib', 'HOC', 'Switch'});
            g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
            g.set_color_options('map',cfg.fig.colormap);
            g.set_text_options('base_size',17);
            g.draw();
            g.facet_axes_handles(1,1).YLim = [0 16];
            g.facet_axes_handles(1,2).YLim = [-3 4];
            g.facet_axes_handles(1,3).YLim = [-12 0.5];
            g.update();
            
            % save figure
            export_fig([savesub '.png'], '-r600', '-nocrop');  % save figure
            close(f);  % close figure
        end
    end %loop subject
end  % if file does not exist

%% EXPONENTIAL BIAS - PATIENT vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_15_exp_bias_corr_cleaned_patient_onevis', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
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
    
    par_corr = parameters(strcmp(parameters.parameter_name,'bias'),:);
    subv = unique(par_corr.subvisit);
    corr = [];
    for i_s = 1:length(subv)
        idx_calib = strcmp(par_corr.subvisit,subv(i_s)) & strcmp(par_corr.part_n,'calib');
        idx_delta = strcmp(par_corr.subvisit,subv(i_s)) & strcmp(par_corr.part_n,'nswitch');
        calib = par_corr.parameter_value(idx_calib,:);
        delta = par_corr.delta_par(idx_delta,:);
        type_cond = par_corr.type_cond(idx_delta,:);
        corr_sub = [calib,delta,type_cond];
        corr = [corr;corr_sub];
    end
    corr=cell2table(corr);
    corr.Properties.VariableNames = {'calib','delta','type_cond'};
    [R,P,RL,RU]= corrcoef(corr.calib,corr.delta);
    diary([savewhere '.txt']);
    disp(strcat('correlation r =', sprintf('%0.3f',R(2)),', p =',sprintf('%0.9f',P(2))));
    diary off
    
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', corr.delta, 'y', corr.calib, 'color', corr.type_cond);  % define data
    g.geom_point();
    g.set_names('column', '', 'x', 'delta bias switch', 'y', 'CALIB bias parameter', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
end  % if file does not exist

%% EXPONENTIAL BIAS - PREOP vs POSTOP vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_16_exp_bias_prepost_delta', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 3;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    parameters = parameters(strcmp(parameters.part_n,'nswitch'),:);
    
    par_vis = [];
    list = unique(parameters.subject_id);
    
    par_pat= parameters(strcmp(parameters.type_cond,'patient'),:);
    par_cont = parameters(strcmp(parameters.type_cond,'control'),:);
    
    
    for i_s = 1 : length(list)
        ch_sub = parameters(strcmp(parameters.subject_id, list{i_s}), :);
        if length(unique(ch_sub.visit)) > 1
            ch_sub(ch_sub.visit > 2,:) = []; %remove visit 3
            if length(unique(ch_sub(strcmp(ch_sub.condition,'preOP'),:).part_n)) == length(unique(ch_sub(strcmp(ch_sub.condition,'postOP'),:).part_n))
                par_vis = [par_vis ; ch_sub];
            else
            end
        else
        end
        
    end
    
    delta_mean = varfun(@nanmean ,par_vis,'InputVariables','delta_par','GroupingVariables',{'subject_id','subvisit','type_cond','condition','visit','part','part_n'});
    delta_mean.delta_par= delta_mean.nanmean_delta_par;
    
    par_pat= delta_mean(strcmp(delta_mean.type_cond,'patient'),:);
    par_cont = delta_mean(strcmp(delta_mean.type_cond,'control'),:);
     ch_prepost = [];
    for i_pt = 1 : 2
        if i_pt == 1
            ch_tmp = par_pat;
            diary([savewhere '.txt']);
            disp('PATIENTS');
            diary off
        elseif i_pt == 2
            ch_tmp = par_cont;
            diary([savewhere '.txt']);
            disp('CONTROLES');
            diary off
        else
        end
        %dispatch according to group
        chSS_grp1 = ch_tmp((ch_tmp.visit ==1), :);
        chSS_grp2 = ch_tmp((ch_tmp.visit ==2), :);
        ch_cmp = chSS_grp1(:,1:10);
        ch_cmp.chSSpre = chSS_grp1.delta_par;
        ch_cmp.chSSpost = chSS_grp2.delta_par;
        ch_prepost = [ch_prepost ; ch_cmp];
        %ttest accordingly
        %[H,P,CI,STATS]
        [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest(chSS_grp2.delta_par, chSS_grp1.delta_par);
        
        % write result to text file
        diary([savewhere '.txt']);
        disp(strcat('part: ', sprintf('%0.f',i_pt)));
        disp(strcat('mean visit 1  (n=', sprintf('%0.f',length(unique(chSS_grp1.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp1.delta_par))));
        disp(strcat('mean visit 2 (n=', sprintf('%0.f',length(unique(chSS_grp2.subvisit))),') = ',sprintf('%0.3f',mean(chSS_grp2.delta_par))));
        
        disp(strcat('test vs retest - repeated ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
        
   
        diary off
        clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
        
    end  
    
    delta_mean.rnvisit(:) = {'pre / test'};
    delta_mean.rnvisit(delta_mean.visit == 2,:)= {'post / retest'};
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', delta_mean.type_cond, 'y', delta_mean.delta_par, 'color', delta_mean.rnvisit);  % define data
    g.geom_point('dodge', 0.6,'alpha',0.4);
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.4,'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.axe_property('XLim', [0.5 2.5],  'XTick', [1:3]);
    g.set_names('column', '', 'x', 'condition', 'y', 'Switch_delta bias', 'color', 'visit');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PATIENT vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_17_exp_bias_allvis_cleaned_frontal_onevis', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    
    %exclusion criteria
    parameters = parameters(strcmp(parameters.exclusion_reason,'RAS'),:);
    parameters_frontomesial = parameters(strcmp(parameters.LocEM,'frontomesial'),:);
    parameters_frontolateral= parameters(strcmp(parameters.LocEM,'frontolateral'),:);
    parameters = [parameters_frontomesial;parameters_frontolateral];
    
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
    
    par_corr = parameters(strcmp(parameters.parameter_name,'bias'),:);
    subv = unique(par_corr.subvisit);
    corr = [];
    for i_s = 1:length(subv)
        idx_calib = strcmp(par_corr.subvisit,subv(i_s)) & strcmp(par_corr.part_n,'calib');
        idx_delta = strcmp(par_corr.subvisit,subv(i_s)) & strcmp(par_corr.part_n,'nswitch');
        calib = par_corr.parameter_value(idx_calib,:);
        delta = par_corr.delta_par(idx_delta,:);
        type_cond = par_corr.type_cond(idx_delta,:);
        corr_sub = [calib,delta,type_cond];
        corr = [corr;corr_sub];
    end
    
    [R,P,RL,RU]= corrcoef(cell2mat(corr(:,1)),cell2mat(corr(:,2)));
    
    % gramm
    g = gramm('x', corr(:,1), 'y', corr(:,2), 'color', corr(:,3));  % define data
    g.geom_point();
    g.set_names('column', '', 'x', 'CALIB', 'y', 'delta switch', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + part_n + LocEM + part_n:LocEM ';
        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'type_cond','ascend'); %put control first (as default)
        resa = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM parameter value')
        disp(resa);
        disp(strcat('R2 =',sprintf('%0.3f',resa.Rsquared.Ordinary)));
        diary off;
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'delta_par ~ 1 + part_n + LocEM + part_n:LocEM ';
        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM delta parameter')
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.delta_par(strcmp(par_tmp.LocEM, 'frontomesial'), :);
            chSS_grp1 = par_tmp.delta_par(strcmp(par_tmp.LocEM, 'frontolateral'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2, chSS_grp1,'Tail','right');
            [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2,chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean frontolateral (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean frontomesial (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('frontolateral vs frontomesial - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('frontolateral vs frontomesial - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('frontolateral vs frontomesial - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
            diary off
            
            
            save_1 = sprintf(strcat(savewhere, '_deltapar_', par_name{i_pn} , '_',num2str(i_pt)));   % construct file name
            if ~exist([save_1 '.png'], 'file')  % only run, when output does not exist
                
                f = figure(1); clf;
                set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
                    'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
                    'PaperUnits','centimeters', 'Visible', 'off');
                
                g = gramm('x', par_tmp.delta_par, 'color', par_tmp.type_cond);  % define data
                g.stat_bin('nbins',6,'geom','bar'); % this does plot the mean point
                g.set_title(strcat(par_name{i_pn},' - delta - ',num2str(i_pt)))
                if i_pn == 1
                    g.axe_property('XLim', [-10 10]);
                elseif i_pn == 2
                    g.axe_property('XLim', [-4 4]);
                elseif i_pn == 3
                    g.axe_property('XLim', [-6 6]);
                else
                end
                g.set_names('column', '', 'x', 'delta distribution', 'y', '', 'color', 'group');
                g.set_text_options('base_size',17);
                g.set_color_options('map',cfg.fig.colormap);
                g.draw();
                
                % save figure
                export_fig([save_1 '.png'], '-r600', '-nocrop');  % save figure
                close(f);  % close figure
                
                save_2 = sprintf(strcat(savewhere, '_parameter_', par_name{i_pn} , '_',num2str(i_pt)));   % construct file name
                
                f = figure(1); clf;
                set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
                    'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
                    'PaperUnits','centimeters', 'Visible', 'off');
                
                g = gramm('x', par_tmp.parameter_value, 'color', par_tmp.LocEM);  % define data
                g.stat_bin('nbins',6,'geom','bar'); % this does plot the mean point
                g.set_title(strcat(par_name{i_pn},' - pardis - ',num2str(i_pt)))
                g.set_names('column', '', 'x', 'parameter distribution', 'y', '', 'color', 'group');
                g.set_text_options('base_size',17);
                g.set_color_options('map',cfg.fig.colormap);
                g.draw();
                
                % save figure
                export_fig([save_2 '.png'], '-r600', '-nocrop');  % save figure
                close(f);  % close figure
                
                
            end
            
            
            % statistics -generalized linear regression model
            % set up formula
            model_formula_lesion = 'delta_par ~ 1 + LocDG + LocEM + vol_preOP';
            model_formula_psysoc = 'delta_par ~ 1 + age + sex + diploma_score + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
            model_formula_treat  = 'delta_par ~ 1 + new_condition + radio + chimio + tt_antiep';
            model_formula_orthoass  = 'delta_par ~ 1 + flu_anx + flu_P + TMT_A + TMT_B + TMT_flex + span_for + span_back ';
            
            % do the analysis
            par_pat = par_tmp(strcmp(par_tmp.type_cond,'patient'),:);
            par_ctrl = par_tmp(strcmp(par_tmp.type_cond,'control'),:);
            par_pat = sortrows(par_pat,'part_n','ascend'); %put calib part first (as default)
            par_pat = sortrows(par_pat,'new_condition','ascend'); %put control first (as default)
            res1 = fitglm(par_pat, model_formula_lesion);
            res2 = fitglm(par_pat, model_formula_psysoc);
            %res2ctrl = fitglm(par_ctrl, model_formula_psysoc);
            res3 = fitglm(par_pat, model_formula_treat);
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
%             disp('CONTROL')
%             disp(res2ctrl);
%            disp(strcat('R2 =',sprintf('%0.3f',res2ctrl.Rsquared.Ordinary)));
            disp('TREATMENT FACTOR')
            disp(res3);
            disp(strcat('R2 =',sprintf('%0.3f',res3.Rsquared.Ordinary)));
            disp('ORTHO ASSESSMENT FACTOR')
            disp(res4);
            disp(strcat('R2 =',sprintf('%0.3f',res4.Rsquared.Ordinary)));
            diary off;
            
            
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
            
        end %loop part
        
        
        clear par
    end %loop parameter name
    
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    
    % gramm
    g = gramm('x', parameters.part, 'y', parameters.parameter_value, 'color', parameters.LocEM);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.set_title('exponential + bias - all visit')
    g.axe_property('XLim', [0.5 3.5], 'XTick', [1:3],'XTickLabel', {'Calib', 'HOC', 'Switch'});
    g.set_names('column', '', 'x', 'block', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [0 10];
    g.facet_axes_handles(1,2).YLim = [-1 2];
    g.facet_axes_handles(1,3).YLim = [-8 -3];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist

%% EXPONENTIAL BIAS - PATIENT vs CONTROL
% =========================================================================
savewhere = sprintf('%s%s_TIME_fatigue_modelling_17_exp_bias_allvis_cleaned_patient_onevis_ns', cfg.dir.modelling_plots, cfg.studyname);   % construct file name

if ~exist([savewhere '.png'], 'file')  % only run, when output does not exist
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 2.5;
    
    
    filename = sprintf('%sVBA_TIME_parameters_fatigue_all_models_allvis_part_outliers_new.mat', cfg.dir.modelling);
    load(filename);  % load parameters
    all_parameters = allparameters_new;
    
    parameters = all_parameters(strcmp(all_parameters.model, 'exponential_plus_bias'), :);
    
    %exclusion criteria
    parameters = parameters(strcmp(parameters.exclusion_reason,'RAS'),:);
    parameters_frontomesial = parameters(strcmp(parameters.LocEM,'frontomesial'),:);
    parameters_frontolateral= parameters(strcmp(parameters.LocEM,'frontolateral'),:);
    parameters = [parameters_frontomesial;parameters_frontolateral];
    
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
    
    par_corr = parameters(strcmp(parameters.parameter_name,'bias'),:);
    subv = unique(par_corr.subvisit);
    corr = [];
    for i_s = 1:length(subv)
        idx_calib = strcmp(par_corr.subvisit,subv(i_s)) & strcmp(par_corr.part_n,'calib');
        idx_delta = strcmp(par_corr.subvisit,subv(i_s)) & strcmp(par_corr.part_n,'nswitch');
        calib = par_corr.parameter_value(idx_calib,:);
        delta = par_corr.delta_par(idx_delta,:);
        type_cond = par_corr.type_cond(idx_delta,:);
        corr_sub = [calib,delta,type_cond];
        corr = [corr;corr_sub];
    end
    
    [R,P,RL,RU]= corrcoef(cell2mat(corr(:,1)),cell2mat(corr(:,2)));
    
    % gramm
    g = gramm('x', corr(:,1), 'y', corr(:,2), 'color', corr(:,3));  % define data
    g.geom_point();
    g.set_names('column', '', 'x', 'CALIB', 'y', 'delta switch', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    
    par_name = unique(parameters.parameter_name);
    for i_pn = 1 : length(par_name)
        par = parameters(strcmp(parameters.parameter_name, par_name{i_pn}), :);
        n_part = unique(par.part);
        
        % write result to text file
        diary([savewhere '.txt']);
        sprintf(par_name{i_pn})
        diary off
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + part_n + LocEM + part_n:LocEM ';
        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'type_cond','ascend'); %put control first (as default)
        resa = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM parameter value')
        disp(resa);
        disp(strcat('R2 =',sprintf('%0.3f',resa.Rsquared.Ordinary)));
        diary off;
        
        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'delta_par ~ 1 + part_n + LocEM + part_n:LocEM ';
        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM delta parameter')
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        
        
        for i_pt = 1 : length(n_part)
            %select the part we want to investigate
            idx_tmp = par.part == i_pt;
            par_tmp = par(idx_tmp, :);
            
            %dispatch according to group
            chSS_grp2 = par_tmp.delta_par(strcmp(par_tmp.LocEM, 'frontomesial'), :);
            chSS_grp1 = par_tmp.delta_par(strcmp(par_tmp.LocEM, 'frontolateral'), :);
            
            %ttest accordingly
            %[H,P,CI,STATS]
            [~, p_interaction_lin_of_exp_1,~,stats_exp_1] = ttest2(chSS_grp2, chSS_grp1);
            [~, p_interaction_lin_of_exp_2,~,stats_exp_2] = ttest2(chSS_grp2, chSS_grp1,'Tail','right');
            [p_wilcoxon,h_wilcoxon,stats_wilcoxon] = ranksum(chSS_grp2,chSS_grp1);
            
            % write result to text file
            diary([savewhere '.txt']);
            disp(strcat('part: ', sprintf('%0.f',i_pt)));
            disp(strcat('mean frontolateral (n=', sprintf('%0.f',length(chSS_grp1)),') = ',sprintf('%0.3f',mean(chSS_grp1))));
            disp(strcat('mean frontomesial (n=', sprintf('%0.f',length(chSS_grp2)),') = ',sprintf('%0.3f',mean(chSS_grp2))));
            
            disp(strcat('frontolateral vs frontomesial - ind. bil. ttest : t(',sprintf('%0.0f',stats_exp_1.df),') = ',sprintf('%0.3f',stats_exp_1.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_1)));
            disp(strcat('frontolateral vs frontomesial - ind. uni. ttest : t(',sprintf('%0.0f',stats_exp_2.df),') = ',sprintf('%0.3f',stats_exp_2.tstat), ', p = ' ,sprintf('%0.3f',p_interaction_lin_of_exp_2)));
            disp(strcat('frontolateral vs frontomesial - wilcoxon rank sum test : h=',sprintf('%0.0f',h_wilcoxon),', p = ' ,sprintf('%0.3f',p_wilcoxon)));
            diary off
            
            
            save_1 = sprintf(strcat(savewhere, '_deltapar_', par_name{i_pn} , '_',num2str(i_pt)));   % construct file name
            if ~exist([save_1 '.png'], 'file')  % only run, when output does not exist
                
                f = figure(1); clf;
                set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
                    'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
                    'PaperUnits','centimeters', 'Visible', 'off');
                
                g = gramm('x', par_tmp.delta_par, 'color', par_tmp.type_cond);  % define data
                g.stat_bin('nbins',6,'geom','bar'); % this does plot the mean point
                g.set_title(strcat(par_name{i_pn},' - delta - ',num2str(i_pt)))
                if i_pn == 1
                    g.axe_property('XLim', [-10 10]);
                elseif i_pn == 2
                    g.axe_property('XLim', [-4 4]);
                elseif i_pn == 3
                    g.axe_property('XLim', [-6 6]);
                else
                end
                g.set_names('column', '', 'x', 'delta distribution', 'y', '', 'color', 'group');
                g.set_text_options('base_size',17);
                g.set_color_options('map',cfg.fig.colormap);
                g.draw();
                
                % save figure
                export_fig([save_1 '.png'], '-r600', '-nocrop');  % save figure
                close(f);  % close figure
                
                save_2 = sprintf(strcat(savewhere, '_parameter_', par_name{i_pn} , '_',num2str(i_pt)));   % construct file name
                
                f = figure(1); clf;
                set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
                    'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
                    'PaperUnits','centimeters', 'Visible', 'off');
                
                g = gramm('x', par_tmp.parameter_value, 'color', par_tmp.LocEM);  % define data
                g.stat_bin('nbins',6,'geom','bar'); % this does plot the mean point
                g.set_title(strcat(par_name{i_pn},' - pardis - ',num2str(i_pt)))
                g.set_names('column', '', 'x', 'parameter distribution', 'y', '', 'color', 'group');
                g.set_text_options('base_size',17);
                g.set_color_options('map',cfg.fig.colormap);
                g.draw();
                
                % save figure
                export_fig([save_2 '.png'], '-r600', '-nocrop');  % save figure
                close(f);  % close figure
                
                
            end
            
            
            % statistics -generalized linear regression model
            % set up formula
            model_formula_lesion = 'delta_par ~ 1 + LocDG + LocEM + vol_preOP';
            model_formula_psysoc = 'delta_par ~ 1 + age + sex + diploma_score + FSS_score + HAD_Ascore + HAD_Dscore + STARK_score + BIS_score';
            model_formula_treat  = 'delta_par ~ 1 + new_condition + radio + chimio + tt_antiep';
            model_formula_orthoass  = 'delta_par ~ 1 + flu_anx + flu_P + TMT_A + TMT_B + TMT_flex + span_for + span_back ';
            
            % do the analysis
            par_pat = par_tmp(strcmp(par_tmp.type_cond,'patient'),:);
            par_ctrl = par_tmp(strcmp(par_tmp.type_cond,'control'),:);
            par_pat = sortrows(par_pat,'part_n','ascend'); %put calib part first (as default)
            par_pat = sortrows(par_pat,'new_condition','ascend'); %put control first (as default)
            res1 = fitglm(par_pat, model_formula_lesion);
            res2 = fitglm(par_pat, model_formula_psysoc);
            %res2ctrl = fitglm(par_ctrl, model_formula_psysoc);
            res3 = fitglm(par_pat, model_formula_treat);
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
%             disp('CONTROL')
%             disp(res2ctrl);
%            disp(strcat('R2 =',sprintf('%0.3f',res2ctrl.Rsquared.Ordinary)));
            disp('TREATMENT FACTOR')
            disp(res3);
            disp(strcat('R2 =',sprintf('%0.3f',res3.Rsquared.Ordinary)));
            disp('ORTHO ASSESSMENT FACTOR')
            disp(res4);
            disp(strcat('R2 =',sprintf('%0.3f',res4.Rsquared.Ordinary)));
            diary off;
            
            
            clear chSS_grp2 chSS_grp1 p_interaction_lin_of_exp_1 par_tmp
            
        end %loop part
        
        
        clear par


        
        % statistics -generalized linear regression model
        % set up formula
        model_formula = 'parameter_value ~ 1 + part_n*LocEM ';
        % do the analysis
        par = sortrows(par,'part_n','ascend'); %put calib part first (as default)
        par = sortrows(par,'type_cond','ascend'); %put control first (as default)
        res1 = fitglm(par, model_formula);
        % write result to text file
        diary([savewhere '.txt']);
        disp('FITGLM')
        disp(res1);
        disp(strcat('R2 =',sprintf('%0.3f',res1.Rsquared.Ordinary)));
        diary off;
        
        clear par
    end %loop parameter name
    
    n_x_subplots = 5;  % 4 plus legend
    n_y_subplots = 1.5;
    % open figure
    f = figure(1); clf;
    set(f, 'color', 'w', 'toolbar', 'none', 'units', 'centimeters', ...
        'position', [0.5 20 cfg.fig.width_per_subplot.*n_x_subplots  cfg.fig.height_per_subplot.*n_y_subplots], 'Papertype', 'usletter', ...
        'PaperUnits','centimeters', 'Visible', 'off');
    parameters = parameters(strcmp(parameters.part_n,'nswitch'),:);
    % gramm
    g = gramm('x', parameters.part_n, 'y', parameters.parameter_value, 'color', parameters.LocEM);  % define data
    g.stat_summary('geom', {'point', 'errorbar'},'dodge',0.5 , 'setylim', false, 'type', 'sem');  % this does plot the mean point
    g.set_point_options('base_size',9);
    g.axe_property('XTickLabel', {'Switch'});
    g.facet_grid('', parameters.parameter_name, 'scale', 'independent');  % split in subplots
    g.set_names('column', '', 'x', '', 'y', '', 'color', 'group');
    g.set_text_options('base_size',17);
    g.set_color_options('map',cfg.fig.colormap);
    g.draw();
    g.facet_axes_handles(1,1).YLim = [3 7.5];
    g.facet_axes_handles(1,2).YLim = [-0.5 1.5];
    g.facet_axes_handles(1,3).YLim = [-8 -4];
    g.update();
    
    % save figure
    export_fig([savewhere '.png'], '-r600', '-nocrop');  % save figure
    close(f);  % close figure
    
    
end  % if file does not exist



%
%
% part_n = unique(par_bias.part_n);
% f3 = figure('Position', fig_position);
% for i_part = 1 : length(part_n)
%     par_tmp = par_bias(strcmp(par_bias.part_n,part_n{i_part}),:);
%     subplot(3, 1, i_part)
%     raincloud_plot(par_beta.parameter_value,'box_on', 1, 'box_dodge', 0, 'box_dodge_amount',0, 'dot_dodge_amount', .3, 'color', cb(i_part*2,:), 'cloud_edge_col', cb(i_part*2,:));
%     title(part_n{i_part});
% end
% box off
%
%
% new_condition = unique(par_bias.new_condition);
% f3 = figure('Position', fig_position);
% for i_part = 1 : length(new_condition)
%     par_tmp = par_bias(strcmp(par_bias.new_condition,new_condition{i_part}),:);
%     subplot(3, 1, i_part)
%     raincloud_plot(par_beta.parameter_value,'box_on', 1, 'box_dodge', 0, 'box_dodge_amount',0, 'dot_dodge_amount', .3, 'color', cb(i_part*2,:), 'cloud_edge_col', cb(i_part*2,:));
%     title(new_condition{i_part});
% end
% box off