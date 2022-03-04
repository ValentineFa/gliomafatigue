function RUN_01_FATIGUE_MOTIVATION_v02
%function RUN_01_FATIGUE_MOTIVATION_v02
%
%   This function is the main wrapper script for the study GLIOMA_FATIGUE.
%
%   It runs the following tasks:
%
%       1) Calibration of choices (3x)
%       2) Creativity tasks mixed with choices (4x)
%       3) Rating of rewards
%       4) Rating of efforts
%       5) Choices of combined rewards and efforts
%       6) Grip force (physical fatigue)
%       7) Training of n-switch task
%       8) n-switch mixed with choices (cognitive fatigue)
%
%   Antonius Wiehler <antonius.wiehler@gmail.com>
%   Created:  2017-12-10
%   Valentine Facque <facque.valentine@gmail.com>
%   Modified: 2019-02-28




%% SET CONFIGURATION PARAMETERS
% =========================================================================
clear all; close all; clc; %#ok<CLALL>

pr = [];  % anything parameter related
tm = [];  % anything timing related
tr = [];  % anything trial related

pr.do.choice_calibration   = 1;  % do we want to run the choice calibrations?

pr.do.choice_trial_genCREA = 1;  % do we want to generate choice trials for the CREAfatigue?
pr.do.choice_trial_gen     = 1;  % do we want to generate choice trials for the fatigue experiment?
pr.do.creativity           = 1;  % do we want to run the creativity tasks?

pr.do.rating_reward        = 1;  % do we want to run the reward ratings?
pr.do.rating_efforts       = 1;  % do we want to run the effort ratings?
pr.do.weight_reward_effort = 1;  % do we want to run the yes/no choices about reward effort combinations?

pr.do.grip_effort          = 1;  % do we want to run the grip force task?


pr.do.nswitch_training     = 1;  % do we want to run the n-switch training?
pr.do.nswitch_fatigue      = 1;  % do we want to run the n-switch fatigue task?


pr.debug = 0;  % PTB debug configuration

% output directories
pr.path.root           = pwd;  % save the current directory
pr.path.meta           = 'outputs/meta/';  % path to save meta data
pr.path.structure      = 'outputs/day_structure/';
pr.path.calib          = 'outputs/choice/calibrations/';  % path to save choice data
pr.path.choiceprep     = 'outputs/choice/choiceprep/';
pr.path.mainexp_choice = 'outputs/choice/mainexp/';
pr.path.CREAexp_choice = 'outputs/choice/CREAexp/';
pr.path.nswitch        = 'outputs/nswitch/';
pr.path.rating         = 'outputs/ratings/';
create_missing_directories(pr.path);

% tool directories
pr.path.battery    = 'motiscan-battery';
addpath(genpath(pr.path.battery));


rng('default'); rng('shuffle');  % reset rng to a seed based on current time

[~, pr.hostname]   = system('hostname'); % get the name of the PC
pr.hostname        = deblank(pr.hostname);
pr.scriptname      = mfilename('fullpath');  % save the name of this script

pr.name            = input('Subject initials : ','s');  % ask for intials
pr.subid           = input('Subject number : ');  % ask for sub ID
pr.group = mod(pr.subid, 2) + 1;  % which group is the subject depending on number

pr.timestamp       = datestr(now, 30); % start time
tm.start.all       = GetSecs;

pr.conditions           = {'GLIOMAFATIGUE01control', 'GLIOMAFATIGUE01preOP', 'GLIOMAFATIGUE01postOP'};
pr.condition            = input('Is the subject a control[1] or, are we before[2] or after[3] the glioma operation? Please enter the corresponding number: ');

if pr.condition == 3
    pr.time_since_operation = input('How many days since the operation? Please enter: ');
else
    pr.time_since_operation = nan;  % if we are before the operation, we cannot ask
end

%
pr.date            = extractBefore(pr.timestamp, 'T');
pr.fnameshort      = [num2str(pr.subid) '_' pr.name '_' pr.date];
pr.fname           = [num2str(pr.subid) '_' pr.name '_' pr.timestamp];
pr.savewhere       = [pr.path.meta 'meta_' mfilename '_' pr.fname];

tm.log = [];  % initialize log structure as empty

% trial generation creativity task
% the following is used per session!
pr.choice_tasks.n_trialsCREA(1, 1)       = 50;  % total number of trials per condition, proportion must give integer trials
pr.choice_tasks.n_trialsCREA(1, 2)       = 40;  % total number of trials per condition, proportion must give integer trials
pr.choice_tasks.n_catch_trialsCREA(1, 1) = 1; % how many catch trials for this task in the whole experiment?
pr.choice_tasks.n_catch_trialsCREA(1, 2) = 1; % how many catch trials for this task in the whole experiment?
pr.choice_tasks.n_total_trialsCREA = pr.choice_tasks.n_trialsCREA + pr.choice_tasks.n_catch_trialsCREA;  % total number of trials

pr.n_blocksCREA = 4; % 3 arr�t pendant la t�che CREA

pr.p.main = [0.2 0.6 0.2];  % proportion of trials, must sum to 1
% pr.p.main(1) = low DV trials
% pr.p.main(2) = mid DV trials
% pr.p.main(3) = high DV trials

% choice tasks
pr.choice_tasks.avialable(1, 1) = {'TD_ID'};  % time discounting imm-delayed
pr.choice_tasks.avialable(1, 2) = {'TD_DD'};  % time discounting delayed-delayed
pr.warmup_script_name = 'meta_RUN_01_FATIGUE_MOTIVATION_v01';

% fatigue task
pr.n_blocks                = 23;  % how many blocks do we want to have within each session?
pr.nswitch.switches        = 8;   % how many rule switches in 24 letters?
pr.nswitch.ISI             = 20;   % this is the ITI before we have response data
pr.time.show_task_reminder = 2.5;  % how long to show the reminder before task bloc?
pr.time.wait_task_reminder = 0.1;  % how long to show the reminder before task bloc?

% trial generation fatigue task
% the following is used per session!
pr.choice_tasks.n_trials(1, 1)       = 50;  % total number of trials per condition, proportion must give integer trials  must be divisible by 5
pr.choice_tasks.n_trials(1, 2)       = 40;  % total number of trials per condition, proportion must give integer trials  must be divisible by 4
pr.choice_tasks.n_catch_trials(1, 1) = 1; % how many catch trials for this task in the whole experiment?
pr.choice_tasks.n_catch_trials(1, 2) = 1; % how many catch trials for this task in the whole experiment?
pr.choice_tasks.n_total_trials = pr.choice_tasks.n_trials + pr.choice_tasks.n_catch_trials;  % total number of trials

pr.p.main = [0.2 0.6 0.2];  % proportion of trials, must sum to 1
% pr.p.main(1) = low DV trials
% pr.p.main(2) = mid DV trials
% pr.p.main(3) = high DV trials


save(pr.savewhere, 'pr', 'tr', 'tm');  % save meta config

%% RATING AT BEGIN OF GLIOMA_EXP TEST
% ========================================================================
pr.ptb = O_open_screen_glioma_v01(pr.debug);
save(pr.savewhere, 'pr', 'tr', 'tm');  % save config

% set rating questions
question     = pr.ptb.accent.question(1 : 3);
questionname = {'Fatigue'; 'Stress'; 'Hunger'};

tm.log = O_Rating_GLIOMA_v01(sprintf('gameday_begin_of_test_%03i', 1), pr, tm, question, questionname);  % rating at begin of test
save(pr.savewhere, 'pr', 'tm', 'tr');  % save main exp config

%% CALIBRATION OF CHOICES
% =========================================================================

if pr.do.choice_calibration
    % randomize choice tasks
    % -------------------------------------------------------------------------
    pr.choice_tasks.order = randperm(size(pr.choice_tasks.avialable, 1));  % randomize task order
    
    for i_t = 1 : size(pr.choice_tasks.avialable, 1)
        pr.choice_tasks.suborder(i_t, :) = randperm(size(pr.choice_tasks.avialable(i_t, :), 2));  % randomize sub task order
    end
    
    tr.choice_tasks.finished = zeros(size(pr.choice_tasks.avialable));  % set flag to zero
    save(pr.savewhere, 'pr', 'tr', 'tm');  % save calibration config
    
    % calibrate choice tasks
    % -------------------------------------------------------------------------
    for i_t = 1 : size(pr.choice_tasks.avialable, 1)  % loop tasks
        
        for i_st = 1 : size(pr.choice_tasks.avialable(i_t, :), 2)  % loop subtasks
            
            task = pr.choice_tasks.avialable{pr.choice_tasks.order(i_t), pr.choice_tasks.suborder(i_t, i_st)};  % which task do we want to calibrate
            tm.start.choicetask(pr.choice_tasks.order(i_t), pr.choice_tasks.suborder(i_t, i_st))  = GetSecs; % start time
            
            n_bisections = 3;  % we do this here for a quick calibration
            DM_calibrate_v11(task, pr.name, pr.subid, n_bisections);  % calibrate choice task
            
            tr.choice_tasks.finished(pr.choice_tasks.order(i_t), pr.choice_tasks.suborder(i_t, i_st)) = 1;  % set flag if task is finished
            tm.finish.choicetask(pr.choice_tasks.order(i_t), pr.choice_tasks.suborder(i_t, i_st))  = GetSecs; % end time
            save(pr.savewhere, 'pr', 'tr', 'tm');  % save calibration
            
        end  % loop sub tasks
    end  % loop main tasks
    
end  % calibration choice

%% GENERATE DAY STRUCTURE AND CHOICE TRIALS (CREATIVITY FATIGUE)
% =========================================================================
if pr.do.choice_trial_genCREA
    
    try
        
        
        clear tr
        
        tr.savewhere = sprintf('%s%s_%s_day_structureCREA', pr.path.structure, pr.conditions{pr.condition}, pr.fname);
        
        % ORDER CHOICE TASKS PER SESSION
        % ---------------------------------------------------------------------
        
        % TASK ORDER
        % ---------------------------------------------------------------------
        
        taskorder = [];
        
        for i_t = 1 : pr.n_blocksCREA
            
            % https://fr.mathworks.com/matlabcentral/answers/313282-repeating-a-vector-up-to-a-given-length
            X   = (1 : size(pr.choice_tasks.avialable, 1))';  % this are indices of avialable task domains
            len = sum(pr.choice_tasks.n_total_trialsCREA(:)) / pr.n_blocksCREA;  % this is the target length of trials per block
            extendedX = [repmat(X, floor(len / numel(X)), 1); ...
                X(1:mod(len, numel(X)))];  % is is the repetition
            extendedX = extendedX';
            extendedX = extendedX(randperm(length(extendedX))); % randomize order
            
            taskorder = [taskorder; extendedX];  %#ok<AGROW> % randomize task domains within each block
        end  % for loop blocks
        
        
        % SUB TASK ORDER
        % ---------------------------------------------------------------------
        subtaskorder = [];
        
        for i_t = 1 : size(pr.choice_tasks.avialable, 1)  % loop domains
            
            subindexvector = 1 : size(pr.choice_tasks.avialable(i_t, :), 2);  % build index of subtasks, one element per subtask
            ntrials = pr.choice_tasks.n_total_trialsCREA(i_t, :);  % how many trials per subtask in this domain?
            subindexvector = repelem(subindexvector, ntrials);  % duplicate index vector to match number of tasks
            subindexvector = subindexvector(randperm(length(subindexvector)));  % randomize subtask index
            subtaskorder(taskorder == i_t) = subindexvector; %#ok<AGROW>
            
        end  % for loop domains
        
        subtaskorder = reshape(subtaskorder, pr.n_blocksCREA, []);
        
        
        % JOIN TASK AND SUB TASK ORDER
        % ---------------------------------------------------------------------
        tasks_ordered = cell(pr.n_blocksCREA, sum(pr.choice_tasks.n_total_trialsCREA(:)) ./ pr.n_blocksCREA);
        
        for i_t = 1 : size(pr.choice_tasks.avialable, 1)  % loop task domains
            for i_st = 1 : size(pr.choice_tasks.avialable(i_t, :), 2)  % loop sub tasks
                tasks_ordered(taskorder== i_t & subtaskorder == i_st) = pr.choice_tasks.avialable(i_t, i_st);  % the current task
            end  % sub task
        end  % task domain
        
        tr.choice_tasks{1} = tasks_ordered;  % save ordered tasks
        
        % WHICHTRIALS
        % ---------------------------------------------------------------------
        tasks_ordered = tasks_ordered'; % copy task order
        tasks_ordered = tasks_ordered(:); % make into vector
        whichtrials   = nan(size(tasks_ordered)); % preallocate whichtrials
        
        for i_t = 1 : size(pr.choice_tasks.avialable, 1)  % loop task domains
            for i_st = 1 : size(pr.choice_tasks.avialable(i_t, :), 2)  % loop sub tasks
                
                task                   = pr.choice_tasks.avialable{i_t, i_st};  % the current task
                n_total_trials         = pr.choice_tasks.n_total_trialsCREA(i_t, i_st);  % how many trials of the current task do we want?
                taskindex              = strcmp(tasks_ordered, task);  % find task
                whichtrials(taskindex) = 1 : n_total_trials;  % save increasing trial indix at the right positions
                
            end  % sub task
        end  % task domain
        
        % reshape and transpose to fit into block scheme
        whichtrials = reshape(whichtrials, [], pr.n_blocksCREA);
        whichtrials = whichtrials';
        tr.whichtrials{1} = whichtrials;
        
        % SAVE
        % ---------------------------------------------------------------------
        tr.when = datestr(now);
        save(tr.savewhere, 'tr');
        
    catch
        sca;  % if an error occures during this script, close the screen
        PsychPortAudio('Close', pr.ptb.sound.pahandle);
        psychrethrow(psychlasterror);  % and show last error
    end  % try catch
    
    %% TRIAL GENERATION CHOICE
    % =========================================================================
    
    for i_t = 1 : size(pr.choice_tasks.avialable, 1)  % loop task domains
        for i_st = 1 : size(pr.choice_tasks.avialable(i_t, :), 2)  % loop sub tasks
            
            task           = pr.choice_tasks.avialable{i_t, i_st};  % the current task
            n_total_trials = pr.choice_tasks.n_trialsCREA(i_t, i_st);  % how many trials of the current task do we want?
            n_catch_trials = pr.choice_tasks.n_catch_trialsCREA(i_t, i_st);  % how many trials of the current task do we want?
            
            
            files = dir(strcat(pr.path.calib, task, '_calib_', pr.fnameshort, '*.mat'));
            fname = files(end).name;  % file name of the last calibration
            calib = load(sprintf('%s%s', pr.path.calib, fname));  % load calibration
            
            
            fprintf('Generate game day trials for subject %s, task %s, condition %s and session %02i... (this might take a while.)\n', pr.fname, task, pr.conditions{pr.condition}, 1);  % display to console
            DM_generate_choice_trials_v02(pr.conditions{pr.condition}, 1, pr.path.choiceprep, n_total_trials, n_catch_trials, calib, pr.p, pr.n_blocksCREA);  % generate trials for main experiment
            fprintf('done.\n');
        end  % sub task loop
    end  % task domain loop
    
end  % trial gen


%% GENERATE DAY STRUCTURE AND CHOICE TRIALS (COGNITIVE FATIGUE)
% =========================================================================
if pr.do.choice_trial_gen
    
    
    clear tr
    
    tr.savewhere = sprintf('%s%s_%s_day_structure', pr.path.structure, pr.conditions{pr.condition}, pr.fname);
    
    
    
    % ORDER CHOICE TASKS PER SESSION
    % ---------------------------------------------------------------------
    
    % TASK ORDER
    % ---------------------------------------------------------------------
    
    taskorder = [];
    
    for i_t = 1 : pr.n_blocks
        
        % https://fr.mathworks.com/matlabcentral/answers/313282-repeating-a-vector-up-to-a-given-length
        X   = (1 : size(pr.choice_tasks.avialable, 1))';  % this are indices of avialable task domains
        len = sum(pr.choice_tasks.n_total_trials(:)) / pr.n_blocks;  % this is the target length of trials per block
        extendedX = [repmat(X, floor(len / numel(X)), 1); ...
            X(1:mod(len, numel(X)))];  % is is the repetition
        extendedX = extendedX';
        extendedX = extendedX(randperm(length(extendedX))); % randomize order
        
        taskorder = [taskorder; extendedX];  %#ok<AGROW> % ranomize task domains within each block
    end  % for loop blocks
    
    
    % SUB TASK ORDER
    % ---------------------------------------------------------------------
    subtaskorder = [];
    
    for i_t = 1 : size(pr.choice_tasks.avialable, 1)  % loop domains
        
        subindexvector = 1 : size(pr.choice_tasks.avialable(i_t, :), 2);  % build index of subtasks, one element per subtask
        ntrials = pr.choice_tasks.n_total_trials(i_t, :);  % how many trials per subtask in this domain?
        subindexvector = repelem(subindexvector, ntrials);  % duplicate index vector to match number of tasks
        subindexvector = subindexvector(randperm(length(subindexvector)));  % ranomize subtask index
        subtaskorder(taskorder == i_t) = subindexvector; %#ok<AGROW>
        
    end  % for loop domains
    
    subtaskorder = reshape(subtaskorder, pr.n_blocks, []);
    
    
    % JOIN TASK AND SUB TASK ORDER
    % ---------------------------------------------------------------------
    tasks_ordered = cell(pr.n_blocks, sum(pr.choice_tasks.n_total_trials(:)) ./ pr.n_blocks);
    
    for i_t = 1 : size(pr.choice_tasks.avialable, 1)  % loop task domains
        for i_st = 1 : size(pr.choice_tasks.avialable(i_t, :), 2)  % loop sub tasks
            tasks_ordered(taskorder== i_t & subtaskorder == i_st) = pr.choice_tasks.avialable(i_t, i_st);  % the current task
        end  % sub task
    end  % task domain
    
    
    
    tr.choice_tasks{1} = tasks_ordered;  % save ordered tasks
    
    
    % WHICHTRIALS
    % ---------------------------------------------------------------------
    tasks_ordered = tasks_ordered'; % copy task order
    tasks_ordered = tasks_ordered(:); % make into vector
    whichtrials   = nan(size(tasks_ordered)); % preallocate whichtrials
    
    for i_t = 1 : size(pr.choice_tasks.avialable, 1)  % loop task domains
        for i_st = 1 : size(pr.choice_tasks.avialable(i_t, :), 2)  % loop sub tasks
            
            task                   = pr.choice_tasks.avialable{i_t, i_st};  % the current task
            n_total_trials         = pr.choice_tasks.n_total_trials(i_t, i_st);  % how many trials of the current task do we want?
            taskindex              = strcmp(tasks_ordered, task);  % find task
            whichtrials(taskindex) = 1 : n_total_trials;  % save increasing trial indix at the right positions
            
        end  % sub task
    end  % task domain
    
    % reshape and transpose to fit into block scheme
    whichtrials = reshape(whichtrials, [], pr.n_blocks);
    whichtrials = whichtrials';
    tr.whichtrials{1} = whichtrials;
    
    
    
    % SAVE
    % ---------------------------------------------------------------------
    tr.when = datestr(now);
    save(tr.savewhere, 'tr');
    
    
    
    
    %% TRIAL GENERATION CHOICE
    % =========================================================================
    
    for i_t = 1 : size(pr.choice_tasks.avialable, 1)  % loop task domains
        for i_st = 1 : size(pr.choice_tasks.avialable(i_t, :), 2)  % loop sub tasks
            
            task           = pr.choice_tasks.avialable{i_t, i_st};  % the current task
            n_total_trials = pr.choice_tasks.n_trials(i_t, i_st);  % how many trials of the current task do we want?
            n_catch_trials = pr.choice_tasks.n_catch_trials(i_t, i_st);  % how many trials of the current task do we want?
            
            
            files = dir(strcat(pr.path.calib, task, '_calib_', pr.fnameshort, '*.mat'));
            fname = files(end).name;  % file name of the last calibration
            calib = load(sprintf('%s%s', pr.path.calib, fname));  % load calibration
            
            
            fprintf('Generate game day trials for subject %s, task %s, condition %s and session %02i... (this might take a while.)\n', pr.fname, task, pr.conditions{pr.condition}, 1);  % display to console
            DM_generate_choice_trials_v02(pr.conditions{pr.condition}, 2, pr.path.choiceprep, n_total_trials, n_catch_trials, calib, pr.p, pr.n_blocks);  % generate trials for main experiment using session2 to adress it
            fprintf('done.\n');
            
        end  % sub task loop
    end  % task domain loop
    
end  % trial gen


%% CREATIVITY MIXED WITH CHOICES
%===============================================================================

if pr.do.creativity
    try
        sca;
        % load everything structure related
        % -------------------------------------------------------------------------
        
        %structure       = load(strcat(pr.path.structure, pr.conditions{pr.condition}, '_', pr.fname, '_day_structureCREA.mat'));
        file = dir(sprintf('%s%s_%sT*_day_structureCREA.mat',pr.path.structure, pr.conditions{pr.condition}, pr.fnameshort));
        structure = load([pr.path.structure file(1).name]);
        
        tr.choice_tasks = structure.tr.choice_tasks;  % order of choice tasks
        tr.whichtrials  = structure.tr.whichtrials;  % order of choice tasks
        
        % Follow CREATIVITY task unrollement
        diary CREAdiary     %create a diary to save the command window
        diary on            %start to save inside the diary
        fprintf ('Please proceed to TAC40')
        pr.CREA = input('Is the Creativity task completed ? YES [1] or NO [2]: ');
        if pr.CREA == 2
            pr.CREA = input('Is the Creativity task completed ? YES [1] or NO [2]: ');
        else
            fprintf('Let s continue...')
        end
        
        for blocknum = 1 : pr.n_blocksCREA
            
            tm.log = putLog(tm.log, GetSecs, sprintf('start_session_%03i_%s_block_%03i', 1, pr.conditions{pr.condition}, blocknum));
            fprintf('Running subject %s in session %i, condition %s and block %i\n', pr.fname, 1, pr.conditions{pr.condition}, blocknum);  % print to console for debug
            pr.ptb = O_open_screen_glioma_v01(pr.debug); %open screen
            save(pr.savewhere, 'pr', 'tr', 'tm');  % save config
            % CHOICE TASK MAIN EXP
            % -------------------------------------------------------------------------
            DrawFormattedText(pr.ptb.PTBwindow, [pr.ptb.accent.Tache, 's de choix'], 'center', 'center', pr.ptb.color.white);
            fliptime = Screen('Flip', pr.ptb.PTBwindow);
            WaitSecs(pr.time.show_task_reminder);
            tm.log = putLog(tm.log, fliptime, sprintf('MiniReminderChoice_%03i_%s_block_%03i', 1, pr.conditions{pr.condition}, blocknum));
            tm.log = putLog(tm.log, GetSecs, sprintf('StartChoice_sess_%03i_%s_block_%03i', 1, pr.conditions{pr.condition}, blocknum));
            
            for i_ct = 1 : size(tr.choice_tasks{1}, 2)  % loop the choice trials of this session
                
                task = tr.choice_tasks{1}{blocknum, i_ct};  % copy current task name
                whichtrials = tr.whichtrials{1}(blocknum, i_ct);  % whichtrials
                
                tm.log = putLog(tm.log, GetSecs, sprintf('Start_%s_sess_%03i_%s_block_%03i', task, 1, pr.conditions{pr.condition}, blocknum));
                %%%%%
                [tm.log] = CREA_present_choice_GLIOMA_v01(pr, tm.log, pr.conditions{pr.condition}, 1, blocknum, task, whichtrials);
                %%%%%
                tm.log = putLog(tm.log, GetSecs, sprintf('Finish_%s_sess_%03i_%s_block_%03i', task, 1, pr.conditions{pr.condition}, blocknum));
                
            end   % for loop choice tasks
            tm.log = putLog(tm.log, GetSecs, sprintf('FinishChoice_sess_%03i_%s_block_%03i', 1, pr.conditions{pr.condition}, blocknum));
            
            sca;
            PsychPortAudio('Close', pr.ptb.sound.pahandle);
            
            if blocknum == 1
                fprintf ('Please proceed to TAR first & distant')
                pr.CREA = input('Is the Creativity task completed ? YES [1] or NO [2]: ');
                if pr.CREA == 2
                    pr.CREA = input('Is the Creativity task completed ? YES [1] or NO [2]: ');
                else
                    fprintf('Let s continue...')
                end
            elseif blocknum == 2
                fprintf ('Please proceed to AUT ')
                pr.CREA = input('Is the Creativity task completed ? YES [1] or NO [2]: ');
                if pr.CREA == 2
                    pr.CREA = input('Is the Creativity task completed ? YES [1] or NO [2]: ');
                else
                    fprintf('Let s continue...')
                end
            elseif blocknum == 3
                fprintf ('Please proceed to Analogy')
                pr.CREA = input('Is the Creativity task completed ? YES [1] or NO [2]: ');
                if pr.CREA == 2
                    pr.CREA = input('Is the Creativity task completed ? YES [1] or NO [2]: ');
                else
                    fprintf('Let s continue...')
                end
                
                diary off   % stop the saving inside the diary
                
                pr.ptb = O_open_screen_glioma_v01(pr.debug); %open screen
                save(pr.savewhere, 'pr', 'tr', 'tm');  % save config
            else
            end
            
        end
    catch
        sca;  % if an error occures during this script, close the screen
        
        psychrethrow(psychlasterror);  % and show last error
    end  % try catch
end % if creativty


%% RATING OF REWARDS
% =========================================================================

if pr.do.rating_reward
    cd(pr.path.battery);
    taskRatingR(pr.subid, pr.condition, 'training', 1, 'testing', 1, 'money', 1, 'fullscreen', 1, 'gaze', 0);  % Reward rating
    cd(pr.path.root);
end  % if reward rating


%% RATING OF EFFORTS
% =========================================================================
if pr.do.rating_efforts
    cd(pr.path.battery);
    taskRatingE(pr.subid, pr.condition, 'training', 1, 'testing', 1, 'fullscreen', 1, 'gaze', 0);  % Effort rating
    cd(pr.path.root);
end  % if reward rating



%% CHOICES OF COMBINED REWARDS AND EFFORTS
% =========================================================================
if pr.do.weight_reward_effort
    cd(pr.path.battery);
    taskWeightRE(pr.subid, pr.condition, 'training', 1, 'testing', 1, 'fullscreen', 1, 'gaze', 0);  % Choice combined rewards effort
    cd(pr.path.root);
end  % if reward rating

%% GRIP FORCE (physical fatigue)
% =========================================================================
if pr.do.grip_effort
    cd(pr.path.battery);
    taskGripR(pr.subid, pr.condition, 'calibration', 1, 'training', 1, 'testing', 1, 'rating', 0, 'fmax', [], 'morpho', 0, 'fullscreen', 1, 'gripdevice', 'vernier');
    cd(pr.path.root);
end



%% TRAINING N-SWITCH TASK
% =========================================================================
pr.ptb = O_open_screen_glioma_v01(pr.debug);
save(pr.savewhere, 'pr', 'tr', 'tm');  % save config

if pr.do.nswitch_training
    NS_Training_nswitch_GLIOMA_v01(pr, tm, pr.conditions{pr.condition}); % Training nSwitch
end



%% N-SWITCHES MIXED WITH CHOICES (COGNITIVE FATIGUE)
% =========================================================================
if pr.do.nswitch_fatigue
    
    try
        
        pr.ptb = O_open_screen_glioma_v01(pr.debug);
        save(pr.savewhere, 'pr', 'tr', 'tm');  % save config
        
        
        % set rating questions
        question     = pr.ptb.accent.question(1 : 3);
        questionname = {'Fatigue'; 'Stress'; 'Hunger'};
        
        
        % RATING AT BEGIN OF SESSION
        % -------------------------------------------------------------------------
        
        tm.log = O_Rating_GLIOMA_v01(sprintf('gameday_begin_of_session_%03i', 1), pr, tm, question, questionname);  % rating at begin of tasks
        save(pr.savewhere, 'pr', 'tm', 'tr');  % save main exp config
        
        % load everything structure related
        % -------------------------------------------------------------------------
        %structure       = dir(sprintf('%s%s_%s_day_structure', pr.path.structure, pr.conditions{pr.condition}, pr.fname));
        %structure       = load(strcat(pr.path.structure, pr.conditions{pr.condition}, '_', pr.fnameshort, '_day_structure.mat'));
        
        file = dir(sprintf('%s%s_%sT*_day_structure.mat',pr.path.structure, pr.conditions{pr.condition}, pr.fnameshort));
        structure = load([pr.path.structure file(1).name]);
        tr.choice_tasks = structure.tr.choice_tasks;  % order of choice tasks
        tr.whichtrials  = structure.tr.whichtrials;  % order of choice tasks
        
        % REPEAT BLOCKS
        % -------------------------------------------------------------------------
        
        for blocknum = 1 : pr.n_blocks
            
            tm.log = putLog(tm.log, GetSecs, sprintf('start_session_%03i_%s_block_%03i', 1, pr.conditions{pr.condition}, blocknum));
            fprintf('Running subject %s in session %i, condition %s and block %i\n', pr.fname, 1, pr.conditions{pr.condition}, blocknum);  % print to console for debug
            
            % N SWITCH
            % --------------------------------------------------------------------
            
            tm.log = putLog(tm.log, GetSecs, sprintf('StartNSwitch_sess_%03i_%s_block_%03i', 1, pr.conditions{pr.condition}, blocknum));
            [~, tm.log] = NS_task_nswitch_GLIOMA_v01(pr, tm, pr.conditions{pr.condition}, 1, blocknum, pr.nswitch.ISI, pr.nswitch.switches);
            tm.log = putLog(tm.log, GetSecs, sprintf('FinishNSwitch_sess_%03i_%s_block_%03i', 1, pr.conditions{pr.condition}, blocknum));
            
            % CHOICE TASK MAIN EXP
            % -------------------------------------------------------------------------
            DrawFormattedText(pr.ptb.PTBwindow, [pr.ptb.accent.Tache, 's de choix'], 'center', 'center', pr.ptb.color.white);
            fliptime = Screen('Flip', pr.ptb.PTBwindow);
            WaitSecs(pr.time.show_task_reminder);
            tm.log = putLog(tm.log, fliptime, sprintf('MiniReminderChoice_%03i_%s_block_%03i', 1, pr.conditions{pr.condition}, blocknum));
            tm.log = putLog(tm.log, GetSecs, sprintf('StartChoice_sess_%03i_%s_block_%03i', 1, pr.conditions{pr.condition}, blocknum));
            
            for i_ct = 1 : size(tr.choice_tasks{1}, 2)  % loop the choice trials of this session
                
                task = tr.choice_tasks{1}{blocknum, i_ct};  % copy current task name
                whichtrials = tr.whichtrials{1}(blocknum, i_ct);  % whichtrials
                
                tm.log = putLog(tm.log, GetSecs, sprintf('Start_%s_sess_%03i_%s_block_%03i', task, 1, pr.conditions{pr.condition}, blocknum));
                %%%%%
                [tm.log] = DM_present_choice_glioma_v01(pr, tm.log, pr.conditions{pr.condition}, 2, blocknum, task, whichtrials); %session 2 = choices for Nswitch
                %%%%%
                tm.log = putLog(tm.log, GetSecs, sprintf('Finish_%s_sess_%03i_%s_block_%03i', task, 1, pr.conditions{pr.condition}, blocknum));
                
            end   % for loop choice tasks
            
            tm.log = putLog(tm.log, GetSecs, sprintf('FinishChoice_sess_%03i_%s_block_%03i', 1, pr.conditions{pr.condition}, blocknum));
            
            
            % RATING AT MID OF SESSION
            % -------------------------------------------------------------------------
            if blocknum == floor(pr.n_blocks / 2)
                tm.log = O_Rating_GLIOMA_v01(sprintf('gameday_mid_of_session_%03i', 1), pr, tm, question, questionname);  % rating at begin of tasks
                save(pr.savewhere, 'pr', 'tm', 'tr');  % save main exp config
            end  % if half of experiment
            
            tm.log = putLog(tm.log, GetSecs, sprintf('finish_session_%03i_%s_block_%03i', 1, pr.conditions{pr.condition}, blocknum));
            
            save(pr.savewhere, 'pr', 'tm', 'tr');  % save main exp config
            
        end  % for loop blocks
        
        
        % RATING AT END OF SESSION
        % -------------------------------------------------------------------------
        tm.log = O_Rating_GLIOMA_v01(sprintf('gameday_end_of_session_%03i', 1), pr, tm, question, questionname);  % rating at begin of tasks
        save(pr.savewhere, 'pr', 'tm', 'tr');  % save main exp config
        
        
        sca; % close PTB screen
        PsychPortAudio('Close', pr.ptb.sound.pahandle);
        
    catch
        sca;  % if an error occures during this script, close the screen
        PsychPortAudio('Close', pr.ptb.sound.pahandle);
        psychrethrow(psychlasterror);  % and show last error
    end  % try catch
    
end  % cognitive fatigue

end  % function


