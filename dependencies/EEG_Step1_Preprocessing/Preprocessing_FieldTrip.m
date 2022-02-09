function Preprocessing_FieldTrip(sub_list,cond_list,step,cfg_pre,path_data)
% Preprocessing EEG data with fieldtrip. It assumes that the data was
% acquired under different blocks/runs and it has different conditions.
%
%   INPUTS: 
%       sub_list (cell array)       - list of subject we want to run
%       cond_list (cell array)      - list of conditions used (e.g. Left and
%           Right visual field stimulations)
%       step (string)               - step to run (see below)
%       cfg_pre (structure)         - structure with the information for preprocessing
%       path_data (string)          - path to the subjects data folder 
%
%   OUPUTS:
%       none
%
%   DESCRIPTION OF EACH STEP (and files created):

    % ------------------------  Intial preprocessing ---------------------
    
    % readdata: the raw data is read. It reads the data inside each subject's folder, in a folder called 'Originals', then creates a new one called 'Preprocessed' to put the
    % outputs of this analysis. My raw data is always named like: block_1.trc, block_2.trc, ....
    % If you have different names you need to change the code.
        % sub01_run1.mat (fieldtrip structure with data)
        % sub01_run1_triggers_orig.mat (original times when each trigger happened)
        % sub01_run1_triggers.mat (times when each trigger happened after removing unwanted triggers)
        % sub01_run1_triggers_names_orig.mat (original names of triggers)
        % sub01_run1_triggers_names.mat (good names of triggers after after removing unwanted triggers and renaming them to match my cond_list)
        % ...
        
    % filtering: filtering of each run/block file:
        % sub01_run1_Filtered.mat
        % sub01_run2_Filtered.mat
        % ...
        
    % segmentation: the data is segmented into trials using the trigger information and combine the the trials from the same condition together (example LEFT condition). We stop having one file per each run/block and start having files like:
        % sub01_Filtered_LEFT.mat
        % ...
           
    % ------------------  Cleanig (User input needed) ---------------------
    
    % remove_trials_eog: user can go through each trial and chose which
    % trial to remove by writing it in a pop window
        % sub01_Filtered_LEFT_eog.mat
        % sub01_LEFT_rejectedtrials_eog.mat - to save wich trials were removed (not important)
        % ...
        
    % ica: performes ICA. User needs to chose each component to remove by
    % writing it in a pop window
        % sub01_Filtered_LEFT_eog_ica.mat
        % sub01_LEFT_ica_components_rej.mat - to save wich trials were removed (not important)
        % sub01_LEFT_ica_components.mat - saves the components obtained with ica (not important)
        % ...
    
    % clean_voltage: user needs to choose which trials and electrodes to
    % remove based on the voltage window from fieldtrip
        % sub01_Filtered_LEFT_eog_ica_volt.mat
        % sub01_LEFT_rejectedtrials_visual.mat - to save wich trials were removed (not important)
        % ...
        
    % clean_visual: user can go through each trial and chose which trial to
    % remove by writing it in a pop window. 
        % sub01_Filtered_LEFT_Clean.mat
        % ...
        
    %---------  Last steps (Recovering bad channels and Average) ---------
    
    % interpolate: interpolates the bad channels with the neighbours
        % sub01_Filtered_LEFT_Clean_Int.mat
        % ...

    % average: averages the trials and re-references the data to the
    % average ref
        % sub01_avg_LEFT.mat
        % ...
        
    % average_savetrials: saves the all trials and re-references the data to the
    % average ref
        % sub01_trials_LEFT.mat
        % ...

%   EXAMPLE TO RUN:

    % >> Example cfg_pre:
    % cfg_pre             = [];
    % cfg_pre.sr          = 1024;                   % sampling rate Hz
    % cfg_pre.poststim    = round(0.3*cfg_pre.sr);  % time frames (0.3 is in seconds)
    % cfg_pre.baseline    = round(0.1*cfg_pre.sr);  % time frames (0.1 is in seconds)
    % cfg_pre.nchans      = 130;                    % number of channels that we want (should match the channames)
    % cfg_pre.channames   = channames;              % channames structure
    % cfg_pre.demean      = 'yes';                  % when averaging across trials, remove the mean of the whole trial segment? 
    % cfg_pre.bandpass    = [0.1 40];               % filter frequencies Hz
    % cfg_pre.layout      = layout_study;           % layout structure

    % >> Run steps
    % Preprocessing_FieldTrip(sub_list,cond_list,'readdata',cfg_pre,path_data) 
    % Preprocessing_FieldTrip(sub_list,cond_list,'filtering',cfg_pre,path_data)
    % Preprocessing_FieldTrip(sub_list,cond_list,'segmentation',cfg_pre,path_data)
    % addpath /home/rita/Documents/MATLAB/fieldtrip-20191206/external/eeglab % needed because ICA requires the eeglab in fieldtrip and not the other one external
    % Preprocessing_FieldTrip(sub_list,cond_list,'ica',cfg_pre,path_data)  
    % Preprocessing_FieldTrip(sub_list,cond_list,'clean_voltage',cfg_pre,path_data)  
    % Preprocessing_FieldTrip(sub_list,cond_list,'clean_visual',cfg_pre,path_data) % not necessary, I use it just as a final confirmation that no strange trial are in my data
    % Preprocessing_FieldTrip(sub_list,cond_list,'interpolate',cfg_pre,path_data)
    % Preprocessing_FieldTrip(sub_list,cond_list,'average',cfg_pre,path_data)
    % Preprocessing_FieldTrip(sub_list,cond_list,'average_savetrials',cfg_pre,path_data)

% -------------------------------------------------------------------------
% Author: Rita Oliveira
% Email: rita.oliveira.uni@gmail.com
% Laboratory for Research in Neuroimaging (LREN)
% Department of Clinical Neuroscience, Lausanne  University Hospital and University of Lausanne
% Mont-Paisible 16, CH-1011 Lausanne, Switzerland
%
% Last updated: 28/07/2021
%
%------------- BEGIN CODE --------------


%% Parameters 

sr              = cfg_pre.sr;         % sampling rate
poststim        = cfg_pre.poststim;   % time frames of post-stimulus epoch
baseline        = cfg_pre.baseline;   % time frames of baseline epoch
nchans          = cfg_pre.nchans;     % number of channels
channames       = cfg_pre.channames;  % channel names
cnt.bandpass    = cfg_pre.bandpass;   % band pass filtering (Hz)
layout          = cfg_pre.layout;     % layout of fieldtrip montage
addpath(path_data)

%% Run

% LOOP FOR SUBJECTS
for sub=1:length(sub_list)
    sub_name    = sub_list{sub};
    fprintf('Working on subject %s\n', sub_name);
    
    % define paths for subject
    pathdata    = fullfile(path_data,sub_name,'EEG','Originals');
    pathoutput  = fullfile(path_data,sub_name,'EEG','Preprocessed');
    pathqc      = fullfile(pathoutput,'QualityCheck');
    
    % create output and quality check folders
    if~exist(pathoutput,'dir')
        mkdir(pathoutput)
    end
    if~exist(pathqc,'dir')
        mkdir(pathqc)
    end
    cd(pathdata)
    
    % get names of the eeg files and how many runs/blocks we have
    % WARNING: USER NEEDS TO CHANGE THE NAME TO LOOK FOR, IF IT'S NOT THE SAME AS THIS ONE
    run_list    = dir('block_*.TRC');
    max_run     = length(run_list);
    cd(pathoutput)
    cond = 1;
    
    % LOOP FOR CONDITIONS
    while cond <= length(cond_list) 
        close all
        cond_name  = cond_list{cond};
        if strcmp(step,'preprocessing')
        disp('Preprocessing')
        else
        fprintf('Working on condition %s\n', cond_name);   
        end

        % STEP TO RUN (only runs what is specified in the string step)
        switch step
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'readdata'
           
            fprintf('READING DATA \n')
            cond = length(cond_list)+1; % here we don't have files with different conditions yet, so this is just to avoid doing the loop over the conditions
            
            % LOOP FOR RUNS
            for run=1:max_run
                close all

                % >> Read .TRC data and create Fieldtrip structure
                datafilename        = run_list(run).name;     % data file
                [datastr,TrigsTime,TrigsName]  = Read_trcdata(pathdata,datafilename);

                % WARNING: SPECIFIC TO MY DATA 
                % >> Rename electrodes micromed to be consistent with
                % .elc file and "oficial" names (because micromed system 
                % at chuv has slightly different names)
                datastr = RenameElectrodesMicromed(datastr);
                
                % >> Set num channels names to remove other strange 
                % channels that are sometimes saved.
                cfg                     = [];
                cfg.channel             = channames;
                dataV                   = ft_preprocessing(cfg, datastr);
                
                save(fullfile(pathoutput,[sub_name,'_run',num2str(run),'_triggers_orig','.mat']),'TrigsTime')
                save(fullfile(pathoutput,[sub_name,'_run',num2str(run),'_triggers_names_orig','.mat']),'TrigsName')

                % WARNING: SPECIFIC TO MY DATA 
                % >> Rename triggers (change from numerical info to
                % categorical tag) &&
                % >> First removal of triggers 
                % (due to the experimental set up itself, it's better to
                % remove some triggers that don't correspond to the
                % idealized trial)
                [TrigsName,TrigsTime] = RenameTriggers_FirstTriggerRemoval(TrigsName,TrigsTime,pathdata,sub_name,run);               
                
                
                % Plot
                if run==1
                    % Plot montage
                    figure; ft_plot_layout(layout);
                    saveas(gcf,fullfile(pathqc,'0_layout.png'))
                    
                    % Plot data for some channels
                    figure;
                    plot(datastr.time{1}/60, datastr.trial{1}(15:25,:));
                    legend(datastr.label(15:25));
                    title('Signal for some channels')
                    ylabel('Amplitude')
                    xlabel('Time (min)')
                    saveas(gcf,fullfile(pathqc,'1_Original.png'))
                    
                    % Plot just one channel
                    pos = randi(length(dataV.label));
                    channel = dataV.label(pos);
                    figure; 
                    plot(datastr.time{1}/60,datastr.trial{1}(pos,:))
                    hold on
                    plot(dataV.time{1}/60,dataV.trial{1}(pos,:))
                    legend('normal','filtered')
                    xlabel('Time (min)')
                    ylabel('Amplitude')
                    title(sprintf('Signal for channel %s',channel{1}))
                    legend('normal','demeaned')
                    saveas(gcf,fullfile(pathqc,'2_Normal_Demean.png'))

                    %Plot data for some channels
                    figure;
                    plot(dataV.time{1}/60, dataV.trial{1}(15:25,:));
                    legend(dataV.label(15:25));
                    title('Signal for some channels')
                    ylabel('Amplitude')
                    xlabel('Time (min)')
                    saveas(gcf,fullfile(pathqc,'3_Demean.png'))
                end
                
                % Save
               save(fullfile(pathoutput,[sub_name,'_run',num2str(run),'.mat']),'dataV')
               save(fullfile(pathoutput,[sub_name,'_run',num2str(run),'_triggers','.mat']),'TrigsTime')
               save(fullfile(pathoutput,[sub_name,'_run',num2str(run),'_triggers_names','.mat']),'TrigsName')
            end
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'filtering'

            fprintf('FILTERING \n')
            cond = length(cond_list)+1; % here we don't have files with different conditions yet, so this is just to avoid doing the loop over the conditions
                      
            % LOOP FOR RUNS
            for run=1:max_run
                close all
                
                % Save
                load(fullfile(pathoutput,[sub_name,'_run',num2str(run),'.mat']),'dataV')

                % >> Filtering
                type            = 1; %1=eeglab 2=fieldtrip
                dataVFilt       = Filter_data(dataV, cnt, type);

                % Save
                save(fullfile(pathoutput,[sub_name,'_run',num2str(run),'_Filtered','.mat']),'dataVFilt')
            
                if run==1                  

                    % Plot one channel filtered vs non filtered
                    pos = randi(length(dataV.label));
                    channel = dataV.label(pos);
                    figure; 
                    plot(dataVFilt.time{1}/60,dataV.trial{1}(pos,:))
                    hold on
                    plot(dataVFilt.time{1}/60,dataVFilt.trial{1}(pos,:))
                    legend('normal','filtered')
                    xlabel('Time (min)')
                    ylabel('Amplitude')
                    title(sprintf('Signal for channel %s',channel{1}))
                    saveas(gcf,fullfile(pathqc,'4_Normal_Filtered.png'))

                    % Plot filtered signal 
                    figure;
                    plot(dataVFilt.time{1}/60, dataVFilt.trial{1}(15:25,:));
                    legend(dataVFilt.label(15:25));
                    title('Signal for some channels Filtered')
                    ylabel('Amplitude')
                    xlabel('Time (min)')
                    saveas(gcf,fullfile(pathqc,'5_Filtered.png'))
                end
            end
  
            clear DataEEG Header TrigsTime TrigsName dataV dataVFilt cfg ...
                run_list max_run sub_name pathoutput ALLCOM ALLEEG ...
                CURRENTSET CURRENTSTUDY LASTCOM PLUGINLIST STUDY EEG e eeglabUpdater
            close all
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        case 'segmentation'

            fprintf('SEGMENTING INTO TRIALS \n')
            
            datatojoin = {};
            cd(pathoutput)
            
            % LOOP FOR RUNS
            for run=1:max_run
                close all
                
                % Get data
                dataVFilt   = load(fullfile(pathoutput,[sub_name,'_run',num2str(run),'_Filtered','.mat']),'dataVFilt');
                TrigsTime   = load(fullfile(pathoutput,[sub_name,'_run',num2str(run),'_triggers','.mat']),'TrigsTime');
                TrigsName   = load(fullfile(pathoutput,[sub_name,'_run',num2str(run),'_triggers_names','.mat']),'TrigsName');
                dataVFilt   = dataVFilt.dataVFilt;
                TrigsTime   = TrigsTime.TrigsTime;
                TrigsName   = TrigsName.TrigsName;

                % Take the condition that we want
                Trigs_cond = [];
                for i=1:length(TrigsTime)
                    if strcmp(TrigsName{i},cond_name)
                        Trigs_cond = [Trigs_cond TrigsTime(i)];
                    end
                end

                % >> Segment
                dataVFiltSeg = Segment_trials(dataVFilt,nchans,Trigs_cond,cond,baseline,poststim,sr); % DATA SEGMENTED
                datatojoin{end+1} = dataVFiltSeg;

                % Save
                save(fullfile(pathoutput,[sub_name,'_run',num2str(run),'_Filtered_',cond_name,'.mat']),'dataVFiltSeg')

            end % end runs
            clear dataVFiltSeg

            % >> Join trials all runs
            cfg                      = [];
            cfg.layout               = layout;
            dataVFiltSeg             = ft_appenddata(cfg, datatojoin{:});
            dataVFiltSeg.hdr         = datatojoin{1}.hdr;
            dataVFiltSeg.hdr.nTrials = length(dataVFiltSeg.trial);

            % Save
            save(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'.mat']),'dataVFiltSeg')
            
            % Plot for all trials
            dataVFiltSeg    = rmfield(dataVFiltSeg,'sampleinfo');
            cfg = [];
            cfg.viewmode    = 'butterfly';
            cfg.xlim        = [-baseline/sr, poststim/sr];
            cfg.channel     = {dataVFiltSeg.label{1:length(dataVFiltSeg.label)}};
            ft_databrowser(cfg,dataVFiltSeg)
            saveas(gcf,fullfile(pathqc,['6_Segmented_',cond_name,'.fig']))
            saveas(gcf,fullfile(pathqc,['6_Segmented_',cond_name,'.png']))
                    
            clear cfg dataVFilt dataVFiltSeg datatojoin cond_name ...
                TrigsTime TrigsName 
            close all
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'remove_trials_eog'
            
            % Get data
            dataVFiltSeg   = load(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'.mat']),'dataVFiltSeg');
            dataVFiltSeg   = dataVFiltSeg.dataVFiltSeg;
            
            % >> Remove trials, remove bad channels           
            cfg                 = [];
            cfg.eog             = {'EOG1';'EOG2'};
            cfg.layout          = layout;
            cfg.tasks           = 1;
            aux                 = Artefact_Removal(cfg,dataVFiltSeg);
            dataVFiltSegEOG     = aux{1};
            trials_rej          = aux{2};

            
            % Save
            save(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'_eog','.mat']),'dataVFiltSegEOG')
            save(fullfile(pathoutput,[sub_name,'_',cond_name,'_rejected_trials_eog','.mat']),'trials_rej')
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        case 'ica'

            fprintf('ICA the TRIALS \n')
             
            % Get data
            dataVFiltSegEOG   = load(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'_eog','.mat']),'dataVFiltSegEOG');
            dataVFiltSegEOG   = dataVFiltSegEOG.dataVFiltSegEOG;
            
             % >> Apply ICA      
            cfg                 = [];
            cfg.eeg             = [1:length(dataVFiltSegEOG.label)-2];
            cfg.layout          = layout;
            cfg.tasks           = 4;
            aux                 = Artefact_Removal(cfg,dataVFiltSegEOG);
            dataVFiltSegEOGICA  = aux{1};
            ica_comp            = aux{2};
            comp_rej            = aux{3};
            
            % Save
            save(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'_eog_ica','.mat']),'dataVFiltSegEOGICA')
            save(fullfile(pathoutput,[sub_name,'_',cond_name,'_ica_compoments','.mat']),'ica_comp')
            save(fullfile(pathoutput,[sub_name,'_',cond_name,'_rejected_comp_ica','.mat']),'comp_rej')

            trial=inputdlg('Trial to show:');
            trial=str2num(cell2mat(trial));
            figure; 
            subplot(4,1,1)
            plot(dataVFiltSegEOG.time{1}*10^3,dataVFiltSegEOG.trial{1,trial}(1,:))
            hold on
            plot(dataVFiltSegEOG.time{1}*10^3,dataVFiltSegEOGICA.trial{1,trial}(1,:))
            legend('Before ICA','After ICA')
            xlabel('Time (ms)')
            ylabel('Amplitude')
            title(sprintf('Trial %d for channel FP1',trial))
            subplot(4,1,2)
            plot(dataVFiltSegEOG.time{1}*10^3,dataVFiltSegEOG.trial{1,trial}(3,:))
            hold on
            plot(dataVFiltSegEOG.time{1}*10^3,dataVFiltSegEOGICA.trial{1,trial}(3,:))
            legend('Before ICA','After ICA')
            xlabel('Time (ms)')
            ylabel('Amplitude')
            title(sprintf('Trial %d for channel FP2',trial))
            subplot(4,1,3)
            plot(dataVFiltSegEOG.time{1}*10^3,dataVFiltSegEOG.trial{1,20}(1,:))
            hold on
            plot(dataVFiltSegEOG.time{1}*10^3,dataVFiltSegEOGICA.trial{1,20}(1,:))
            legend('Before ICA','After ICA')
            xlabel('Time (ms)')
            ylabel('Amplitude')
            title('Trial 20 for channel FP1')
            subplot(4,1,4)
            plot(dataVFiltSegEOG.time{1}*10^3,dataVFiltSegEOG.trial{1,20}(3,:))
            hold on
            plot(dataVFiltSegEOG.time{1}*10^3,dataVFiltSegEOGICA.trial{1,20}(3,:))
            legend('Before ICA','After ICA')
            xlabel('Time (ms)')
            ylabel('Amplitude')
            title('Trial 20 for channel FP2')
            saveas(gcf,fullfile(pathqc,['7_Before_After_ICA_',cond_name,'.png']))
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'clean_voltage'

            fprintf('CLEAN BAD CHANNELS AND TRIALS \n') 
            warndlg('INSTRUCTIONS: On the first window please select the trials to be removed and the click quit. On the second window please select the channels to be removed and click quit.')
            
            
            % Get data             
            dataVFiltSegEOGICA   = load(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'_eog_ica.mat']),'dataVFiltSegEOGICA');
            dataVFiltSegEOGICA   = dataVFiltSegEOGICA.dataVFiltSegEOGICA;

            % >> Remove trials, remove bad channels 
            cfg                 = [];
            cfg.eeg             = [1:length(dataVFiltSegEOGICA.label)];
            cfg.layout          = layout;
            cfg.tasks           = 3;
            dataVFiltSegEOGICAVolt = Artefact_Removal(cfg,dataVFiltSegEOGICA);
             
            % Save
            save(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'_eog_ica_volt','.mat']),'dataVFiltSegEOGICAVolt')
          
            clear cfg dataVFiltSegEOGICA dataVFiltSegEOGICAVolt  
            close all
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'clean_visual'
 
            % Get data             
            dataVFiltSegEOGICAVolt   = load(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'_eog_ica_volt.mat']),'dataVFiltSegEOGICAVolt');
            dataVFiltSegEOGICAVolt   = dataVFiltSegEOGICAVolt.dataVFiltSegEOGICAVolt;
            
             % >> Remove trials,checking one by one
            cfg                 = [];
            cfg.eeg             = [1:length(dataVFiltSegEOGICAVolt.label)];
            cfg.layout          = layout;
            cfg.tasks           = 5;
            aux                 = Artefact_Removal(cfg,dataVFiltSegEOGICAVolt);
            dataVFiltSegClean   = aux{1};
            trials_rej          = aux{2};
            
            % Save
            save(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'_Clean','.mat']),'dataVFiltSegClean')
            save(fullfile(pathoutput,[sub_name,'_',cond_name,'_rejected_trials_visual','.mat']),'trials_rej')
 
            % Plot trials without bad channels
            dataVFiltSegClean   = rmfield(dataVFiltSegClean,'sampleinfo');
            cfg = [];
            cfg.viewmode    = 'butterfly';
            cfg.xlim        = [-baseline/sr, poststim/sr];
            cfg.channel     = {dataVFiltSegClean.label{1:length(dataVFiltSegClean.label)}};
            ft_databrowser(cfg,dataVFiltSegClean)
            saveas(gcf,fullfile(pathqc,['8_Cleaned_',cond_name,'.fig']))
            saveas(gcf,fullfile(pathqc,['8_Cleaned_',cond_name,'.png']))

            clear cfg dataVFiltSegICAVolt dataVFiltSegClean  
            close all
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'interpolate'

            fprintf('INTERPOLATE BAD CHANNELS \n')
            
            % Get data
            dataVFiltSegClean   = load(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'_Clean','.mat']),'dataVFiltSegClean');
            dataVFiltSegClean   = dataVFiltSegClean.dataVFiltSegClean;
            dataVFiltSeg        = load(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'.mat']),'dataVFiltSeg');
            dataVFiltSeg        = dataVFiltSeg.dataVFiltSeg;

            % >> Interpolate missing channels
            cfg                         = [];
            cfg.bad_channels            = {};
            cfg.layout                  = layout;
            dataVFiltSegCleanInt        = Interpolate_BadChannels(cfg,dataVFiltSegClean);
            
            % Save
            save(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'_Clean_Int','.mat']),'dataVFiltSegCleanInt')

            % Plot clean interpolated
            dataVFiltSegCleanInt = rmfield(dataVFiltSegCleanInt,'sampleinfo');
            cfg           = [];
            cfg.viewmode  = 'butterfly';
            cfg.xlim      = [-baseline/sr, poststim/sr];
            cfg.channel   = {dataVFiltSegCleanInt.label{1:length(dataVFiltSegCleanInt.label)}};
            ft_databrowser(cfg,dataVFiltSegCleanInt)
            saveas(gcf,fullfile(pathqc,['9_CleanedInterpolated_',cond_name,'.png']))
            saveas(gcf,fullfile(pathqc,['9_CleanedInterpolated_',cond_name,'.fig']))
            
            % Plot channel of interest
            missingchans = {};
            la=dataVFiltSeg.label(1:end-2);
            ld=dataVFiltSegClean.label;
            for i=1:length(la)
                if ~strcmp(la{i},ld)
                    missingchans=[missingchans,{la{i}}];
                end
            end
            if length(missingchans)>0
            pos = randi(length(missingchans)); % chosse a random channel of the interpolated ones
            channel = missingchans(pos);
            idx = find(strcmp(dataVFiltSeg.label,channel));
            figure; 
            plot(dataVFiltSeg.time{1}*10^3,dataVFiltSeg.trial{1,1}(idx,:))
            hold on
            plot(dataVFiltSeg.time{1}*10^3,dataVFiltSegCleanInt.trial{1,1}(idx,:))
            legend('normal','cleaned and interpolated')
            xlabel('Time (ms)')
            ylabel('Amplitude')
            title(sprintf('Trial one for the bad channel %s',channel{1}))
            saveas(gcf,fullfile(pathqc,['10_Normal_CleanedInterpolated_',cond_name,'.png']))
            end
            
            clear cfg dataVFiltSegClean dataVFiltSegCleanInt task
            close all
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'average'

            fprintf('AVERAGING TRIALS \n')
            
            % Get data
            dataVFiltSegCleanInt   = load(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'_Clean_Int','.mat']),'dataVFiltSegCleanInt');
            dataVFiltSegCleanInt   = dataVFiltSegCleanInt.dataVFiltSegCleanInt;

            % Apply average reference
            cfg                     = [];
            cfg.eeg                 = dataVFiltSegCleanInt; %dataint should be the output of the interpolation
            cfg.reref               = 'yes';
            cfg.refchannel          = 'all';
            dataVFiltSegCleanIntRef = ft_preprocessing(cfg,dataVFiltSegCleanInt);

            % >> Average single trials
            cfg                = [];
            cfg.channel        = 'all';
            cfg.trials         = 'all';
            cfg.covariance     = 'no';
            cfg.keeptrials     = 'no';
            cfg.removemean     = 'yes';
            cfg.vartrllength   = 0;
            [dataAverage]      = ft_timelockanalysis(cfg, dataVFiltSegCleanIntRef);

            % >> Demean average
            cfg                 = [];
            cfg.demean          = 'yes';
            dataAverage_demean  = ft_preprocessing(cfg, dataAverage);
            dataAverage         = dataAverage_demean;
            clear dataAverage_deman

            % Save
            save(fullfile(pathoutput,[sub_name,'_avg_',cond_name,'.mat']),'dataAverage')

            % Plot
            cfg             = [];
            cfg.viewmode    = 'butterfly';
            cfg.xlim        = [-baseline/sr, poststim/sr];
            cfg.channel     = {dataAverage.label{1:length(dataAverage.label)}};
            ft_databrowser(cfg,dataAverage)
            saveas(gcf,fullfile(pathqc,['11_Average_',cond_name,'.png']))
            saveas(gcf,fullfile(pathqc,['11_Average_',cond_name,'.fig']))
            saveas(gcf,fullfile(pathoutput,['Avg_',cond_name,'.png']))

            cfg = [];
            cfg.comment = 'xlim';
            cfg.commentpos = 'title';
            cfg.layout = layout;
            interval = linspace(0,0.23,13);
            figure('units','normalized','outerposition',[0 0 1 1]);
            tiledlayout(2,6,'TileSpacing','Compact','Padding','Compact');
            cfg.zlim = [-1.5 1.5];
            for k = 1:length(interval)-1
                nexttile
                cfg.xlim = [interval(k), interval(k+1)];
                ft_topoplotTFR(cfg,dataAverage);
                hold on
            end
            colorbar('horiz','Position',[0.43 0.03 0.16 0.03]);
            sgtitle(['TP plot ',cond_name])
            saveas(gcf,fullfile(pathqc,['12_Average_TP_',cond_name,'.png']))
            saveas(gcf,fullfile(pathoutput,['Avg_TP_',cond_name,'.png']))
            
            clear cfg dataVFiltSegCleanInt dataVFiltSegCleanIntRef ...
                dataAverage task
            close all
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case 'average_savetrials'
             
            % Used to save the individual trials for Brainstorm
            fprintf('AVERAGING TRIALS FOR STATISTICAL PURPOSES \n')
            
            % Get data
            dataVFiltSegCleanInt   = load(fullfile(pathoutput,[sub_name,'_Filtered_',cond_name,'_Clean_Int','.mat']),'dataVFiltSegCleanInt');
            dataVFiltSegCleanInt   = dataVFiltSegCleanInt.dataVFiltSegCleanInt;

            % >> Apply average reference
            cfg                     = [];
            cfg.eeg                 = dataVFiltSegCleanInt; %dataint should be the output of the interpolation
            cfg.reref               = 'yes';
            cfg.refchannel          = 'all';
            dataVFiltSegCleanIntRef = ft_preprocessing(cfg,dataVFiltSegCleanInt);

            % >> Average single trials
            cfg                = [];
            cfg.channel        = 'all';
            cfg.trials         = 'all';
            cfg.covariance     = 'no';
            cfg.keeptrials     = 'yes';
            cfg.removemean     = 'yes';
            cfg.vartrllength   = 0;
            [dataAverage]      = ft_timelockanalysis(cfg, dataVFiltSegCleanIntRef);

            % >> Demean average
            cfg                 = [];
            cfg.demean          = 'yes';
            dataAverage_demean  = ft_preprocessing(cfg, dataAverage);
            dataAverage         = dataAverage_demean;
            clear dataAverage_deman

            % Save
            save(fullfile(pathoutput,[sub_name,'_trials_',cond_name,'.mat']),'dataAverage')                            
                 
        end % END step
        
        cond = cond+1;
        
    end % END conditions
end % END subjects

end % END function
