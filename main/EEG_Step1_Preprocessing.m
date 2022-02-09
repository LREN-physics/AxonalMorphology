function EEG_Step1_Preprocessing()
% Preprocessing EEG data with Fieldtrip. It assumes that the data was
% acquired under different blocks/runs and it has different conditions.
%
%   INPUTS: 
%       none
%
%   OUPUTS:
%       none
%
%   SCRIPTS USED: 
%       Preprocessing_FieldTrip.m, Read_trcdata.m, Filter_data.m, 
%       Segment_trials.m, RenameTriggers_FirstTriggerRemoval.m,
%       RenameElectrodesMicromed.m, Interpolate_BadChannels.m,
%       Filter_data.m, Artefact_Removal.m, inputsdlg.m
%
%   FILES USED: 
%       paths.study/Subjects/SubName/EEG/Originals/block_*.TRC (original recording file: micromed .TRC extention)
%		paths.study/Subjects/SubName/EEG/Matlab_Outputs/SubName_timestamps_block_*.mat (structure with the timestamps of the triggers presentation: structure with timestamps_check.Time_Dur-->duration of each stimulus presentation in ms with [# trials x 1])
%       paths.study/Subjects/SubName/EEG/Matlab_Outputs/SubName_input_configuration.mat (structure with the parameters used in the EEG paradigm: structure with inputs.DURATION_STIM-->stimulus duration in seconds)
%       paths.study/Subjects/SubName/EEG/Matlab_Outputs/SubName_output_block_*.mat  (structure with the timestamps and subject responses: structure with out.RESPONSE-->response time with [# trials x 1])
%       paths.mat, params.mat, config_eeg.mat
% 
%   FILES GENERATED:        
%       paths.study/Subjects/SubName/EEG/Preprocessed/SubName_Filtered_CondName.mat (Data from all runs together and filtered: Fieldtrip structure)
%       paths.study/Subjects/SubName/EEG/Preprocessed/SubName_Filtered_CondName_Clean.mat (Data from all runs together, filtered and clean for bad trials and bad channels: Fieldtrip structure)
%       paths.study/Subjects/SubName/EEG/Preprocessed/SubName_Filtered_CondName_Clean_Int.mat (Data from all runs together, filtered, clean for bad trials and bad channels and with bad trials interpolated: Fieldtrip structure) 
%       paths.study/Subjects/SubName/EEG/Preprocessed/SubName_avg_CondName.mat (Data all processed and averaged across trials: Fieldtrip structure)
%       paths.study/Subjects/SubName/EEG/Preprocessed/SubName_trials_CondName.mat (Data all processed, with the trial structure kept: Fieldtrip structure)
%       paths.study/Subjects/SubName/EEG/Preprocessed/dataTrials_for_BST_SubName_CondName.mat (Data all processed, with the trial structure kept, in a format that can be uploaded to Brainstorm software: Fieldtrip structure with only one time vector instead of one time cell for each trial)
%       paths.study/Subjects/SubName/EEG/Preprocessed/Avg_Left.png (Plot of SubName_avg_CondName.mat - visual evoked response)
%       paths.study/Subjects/SubName/EEG/Preprocessed/Avg_TP_Left.png (Plot of SubName_avg_CondName.mat - topographic map)
%       (CondName = condition, SubName = subject's name)
%
%   EXTERNAL PACKAGES USED:
%       Fieldtrip, EEGlab
%    
%   SYNTAX:
%       EEG_Step1_Preprocessing()
%
%	NOTE:
%		The step in Preprocessing_FieldTrip 'readdata' is highly 
%		dependent on the type of data used in this analysis: reads
%		only .TRC files, assumes the electrode names given by the Micromed
%		system used, uses the matlab files created during the EEG paradigm.
%		If the user wants to adapt this scripts for their own data, this
%		step has to be adapted. 
%
% -------------------------------------------------------------------------
% Author: Rita Oliveira
% Email: rita.oliveira.uni@gmail.com
% Laboratory for Research in Neuroimaging (LREN)
% Department of Clinical Neuroscience, Lausanne  University Hospital and University of Lausanne
% Mont-Paisible 16, CH-1011 Lausanne, Switzerland
%
% Last updated: 10/02/2022
% -------------------------------------------------------------------------

clc
clear all
close all

% Load paths and params
load('paths.mat')
load('params.mat')
load('config_eeg.mat')

sub_list    = params.sub_list;
cond_list   = params.cond_list;
cfg         = config_eeg.cfg_pre;

addpath(genpath((paths.scripts)))
addpath(paths.eeglab);
addpath(paths.fieldtrip);
ft_defaults

% Initial preprocessing
Preprocessing_FieldTrip(sub_list,cond_list,'readdata',cfg,paths.data)
Preprocessing_FieldTrip(sub_list,cond_list,'filtering',cfg,paths.data)
Preprocessing_FieldTrip(sub_list,cond_list,'segmentation',cfg,paths.data)

% Cleaning steps - User Input needed
Preprocessing_FieldTrip(sub_list,cond_list,'remove_trials_eog',cfg,paths.data)
addpath([paths.fieldtrip,'/external/eeglab']) % add eeglab path for ICA otherwise it doesn't recognize 
Preprocessing_FieldTrip(sub_list,cond_list,'ica',cfg,paths.data)
Preprocessing_FieldTrip(sub_list,cond_list,'clean_voltage',cfg,paths.data)
Preprocessing_FieldTrip(sub_list,cond_list,'clean_visual',cfg,paths.data)

% Recovering bad channels
Preprocessing_FieldTrip(sub_list,cond_list,'interpolate',cfg,paths.data)

% Average
Preprocessing_FieldTrip(sub_list,cond_list,'average',cfg,paths.data)

% Prepare trials in a format that brainstorm can read
Preprocessing_FieldTrip(sub_list,cond_list,'average_savetrials',cfg,paths.data)
Prepare_Brainstorm(sub_list,cond_list,paths.data)
clear ALLEEG EEG ALLCOM LASTCOM CURRENTSET CURRENTSTUDY PLUGINLIST STUDY eeglabUpdater

% Clean files
for sub=1:length(sub_list)
    sub_name    = sub_list{sub};
    sub_dir  = fullfile(paths.data,sub_name,'EEG','Preprocessed');
    cd(sub_dir)
    files=dir('*run*.mat');
    delete(files.name);
    delete('*_eog*')
    delete('*_eog_ica*')
    delete('*_eog_ica_volt*')
    delete('*rejected*')
    delete('*ica*')
    rmdir('QualityCheck','s')
end
    
cd(fullfile(paths.scripts,'main'))

end
