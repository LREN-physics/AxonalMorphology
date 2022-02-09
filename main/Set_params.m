function Set_params()
% Creates parameters necessary to run the analysis.
% 1. Defines params.mat with:
%   - the name of the subjects we want to analyse (sub_list)
%   - the name of the EEG conditions used in the stimulation (cond_list)
% 2. Defines config_eeg.mat with:
%   - structure with the configuration parameters for preprocessing (cfg_pre)
%   - structure with the configuration parameters for cluster permutation (cfg_cluster)
%
%   INPUTS: 
%       none
%
%   OUPUTS:
%       none
%
%   SCRIPTS USED: 
%       none
%
%   FILES USED: 
%       paths.mat
%
%   FILES GENERATED: 
%       params.mat
%
%   SYNTAX:
%       Set_params()
%
%   HOW TO USE:
%       Manually edit the script and change the variables to match the
%       ones desired. Then run the script. 
%		If you are not interested in doing the steps:
%		EEG - PREPROCESSING, EEG - SOURCE RECONSTRUCTION or EEG - CLUSTER PERMUTATION 
%	    you can comment lines 64-110
%
% -------------------------------------------------------------------------
% Author: Rita Oliveira
% Email: rita.oliveira.uni@gmail.com
% Laboratory for Research in Neuroimaging (LREN)
% Department of Clinical Neuroscience, Lausanne  University Hospital and University of Lausanne
% Mont-Paisible 16, CH-1011 Lausanne, Switzerland
%
% Last updated: 10/01/2022
% -------------------------------------------------------------------------

clc
clear all
close all
load('paths.mat')

%% PARAMS.MAT 

% -------------------- EDIT HERE --------------------
sub_list                = {'sub_11'};       % (cell array) - list of subject we want to run
cond_list               = {'Left'};         % (cell array) - list of conditions used (in my case I have conditions 'Right' and 'Left', corresponding to stimuli presented on the Left or Right visual field).

% -----------------------------------------------------

% Define structure
params              = [];
params.sub_list     = sub_list;
params.cond_list    = cond_list;

% Save structure
save(fullfile(paths.scripts,'main','params.mat'),'params')

%% CONFIG_EEG.MAT

% Load names of the channels and layout for fieldtrip 
if exist(fullfile(paths.study,'layout_study.mat'),'file')
    load(fullfile(paths.study,'layout_study'))
else
    warndlg('No layout structure exists for your cap. Please save a layout_study.mat in your EEG study folder. See Fieldtrip website for more information')

end

if exist(fullfile(paths.study,'channames.mat'),'file')
    load(fullfile(paths.study,'channames'))
else
    warndlg('No channames structure exists for your cap. Please write and save a vector with the names of the channels in your cap in your EEG study folder')    
end

% -------------------- EDIT HERE --------------------
% Preprocessing configuration structure
cfg_pre             = [];
cfg_pre.sr          = 1024;                   % sampling rate Hz
cfg_pre.poststim    = round(0.3*cfg_pre.sr);  % time frames for post stimulus period (0.3 is in seconds)
cfg_pre.baseline    = round(0.1*cfg_pre.sr);  % time frames for baseline period  (0.1 is in seconds)
cfg_pre.nchans      = 130;                    % number of channels that we want (should match the channames)
cfg_pre.channames   = channames;              % channames structure
cfg_pre.demean      = 'yes';                  % when averaging across trials, remove the mean of the whole trial segment? 
cfg_pre.bandpass    = [0.1 40];               % filter frequencies Hz
cfg_pre.layout      = layout_study;           % layout structure
% -----------------------------------------------------


% -------------------- EDIT HERE --------------------
% Cluster permutation configuration structure
cfg_cluster             = [];                     % initiate structure
cfg_cluster.sr          = 1024;                   % sampling rate Hz
cfg_cluster.baseline    = 98;                     % number of samples of baseline
cfg_cluster.begwindow   = 185;                    % sample corresponding to the beginning of the time window we want to do the analysis (RO: here it corresponds to ~80ms)
cfg_cluster.endwindow   = 315;                    % sample corresponding to the beginning of the time window we want to do the analysis (RO: here it corresponds to ~207ms)
cfg_cluster.layout      = layout_study;           % layout structure
% -----------------------------------------------------

% Define structure
config_eeg              = [];
config_eeg.cfg_pre      = cfg_pre;
config_eeg.cfg_cluster  = cfg_cluster;

% Save structure
save(fullfile(paths.scripts,'main','config_eeg.mat'),'config_eeg')

end
