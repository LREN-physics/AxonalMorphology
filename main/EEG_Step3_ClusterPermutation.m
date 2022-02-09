function EEG_Step3_ClusterPermutation()
% Performs cluster permutation (CP) statistics on the current source densities (CD)
% values obtained within the visual cortical ROIs at the source level
%
%   INPUTS: 
%       none
%
%   OUPUTS:
%       none
%
%   SCRIPTS USED: 
%       Clust_Perm_Source.m, Visualization_ClusterPermutation.m
%
%   FILES USED: 
%       paths.study/Subjects/SubName/EEG/V1V2_vertices.mat
%       paths.study/Subjects/SubName/EEG/cortex.mat
%       paths.study/Subjects/SubName/EEG/Preprocessed/SubName_trials_CondName.mat, 
%       paths.brainstorm_func/SubName/CondName/matrix_scout_*.mat
% 		paths.mat, params.mat, config_eeg.mat      
%
%   FILES GENERATED: 
%       paths.study/Subjects/SubName/EEG/Statistics_Source/Trials_for_stats_Source_SubName_CondName_V1V2.mat (Structure with the trials inputted on the CP analysis)
%       paths.study/Subjects/SubName/EEG/Statistics_Source/Stats_Source_SubName_CondName_V1V2.mat (result of the cluster permutation)
%       paths.study/Subjects/SubName/EEG/Statistics_Source/Neighbours_SubName_CondName_V1V2.mat (neighbour vertices on the cortex)
%       paths.study/Subjects/SubName/EEG/Statistics_Source/GA_activation_SubName_CondName_V1V2.mat (average over trials for each subject of the CD post-stimulus) 
%       paths.study/Subjects/SubName/EEG/Statistics_Source/GA_baseline_SubName_CondName_V1V2.mat  (average over trials for each subject of the CD baseline) 
%       paths.study/Subjects/SubName/EEG/Statistics_Source/GA_diff_SubName_CondName_V1V2.mat (average over trials for each subject of the CD difference between baseline and post-stimulus)   
%       paths.study/Subjects/SubName/EEG/Statistics_Source/separate_complete_L/R_am_SubName_CondName_V1V2.mat (raw average CD during the peirod tested in CP on the left brain V1V2 or right V1V2)
%       paths.study/Subjects/SubName/EEG/Statistics_Source/separate_L/R_am_SubName_CondName_V1V2.mat (significant average CD during the peirod tested in CP on the left brain V1V2 or right V1V2)
%       paths.study/Study_Group_Results/Source_group_CP_V1V2_CondName_1.png (CD plot for each individual subject)
%       paths.study/Study_Group_Results/Source_group_CP_V1V2_CondName_2.png (CD plot for mean of all subjects)
%       paths.study/Study_Group_Results/Source_group_CP_V1V2_CondName_3.png (# significant vertices plot for each individual subject)
%       paths.study/Study_Group_Results/Source_group_CP_V1V2_CondName_4.png (# significant vertices plot for mean of all subjects)
%       (CondName = condition, SubName = subject's name)
%
%   EXTERNAL PACKAGES USED:
%       Fieldtrip
%    
%   SYNTAX:
%       EEG_Step3_ClusterPermutation()
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

% Set paths
load('paths.mat')
load('params.mat')
load('config_eeg.mat')
addpath(paths.fieldtrip);
ft_defaults

sub_list    = params.sub_list;
cond_list   = params.cond_list;
cfg         = config_eeg.cfg_cluster;

addpath(genpath((paths.scripts)))

% Load data and transform the matrices in the format necessary for cluster
% permutation
Clust_Perm_Source(sub_list,cond_list,'prepare',cfg,paths.data,paths.group,paths.brainstorm_func,paths.brainstorm_anat,{'V1V2'})

% Actually run cluster permutation analsyis
Clust_Perm_Source(sub_list,cond_list,'run',cfg,paths.data,paths.group,paths.brainstorm_func,paths.brainstorm_anat,{'V1V2'})

% Plot the results 
Clust_Perm_Source(sub_list,cond_list,'get metrics',cfg,paths.data,paths.group,paths.brainstorm_func,paths.brainstorm_anat,{'V1V2'})

% Plot results for all subjects
Visualization_ClusterPermutation(sub_list,cond_list,{'V1V2'},paths.data,paths.group)

cd(fullfile(paths.scripts,'main'))

end
