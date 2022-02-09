function EEG_Step4_IHTT()
% Interhemispheric transfer time (IHTT) estimation from a group average
% current source density (CD). Calculation of an interval of confidence 
% around the value obtained for IHTT.
%
%   INPUTS: 
%       none
%
%   OUPUTS:
%       none
%
%   SCRIPTS USED: 
%       Visualization_Source_Group.m, Get_confidence_interval_delay.m,
%       peakseek.m
%
%   FILES USED: 
%        paths.study/Subjects/SubName/EEG/Delays/CD_CondNameVisualField_LeftBrainV1V2.mat (Current source densities (pA.m) of each brain vertice, trial and time point for the left brain visual cortex: [#trials x #vertices x #timepoints])
%        paths.study/Subjects/SubName/EEG/Delays/CD_CondNameVisualField_RightBrainV1V2.mat (Current source densities (pA.m) of each brain vertice, trial and time point for the right brain visual cortex: [#trials x #vertices x #timepoints])
%		 paths.study/Study_Group_Results/time_vec.mat (Vector of the time sample of the EEG epochs [1 x #time points])
%        paths.mat, params.mat
%
%   FILES GENERATED: 
%       paths.study/Study_Group_Results/Source_group_V1V2_CondName_1.png (CD plot for each individual subject)
%       paths.study/Study_Group_Results/Source_group_V1V2_CondName_2.png (CD plot for mean of all subjects)
%       paths.study/Study_Group_Results/Source_group_V1V2_CondName_3.png (CD plot for confidence interval estimation)
%       paths.study/Study_Group_Results/Delay_IHTT_CondNameVF.mat (IHTT estimation in ms for CondName: double)
%       paths.study/Study_Group_Results/ConfidenceInterval_IHTT_CondNameVF.mat (standard deviation of the IHTT estimation in ms for CondName: double)
%       (CondName = condition)
%
%   EXTERNAL PACKAGES USED:
%       Statistical and Machine Learning Toolbox from MATLAB
%
%   SYNTAX:
%       EEG_Step4_IHTT()
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

addpath(genpath((paths.scripts)))

sub_list    = params.sub_list;
cond_list   = params.cond_list;

% Visualization of the group source
Visualization_Source_Group(sub_list,cond_list,{'V1V2'},paths.data,paths.group)

% Subdivide the trials of the group of subjects to great a confidence
% interval on the delay estimation
Get_confidence_interval_delay(sub_list,cond_list,{'V1V2'},paths.data,paths.group)

close all
end
