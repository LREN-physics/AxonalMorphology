function MODEL()
% From G-ratio and interhemispheric transfer time (IHTT) measures, we
% estimate morphological features of the white matter tracts:
% 1) the axonal radius distribution (P(r)) 
% 2) a measure of myelin thickness (G(r)). 
% The way this script is implemented it assumes the tract is the 
% VISUAL TRANSCALLOSAL TRACT connecting the primary and secondary visual 
% areas (V1 and V2) in both hemispheres. It also assumes that the IHTT was
% estimated based on a LEFT visual field stimulation.
%
%   INPUTS: 
%       none
%
%   OUPUTS:
%       none
%
%   SCRIPTS USED: 
%       EstimateMicrostructure_GroupDelay.m
%
%   FILES USED: 
%       paths.study/Study_Group_Results/Delay_IHTT_LeftVF.mat (IHTT estimation in ms for CondName: double)
%		paths.study/Study_Group_Results/ConfidenceInterval_IHTT_LeftVF.mat (standard deviation of the IHTT estimation in ms for CondName: double)
%       paths.study/Subjects/SubName/MRI/G_ratio_samples.mat (G-ratio values sampled along the transcallosal visual tract: [# MRI_gratio samples x 1])
%       paths.study/Subjects/SubName/MRI/Tract_length.mat (Length of the transcallosal visual tract: double)
%		paths.mat, params.mat
%
%   FILES GENERATED:
%       paths.study/Subjects/SubName/MRI/MicrostructureEstimation/AxonalDistribution_FixM.mat (Results of the model for subject SubName: structure with the estimated and used parameters)
%       paths.study/Subjects/SubName/MRI/MicrostructureEstimation/P_r.png (Plot of P(r) for subject SubName)
%       paths.study/Subjects/SubName/MRI/MicrostructureEstimation/G_r.png (Plot of P(r) for subject SubName)
%	    paths.study/Subjects/SubName/MRI/Velocity.mat (estimated velocity with Delay_IHTT_LeftVF.mat and Tract_length.mat: double)
%       paths.study/Study_Group_Results/G_r_group.png (Plot of G(r) for all subjects)
%
%   EXTERNAL PACKAGES USED:
%       spm, Optimization Toolbox from MATLAB
%  
%   SYNTAX:
%       MODEL()
%
% -------------------------------------------------------------------------
% Author: Rita Oliveira
% Email: anaritaveigaoliveira@gmail.com
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
addpath(paths.spm);
addpath(genpath((paths.scripts)))

sub_list    = params.sub_list;
cond_name   = 'Left';

% Load the IHTT estimation
load(fullfile(paths.group,sprintf('Delay_IHTT_%sVF.mat',cond_name)))
load(fullfile(paths.group,sprintf('ConfidenceInterval_IHTT_%sVF.mat',cond_name)))

% Run model
EstimateMicrostructure_GroupDelay(sub_list,paths.data,paths.group,delay,ci)

cd(fullfile(paths.scripts,'main'))

end
