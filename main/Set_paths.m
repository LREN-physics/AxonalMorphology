function Set_paths()
% Creates the path structure used in all the other scripts.
%
%   INPUTS: 
%       none
%
%   OUPUTS:
%       none
%
%   FILES USED: 
%       none
%
%   FILES GENERATED: 
%       paths.mat
%
%   SYNTAX:
%       Set_paths()
%
%   HOW TO USE:
%       Manually edit the script and change the path variables to match the
%       ones on your local folder structure. Then run the script.
%			paths.study             (Path where the data for this project is located)
%			paths.scripts           (Path where the scripts are located)
%			paths.eeglab            (Path where the EEGlab software folder is located)
%			paths.fieldtrip         (Path where the Fieldtrip software folder is located)
%			paths.brainstorm_func   (Path where the Brainstorm software database is located, subfolder of the EEG/functional data)
%			paths.brainstorm_anat   (Path where the Brainstorm software database is located, subfolder of the MRI/structural data)
%			paths.brainstorm        (Path where the Brainstorm software folder is located)
%			paths.spm               (Path where the SPM software folder is located)
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

paths = [];


% -------------------- EDIT HERE --------------------
paths.study             = '.../AxonalMorphologyProject/example_data';
paths.scripts           = '.../AxonalMorphologyProject';
paths.eeglab            = 'yourpath/eeglab13_4_4b';
paths.fieldtrip         = 'yourpath/fieldtrip-20191206';
paths.brainstorm_func   = 'yourpath/brainstorm_db/data';
paths.brainstorm_anat   = 'yourpath/brainstorm_db/anat';
paths.brainstorm        = 'yourpath/brainstorm3';
paths.spm               = 'yourpath/spm12';

% -----------------------------------------------------

% Generate remaining paths
paths.data          = fullfile(paths.study,'Subjects');
paths.group         = fullfile(paths.study,'Study_Group_Results');

% Save structure
save(fullfile(paths.scripts,'main','paths.mat'),'paths')

% Create folders
if~exist(paths.data,'dir')
    mkdir(paths.data)
end
if~exist(paths.group,'dir')
    mkdir(paths.group)
end     

end
