function Source_Brainstorm_Trials_kernel(sub_list,cond_list,step,db_path)
% Performs source analysis in Brainstorm.
% These are scripts to be used with Brainstorm. 
% The GUI has to be open for these functions to work.
%
%   INPUTS: 
%       sub_list (cell array)       - list of subject we want to run
%       cond_list (cell array)      - list of conditions used (e.g. Left and
%           Right visual field stimulations)
%       step (string)               - step to run (see below)
%       db_path (string)            - path to brainstorm database folder
%
%   OUPUTS:
%       none
%       
%   DESCRIPTION OF EACH STEP:

    % matlab: when importing a file from fieldtrip the units of the
    % voltages are in microVolts, while brainstorm works with SI Units, so
    % I use a code to transform the units. The original files are
    % overwriten
    
    % source: creates a kernel for the source analysis (computational
    % easier than computing a source grid for each trial). Here it also
    % creates the noise matrix  
    
    % scout all V1V2: get time series of all voxels for all trials in V1V2_Left
    % and V1V2_Right areas together and separately (areas have to be defined
    % a priori for each subject in the GUI)
        
    % average: average all the trials, to be easier to see the mean in the
    % ROIs of interest on the average signal
       
% Author: Rita Oliveira
% Email: rita.oliveira.uni@gmail.com
% Laboratory for Research in Neuroimaging (LREN)
% Department of Clinical Neuroscience, Lausanne  University Hospital and University of Lausanne
% Mont-Paisible 16, CH-1011 Lausanne, Switzerland
%
% Last updated: 29/07/2021
%
%------------- BEGIN CODE --------------


% LOOP FOR SUBJECTS
for sub=1:length(sub_list)
    sub_name    = sub_list{sub};
    fprintf('Working on subject %s\n', sub_name);
    
    % LOOP FOR CONDITIONS
    for cond=1:length(cond_list)
        cond_name  = cond_list{cond};
        fprintf('Working on condition %s\n', cond_name);
        
        % RUN STEP
        switch step
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'matlab'
                
                % -------------- Clean older files --------------
                cd(fullfile(db_path,sub_name,cond_name))
                files = dir('*_matlab.mat');
                sFiles = {};
                if ~isempty(files)
                    for i=1:length(files)
                        sFiles{end+1} = fullfile(sub_name,cond_name,files(i).name);
                    end
                end
                sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], ...
                    'target', 1);
                
                % ------------------ Get files -----------------
                cd(fullfile(db_path,sub_name,cond_name))
                files = dir('data_dataTrials*.mat'); 
                sFiles ={};
                for i=1:length(files)
                    % Input files
                    file=files(i);
                    sFiles{end+1} = [sub_name,'/',cond_name,'/',file.name];
                end
                
                % --------------------- Run -------------------
                sFiles = bst_process('CallProcess', 'process_matlab_eval', sFiles, [], ...
                    'matlab',      ['% Available variables: Data, TimeVector' 10 '' 10 'Data = Data*1e-6;' 10 ''], ...
                    'sensortypes', 'EEG', ...
                    'overwrite',   1);
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'source'
                
                % -------------- Clean older files --------------
                cd(fullfile(db_path,sub_name,cond_name))
                files = dir('*_KERNEL_*.mat');
                if ~isempty(files)
                    for i=1:length(files)
                        sFiles = {fullfile(sub_name,cond_name,files(i).name)};
                        sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], ...
                            'target', 1);
                    end
                end
                 
                % ------------------ Get files -----------------
                cd(fullfile(db_path,sub_name,cond_name))
                files = dir('data_dataTrials*_matlab.mat');
                sFiles ={};
                for i=1:length(files)
                    % Input files
                    file=files(i);
                    sFiles{end+1} = [sub_name,'/',cond_name,'/',file.name];
                end
                
                % --------------------- Run -------------------
                % Process: Compute covariance (noise or data)
                sFiles = bst_process('CallProcess', 'process_noisecov', sFiles, [], ...
                    'baseline',       [ ], ...
                    'datatimewindow', [ ], ...
                    'sensortypes',    'EEG', ...
                    'target',         1, ...  % Noise covariance     (covariance over baseline time window)
                    'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
                    'identity',       0, ...
                    'copycond',       0, ...
                    'copysubj',       0, ...
                    'copymatch',      0, ...
                    'replacefile',    1);  % Replace
                 
                % ------------------ Get files -----------------
                cd(fullfile(db_path,sub_name,cond_name))
                files = dir('data_dataTrials*_matlab.mat');     
                sFiles ={};
                for i=1:length(files)
                    % Input files
                    file=files(i);
                    sFiles{end+1} = [sub_name,'/',cond_name,'/',file.name];
                end
                
                % --------------------- Run -------------------
                % Process: Compute sources [2018]
                sFiles = bst_process('CallProcess', 'process_inverse_2018', sFiles, [], ...
                    'output',  1, ...  % Kernel only: shared
                    'inverse', struct(...
                    'Comment',        'MN: EEG', ...
                    'InverseMethod',  'minnorm', ...
                    'InverseMeasure', 'amplitude', ...
                    'SourceOrient',   {{'fixed'}}, ...
                    'Loose',          0.2, ...
                    'UseDepth',       1, ...
                    'WeightExp',      0.5, ...
                    'WeightLimit',    10, ...
                    'NoiseMethod',    'reg', ...
                    'NoiseReg',       0.1, ...
                    'SnrMethod',      'fixed', ...
                    'SnrRms',         1e-06, ...
                    'SnrFixed',       3, ...
                    'ComputeKernel',  1, ...
                    'DataTypes',      {{'EEG'}}));
                               
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            case 'scout all V1V2'
                
                %% Scout Hemishpheres together
                % -------------- Clean older files --------------
                cd(fullfile(db_path,sub_name,cond_name))
                files = dir('matrix_scout*.mat');
                if ~isempty(files)
                    for i=1:length(files)
                        load(files(i).name,'Comment');
                        if contains(Comment,'TimeSeriesV1V2All')
                            sFiles = {fullfile(sub_name,cond_name,files(i).name)};
                            sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], ...
                                'target', 1);
                        end
                    end
                end
                
                % ---------------- Get kernel ----------------
                files = dir('data_*_matlab.mat');
                cd(fullfile(db_path,sub_name,cond_name))
                kernel = dir('*_KERNEL_*.mat');
                
                % ------------------ Get files -----------------
                sFiles ={};
                for i=1:length(files)
                    % Input files
                    file=files(i);
                    sFiles{end+1} = ['link|',sub_name,'/',cond_name,'/',kernel.name,'|',sub_name,'/',cond_name,'/',file.name];
                end
                
                % --------------------- Run -------------------
                % Process: get time series of all voxels in V1V2_Left and 
                % V1V2_Right areas together
                sFiles = bst_process('CallProcess', 'process_extract_scout', sFiles, [], ...
                    'timewindow',     [ ], ...
                    'scouts',         {'V1V2_together', {'V1_exvivo L & V1_exvivo R & V2_exvivo L & V2_exvivo R'}}, ...
                    'scoutfunc',      5, ...  % All
                    'isflip',         1, ...
                    'isnorm',         0, ...
                    'concatenate',    1, ...
                    'save',           1, ...
                    'addrowcomment',  1, ...
                    'addfilecomment', 1);               
                sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                    'tag',           'TimeSeriesV1V2All', ...
                    'isindex',       0);
                
                %% Scout Hemishpheres separately - Left
                % -------------- Clean older files --------------
                cd(fullfile(db_path,sub_name,cond_name))
                files = dir('matrix_scout*.mat');
                if ~isempty(files)
                    for i=1:length(files)
                        load(files(i).name,'Comment');
                        if contains(Comment,'TimeSeriesV1V2LAll') || contains(Comment,'TimeSeriesV1V2RAll')
                            sFiles = {fullfile(sub_name,cond_name,files(i).name)};
                            sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], ...
                                'target', 1);
                        end
                    end
                end
                
                % ---------------- Get kernel ----------------
                cd(fullfile(db_path,sub_name,cond_name))
                files = dir('data_*_matlab.mat');
                cd(fullfile(db_path,sub_name,cond_name))
                kernel = dir('*_KERNEL_*.mat');
                
                % ------------------ Get files -----------------
                sFiles ={};
                for i=1:length(files)
                    % Input files
                    file=files(i);
                    sFiles{end+1} = ['link|',sub_name,'/',cond_name,'/',kernel.name,'|',sub_name,'/',cond_name,'/',file.name];
                end
                
                % --------------------- Run -------------------
                % Process: get time series of all voxels in V1V2_Left 
                sFiles = bst_process('CallProcess', 'process_extract_scout', sFiles, [], ...
                    'timewindow',     [ ], ...
                    'scouts',         {'V1V2_separate', {'V1_exvivo L & V2_exvivo L'}}, ...
                    'scoutfunc',      5, ...  % All
                    'isflip',         1, ...
                    'isnorm',         0, ...
                    'concatenate',    1, ...
                    'save',           1, ...
                    'addrowcomment',  1, ...
                    'addfilecomment', 1);  
                sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                    'tag',           'TimeSeriesV1V2LAll', ...
                    'isindex',       0);
                
                %% Scout Hemishpheres separately - Right
                % ------------------ Get files -----------------
                sFiles ={};
                for i=1:length(files)
                    % Input files
                    file=files(i);
                    sFiles{end+1} = ['link|',sub_name,'/',cond_name,'/',kernel.name,'|',sub_name,'/',cond_name,'/',file.name];
                end
                
                % --------------------- Run -------------------
                % Process: get time series of all voxels in V1_Right
                sFiles = bst_process('CallProcess', 'process_extract_scout', sFiles, [], ...
                    'timewindow',     [ ], ...
                    'scouts',         {'V1V2_separate', {'V1_exvivo R & V2_exvivo R'}}, ...
                    'scoutfunc',      5, ...  % All
                    'isflip',         1, ...
                    'isnorm',         0, ...
                    'concatenate',    1, ...
                    'save',           1, ...
                    'addrowcomment',  1, ...
                    'addfilecomment', 1);   
                sFiles = bst_process('CallProcess', 'process_set_comment', sFiles, [], ...
                    'tag',           'TimeSeriesV1V2RAll', ...
                    'isindex',       0);
                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
            case 'average'
                % -------------- Clean older files --------------
                cd(fullfile(db_path,sub_name,cond_name))
                files = dir('*average*.mat');
                if ~isempty(files)
                    for i=1:length(files)
                        sFiles = {fullfile(sub_name,cond_name,files(i).name)};
                        sFiles = bst_process('CallProcess', 'process_delete', sFiles, [], ...
                            'target', 1);
                    end
                end
                
                % ------------------ Get files -----------------
                files = dir('data_dataTrials*_matlab.mat');
                sFiles ={};
                for i=1:length(files)
                    % Input files
                    file=files(i);
                    sFiles{end+1} = [sub_name,'/',cond_name,'/',file.name];
                end
                
                % --------------------- Run -------------------
                % Process: Average: By trial group (folder average)
                sFiles = bst_process('CallProcess', 'process_average', sFiles, [], ...
                    'avgtype',       5, ...  % By trial group (folder average)
                    'avg_func',      1, ...  % Arithmetic average:  mean(x)
                    'weighted',      0, ...
                    'keepevents',    0);
                
        end % END STEP
    end % END CONDITIONS
end % END SUBJECTS
end % END FUNCTION
