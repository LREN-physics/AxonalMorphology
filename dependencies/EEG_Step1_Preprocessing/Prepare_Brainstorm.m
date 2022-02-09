function Prepare_Brainstorm(sub_list,cond_list,path_data)
% Reorganizes the EEG data trials to be used in Brainstorm.
%
%   INPUTS: 
%       sub_list (cell array)       - list of subject we want to run
%       cond_list (cell array)      - list of conditions used (e.g. Left and
%           Right visual field stimulations)
%       path_data (string)          - path to the subjects data folder 
%
%   OUPUTS:
%       none
%       
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
    
    % define subject
    sub_name    = sub_list{sub};
    pathopreproc  = fullfile(path_data,sub_name,'EEG','Preprocessed');
    fprintf('Working on subject %s\n', sub_name);
    
    % LOOP FOR CONDITIONS
    for cond=1:length(cond_list)
        cond_name  = cond_list{cond};
        fprintf('Working on condition %s\n', cond_name);
        
        % load data
        dataAverage     = load(fullfile(pathopreproc,[sub_name,'_trials_',cond_name,'.mat']),'dataAverage');
        dataAverage     = dataAverage.dataAverage;
        
        % transform data
        dataAverage2 = dataAverage;
        dataAverage2.trial = [];
        for i=1:length(dataAverage.sampleinfo)
            dataAverage2.trial{i} = squeeze(dataAverage.trial(i,:,:));
        end
        dataAverage2.dimord = 'chan_time';
        
        % save data
        save(fullfile(pathopreproc,['dataTrials_for_BST_',sub_name,'_',cond_name,'.mat']),'dataAverage2')
   
    end % END CONDITIONS
    
end % END SUBJECTS

end % END FUNTION
