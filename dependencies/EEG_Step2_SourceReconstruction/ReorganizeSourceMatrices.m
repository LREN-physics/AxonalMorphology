function ReorganizeSourceMatrices(sub_list,cond_list,source_list,db_path,path_data,path_group_results)
% Reorganizes the matrices obtained with Brainstorm.
% This script gets the data from the source reconstruction using brainstorm
% and transforms it into matrix easier to manipulate
%
%   INPUTS: 
%       sub_list (cell array)       - list of subject we want to run
%       cond_list (cell array)      - list of conditions used (e.g. Left and
%           Right visual field stimulations)
%       source_list (cell array)    - list of regions of interest to look at
%           (e.g. V1V2 - primary and secondary visual cortices)
%       db_path (string)            - path to brainstorm database folder
%       path_data (string)          - path to the subjects data folder 
%       path_group (string)         - path to the group results folder 
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
% Last updated: 20/10/2021
%
%------------- BEGIN CODE --------------

disp('Reorganizing matrices')

% LOOP FOR CONDITIONS
for cond=1:length(cond_list)
    cond_name  = cond_list{cond};
    fprintf('Working on condition %s\n', cond_name);
    close all
    
    
    % LOOP FOR SOURCES
    for source=1:length(source_list)
        scout_name=source_list{source};
        
        
        % LOOP FOR SUBJECTS
        for sub=1:length(sub_list)
            sub_name    = sub_list{sub};
            fprintf('Working on subject %s\n', sub_name);
            pathoutput  = fullfile(path_data,sub_name,'EEG','Delays');
            if~exist(pathoutput,'dir')
                mkdir(pathoutput)
            end
            
            % get vertices
            cd(pathoutput)
            cd ..
            vertices = load(sprintf('%s_vertices.mat',scout_name));
            vertices = vertices.(scout_name);
            V1_L = vertices(1).Vertices;
            V1_R = vertices(2).Vertices;
            Vertices_Interest = [V1_L, V1_R];
            
            % --------- Left hemisphere ---------
            clear aux left_brain_hemisphere
            
            % get source data
            cd(fullfile(db_path,sub_name,cond_name))
            files           = dir('matrix_scout*.mat');
            for i=1:length(files)
                load(files(i).name,'Comment')
                if contains(Comment,[scout_name,'LAll'])
                    aux        = load(files(i).name);
                    aux.Value = aux.Value;
                end
            end
            
            % manipule data
            k=1;
            for i=1:length(V1_L):length(aux.Value)-length(V1_L)+1
                left_brain_hemisphere(k,:,:)= aux.Value(i:i+length(V1_L)-1,:)*10^12;
                k=k+1;
            end
                       
            
            % --------- Right hemisphere ---------
            clear aux right_brain_hemisphere
            
            % get source data
            for i=1:length(files)
                load(files(i).name,'Comment')
                if contains(Comment,[scout_name,'RAll'])
                    aux        = load(files(i).name);
                    aux.Value = aux.Value;
                end
            end
            
            % manipule data
            k=1;
            for i=1:length(V1_R):length(aux.Value)-length(V1_R)+1
                right_brain_hemisphere(k,:,:)= aux.Value(i:i+length(V1_R)-1,:)*10^12;
                k=k+1;
            end
            
            % --------- Define data and variables ---------
            time_vec=aux.Time*10^3;
               
            save(fullfile(pathoutput,sprintf('CD_%sVisualField_LeftBrainV1V2.mat',cond_name')),'left_brain_hemisphere')
            save(fullfile(pathoutput,sprintf('CD_%sVisualField_RightBrainV1V2.mat',cond_name')),'right_brain_hemisphere')
            
            save(fullfile(path_group_results,'time_vec.mat'),'time_vec')
            
        end % end subjects
        
    end % end source
    
    
end % end conditions


end % end function




