function Get_confidence_interval_delay(sub_list,cond_list,source_list,path_data,path_group_results)
% Calculates a confidence interval on the interhemispheric transfer time 
% (IHTT, in miliseconds).
%
%   INPUTS: 
%       sub_list (cell array)       - list of subject we want to run
%       cond_list (cell array)      - list of conditions used (e.g. Left and
%           Right visual field stimulations)
%       source_list (cell array)    - list of regions of interest to look at
%           (e.g. V1V2 - primary and secondary visual cortices)
%       path_data (string)          - path to the subjects data folder 
%       path_group_results (string) - path to group results folder
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


% LOOP FOR CONDITIONS
for cond=1:length(cond_list)
    cond_name  = cond_list{cond};
    fprintf('Working on condition %s\n', cond_name);
    close all
    
    % LOOP FOR SOURCES
    for source=1:length(source_list)
        scout_name=source_list{source};
        
        %% Get data into several partitions
        Subject_matrix_left =[];
        Subject_matrix_right = [];
        Delay = [];
        nsh = 4; % number of partitions
        
        clear Subject_matrix_left Subject_matrix_right
        % LOOP FOR SUBJECTS
        for sub=1:length(sub_list)
            sub_name    = sub_list{sub};
            fprintf('Working on subject %s\n', sub_name);
            pathoutput  = fullfile(path_data,sub_name,'EEG','Delays');
            
            % Load data
            load(fullfile(pathoutput,sprintf('CD_%sVisualField_LeftBrainV1V2.mat',cond_name')))
            load(fullfile(pathoutput,sprintf('CD_%sVisualField_RightBrainV1V2.mat',cond_name')))
            load(fullfile(path_group_results,'time_vec.mat'),'time_vec')
            
            left_brain_hemisphere  = squeeze(mean(left_brain_hemisphere,2));
            right_brain_hemisphere = squeeze(mean(right_brain_hemisphere,2));
            
            partition=floor(size(left_brain_hemisphere,1)/nsh);
            P=randperm(size(left_brain_hemisphere,1));
            left_brain_hemisphere=left_brain_hemisphere(P,:);
            right_brain_hemisphere=right_brain_hemisphere(P,:);
                
            for pp=1:nsh
                left_pp=left_brain_hemisphere((pp-1)*partition+1:pp*partition,:);
                right_pp=right_brain_hemisphere((pp-1)*partition+1:pp*partition,:);
                
                Subject_matrix_left(sub,:,pp)=abs(mean(left_pp));
                Subject_matrix_right(sub,:,pp)=abs(mean(right_pp));
            end
        end
        
        %% Plot results and find confidence interval
        figure('units','normalized','outerposition',[0 0 1 1]);
        t = tiledlayout(ceil(nsh/2),2,'TileSpacing','Compact','Padding','Compact');
        
        for pp=1:nsh
            GrandAvg_left = mean(Subject_matrix_left(:,:,pp),1);
            GrandAvg_right = mean(Subject_matrix_right(:,:,pp),1);
            
            % define partitions trials
            nexttile
            plot(time_vec,GrandAvg_left,'b','LineWidth',1.5)
            hold on
            plot(time_vec,GrandAvg_right,'r','LineWidth',1.5)
            xlabel('Time (ms)')
            ylabel('Amplitude (pAm)')
            
            % define time of interest to search for peaks
            if strcmp(cond_name,'Left')
                time_of_interest = [145, 160, 100, 145];
            elseif strcmp(cond_name,'Right')
                time_of_interest = [100, 145, 145, 160];
            end
            
            % get time of maximum on left side
            [~,beg_sample] = min(abs(time_vec-time_of_interest(1)));
            [~,end_sample] = min(abs(time_vec-time_of_interest(2)));
            [~,idx_max]= max(GrandAvg_left(beg_sample:end_sample));
            idx_max=idx_max+beg_sample-1;
            time_max_L=time_vec(idx_max);
            
            % get time of maximum on right side
            [~,beg_sample] = min(abs(time_vec-time_of_interest(3)));
            [~,end_sample] = min(abs(time_vec-time_of_interest(4)));
            [~,idx_max]= max(GrandAvg_right(beg_sample:end_sample));
            idx_max=idx_max+beg_sample-1;
            time_max_R=time_vec(idx_max);
            
            % get delay
            delay = time_max_L-time_max_R;
            Delay = [Delay, delay];
            
            % extra 
            xline(time_max_L,':')
            xline(time_max_R,':')
            title(sprintf('%s Visual Field - %s Area \n Sub group %d - delay: %0.02f ms.',cond_name,scout_name,pp,delay));
        
        end
        
        % save
        legend({'Left','Right'},'Location','northwest')
        sgtitle(sprintf('Mean delay: %0.02f; std: %0.02f',mean(Delay),std(Delay)))
        saveas(gcf,fullfile(path_group_results,sprintf('Source_group_%s_%s_3.png',scout_name,cond_name)))
        ci = std(Delay);
        save(fullfile(path_group_results,sprintf('ConfidenceInterval_IHTT_%sVF.mat',cond_name)),'ci')
        
    end % end source
    
    
end % end conditions


end % end function
