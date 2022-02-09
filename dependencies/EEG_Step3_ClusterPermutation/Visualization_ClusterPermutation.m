function Visualization_ClusterPermutation(sub_list,cond_list,scout_list,path_data,path_group_results)
% Plots the source activations in a given region of interest
% after cluster permutation. Plots the number of significant voxels as
% well. The plots are made for each subject and the average fo subjects
%
%   INPUTS: 
%       sub_list (cell array)       - list of subject we want to run
%       cond_list (cell array)      - list of conditions used (e.g. Left and
%           Right visual field stimulations)
%       scout_list (cell array)    - list of regions of interest to look at
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

close all

% LOOP FOR CONDITIONS
for cond=1
    cond_name=cond_list{cond};
    
    % define time of interest to search for peaks
    if strcmp(cond_name,'Left')
        time_of_interest = [140, 160, 100, 150];
    elseif strcmp(cond_name,'Right')
        time_of_interest = [100, 150, 140, 160];
    end
        
    % LOOP FOR SOURCES
    for s=1:length(scout_list)
        scout_name=scout_list{s};
        
        
        %% Figure for individual subjects
        figure('units','normalized','outerposition',[0 0 1 1]);
        tiledlayout(ceil(length(sub_list)/2),2,'TileSpacing','Compact','Padding','Compact');
        
        % LOOP FOR SUBJECTS
        for sub=1:length(sub_list)
            
            % define paths
            sub_name=sub_list{sub};
            path_sub = fullfile(path_data,sub_name,'EEG','Statistics_Source');
            
            % load data
            load(fullfile(path_sub,sprintf('Stats_Source_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'stat');
            load(fullfile(path_sub,sprintf('separate_complete_L_am_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'separate_complete_L_am')
            load(fullfile(path_sub,sprintf('separate_complete_R_am_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'separate_complete_R_am')
            load(fullfile(path_sub,sprintf('separate_L_am_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'separate_L_am')
            load(fullfile(path_sub,sprintf('separate_R_am_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'separate_R_am')
            
            % define time
            if sub==1
                time_vec=stat.time*10^3;
            end
            
            
            % initiate variable
            if sub==1
                Subject_matrix_left = zeros(length(separate_complete_L_am),length(sub_list));
                Subject_matrix_right = zeros(length(separate_complete_L_am),length(sub_list));
                Subject_matrix_left_stat = zeros(length(separate_complete_L_am),length(sub_list));
                Subject_matrix_right_stat = zeros(length(separate_complete_L_am),length(sub_list));
            end
            
            % --- plot subject sources ---
            nexttile
            plot(stat.time*10^3,abs(separate_complete_L_am),'--b','LineWidth',1.5)
            hold on
            plot(time_vec,abs(separate_complete_R_am),'--r','LineWidth',1.5)
            hold on
            plot(time_vec,abs(separate_L_am),'-b','LineWidth',0.7)
            hold on
            plot(time_vec,abs(separate_R_am),'-r','LineWidth',0.7)
            
            % get time of maximum on left side
            [~,beg_sample] = min(abs(time_vec-time_of_interest(1)));
            [~,end_sample] = min(abs(time_vec-time_of_interest(2)));
            [~,idx_max]= max(separate_L_am(beg_sample:end_sample));
            idx_max=idx_max+beg_sample-1;
            time_max_L=time_vec(idx_max);
            
            % get time of maximum on right side
            [~,beg_sample] = min(abs(time_vec-time_of_interest(3)));
            [~,end_sample] = min(abs(time_vec-time_of_interest(4)));
            [~,idx_max]= max(separate_R_am(beg_sample:end_sample));
            idx_max=idx_max+beg_sample-1;
            time_max_R=time_vec(idx_max);
            
            % finish
            xline(time_max_L,':')
            xline(time_max_R,':')
            title(sub_name,'Interpreter','none')
            % save data matrix
            Subject_matrix_left(:,sub)=abs(separate_complete_L_am)';
            Subject_matrix_right(:,sub) = abs(separate_complete_R_am)';
            Subject_matrix_left_stat(:,sub) = abs(separate_L_am)';
            Subject_matrix_right_stat(:,sub) = abs(separate_R_am)';
            
            
        end % end subjects
        clear separate_complete_L_am separate_complete_R_am separate_L_am separate_R_am
        
        ylabel('Activation pAm')
        xlabel('Time (ms)')
        xlabel('Post-stimulus period (ms)')   
        sgtitle(sprintf('%s Visual Field - %s Area - Cluster Permutation',cond_name,scout_name))
        legend({'Left','Right','Left sig','Right sig'},'Location','northwest')
        saveas(gcf,fullfile(path_group_results,sprintf('Source_group_CP_%s_%s_1.png',scout_name,cond_name)))
        
        
        % grand average of subjects
        GrandAvg_left = mean(Subject_matrix_left,2);
        GrandAvg_right = mean(Subject_matrix_right,2);
        GrandAvg_left_stat = mean(Subject_matrix_left_stat,2);
        GrandAvg_right_stat = mean(Subject_matrix_right_stat,2);
        
        %% Figure for mean all subjects
        figure('units','normalized','outerposition',[0 0 1 1]);
        tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        
        plot(time_vec,GrandAvg_left,'--b','LineWidth',1.5)
        hold on
        plot(time_vec,GrandAvg_right,'--r','LineWidth',1.5)
        hold on
        plot(time_vec,GrandAvg_left_stat,'-b','LineWidth',0.7)
        hold on
        plot(time_vec,GrandAvg_right_stat,'-r','LineWidth',0.7)
        ylabel('Activation pAm')
        xlabel('Time (ms)')
        xlabel('Post-stimulus period (ms)')
        
        % get time of maximum on left side
        [~,beg_sample] = min(abs(time_vec-time_of_interest(1)));
        [~,end_sample] = min(abs(time_vec-time_of_interest(2)));
        [~,idx_max]= max(GrandAvg_left_stat(beg_sample:end_sample));
        idx_max=idx_max+beg_sample-1;
        time_max_L=time_vec(idx_max);
        
        % get time of maximum on right side
        [~,beg_sample] = min(abs(time_vec-time_of_interest(3)));
        [~,end_sample] = min(abs(time_vec-time_of_interest(4)));
        [~,idx_max]= max(GrandAvg_right_stat(beg_sample:end_sample));
        idx_max=idx_max+beg_sample-1;
        time_max_R=time_vec(idx_max);
        
        delay = time_max_L - time_max_R;
        
        % add xlines
        xline(time_max_L,':')
        xline(time_max_R,':')
        title(sprintf('Group Average - %s Visual Field - %s Area - Cluster Permutation \n Delay: %0.02f ms',cond_name,scout_name,delay))
        legend({'Left','Right','Left sig','Right sig'},'Location','northwest')
        saveas(gcf,fullfile(path_group_results,sprintf('Source_group_CP_%s_%s_2.png',scout_name,cond_name)))
        
        
        %% Number of significant voxels for each subject
        figure('units','normalized','outerposition',[0 0 1 1]);
        tiledlayout(ceil(length(sub_list)/2),2,'TileSpacing','Compact','Padding','Compact');
        
        N_group_L = [];
        N_group_R = [];
        
        % LOOP FOR SUBJECTS
        for sub=1:length(sub_list)
            sub_name=sub_list{sub};
            path_sub = fullfile(path_data,sub_name,'EEG','Statistics_Source');
            
            % load vertices
            cd(path_sub)
            cd ..
            vertices = load(sprintf('%s_vertices.mat',scout_name));
            vertices = vertices.(scout_name);
            V1_L = vertices(1).Vertices;
            V1_R = vertices(2).Vertices;
            
            % load data
            load(fullfile(path_sub,sprintf('Stats_Source_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'stat');
            load(fullfile(path_sub,sprintf('GA_diff_%s_%s_%s.mat',sub_name,cond_name,scout_name)));
            
            % get positive and negative clusters
            pos_cluster_pvals = [stat.posclusters(:).prob];
            neg_cluster_pvals = [stat.negclusters(:).prob];
            
            pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
            pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
            
            % mask the stats
            neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
            neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
            alltog = pos + neg;
            mask_positive = sign(GA_diff.avg);
            mask_positive(mask_positive~=1)=0;
            result = alltog.*mask_positive;
            
            % get number of vertices active on the left
            pos_V1_L = result(1:length(V1_L),:);
            n_V1_L = [];
            for j=1:size(result,2)
                n_V1_L = [n_V1_L sum(pos_V1_L(:,j))];
            end
            
            % get number of vertices active on the right
            n_V1_L = n_V1_L/length(V1_L);
            pos_V1_R = result(length(V1_L)+1:end,:);
            n_V1_R = [];
            for j=1:size(result,2)
                n_V1_R = [n_V1_R sum(pos_V1_R(:,j))];
            end
            n_V1_R = n_V1_R/length(V1_R);
            
            % --- plot in figure ---
            nexttile
            plot(time_vec,n_V1_L*100,'--b','LineWidth',1.5)
            hold on
            plot(time_vec,n_V1_R*100,'--r','LineWidth',1.5)
            
            title(sub_name,'Interpreter','none')
            xlim([0 220])
            
            N_group_L = [N_group_L; n_V1_L];
            N_group_R = [N_group_R; n_V1_R];
            
        end % end subjects
        
        ylabel('% sig voxels')
        xlabel('Time (ms)')
        legend('Left','Right')
        sgtitle(sprintf('%s Visual Field - %s Area: Percentage significant voxels',cond_name,scout_name))
        
        saveas(gcf,fullfile(path_group_results,sprintf('Source_group_CP_%s_%s_3.png',scout_name,cond_name)))
        
        %% Number of significant voxels for mean subjects
        figure('units','normalized','outerposition',[0 0 1 1]);
        tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
        plot(time_vec,mean(N_group_L,1)*100,'--b','LineWidth',1.5)
        hold on
        plot(time_vec,mean(N_group_R,1)*100,'--r','LineWidth',1.5)
        legend('Left','Right')
        xlim([0 220])
        ylabel('% sig voxels')
        xlabel('Time (ms)')
        
        sgtitle(sprintf('Group Average - %s Visual Field - %s Area: Percentage significant voxels mean group',cond_name,scout_name))
        saveas(gcf,fullfile(path_group_results,sprintf('Source_group_CP_%s_%s_4.png',scout_name,cond_name)))
        
    end % end sources

end % end conditions

end % end functions
