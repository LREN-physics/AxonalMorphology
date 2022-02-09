function Clust_Perm_Source(sub_list,cond_list,step,cfg_input,path_data,path_group_results,db_path,db_path_anat,source_list)
% Applies cluster permutation to the time-course of the source activations
% in a region of interest.
% This script uses data from a patch containing both regions
% of interest. The data was taken as flipping the values to just one
% orientation that represents that patch. 
%
%   INPUTS: 
%       sub_list (cell array)       - list of subject we want to run
%       cond_list (cell array)      - list of conditions used (e.g. Left and
%           Right visual field stimulations)
%       step (string)               - step to run (see below)
%       cfg_input (structure)       - structure with the information
%       regarding the cluster permutation parameters
%       path_data (string)          - path to the subjects data folder 
%       db_path (string)            - path to brainstorm database functional folder
%       db_path_anat (string)       - path to brainstorm database anatamoical folder
%       source_list (cell array)    - list of regions of interest to look at
%           (e.g. V1V2 - primary and secondary visual cortices)
%
%   OUPUTS:
%       none
% 
%   DESCRIPTION OF EACH STEP (and files created):
    
    % prepare: load data and transform it to prepare for cluster
    % permutation. Also prepares the neighbours of the vertices and the
    % average of the trials to plot later
    
    % run: runs the cluster permutation analysis
    
    % get metrics: reorganize the matrices of the results for later use       
       
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
    
    % define subject and paths
    sub_name    = sub_list{sub};
    pathoutput  = fullfile(path_data,sub_name,'EEG','Statistics_Source');
    if~exist(pathoutput,'dir')
        mkdir(pathoutput)
    end
    fprintf('Working on subject %s\n', sub_name);
    
    % LOOP FOR CONDITIONS
    for cond=1:length(cond_list)
        cond_name  = cond_list{cond};
        fprintf('Working on condition %s\n', cond_name);
        close all
        
        % STEP TO RUN (only runs what is specified in the string step)
        switch step
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'prepare'
                
                % Load data
                dataAverage     = load(fullfile(path_data,sub_name,'EEG','Preprocessed',sprintf('%s_trials_%s.mat',sub_name,cond_name)),'dataAverage');
                dataAverage     = dataAverage.dataAverage;
                
                % LOOP FOR SOURCES
                for source=1:length(source_list)
                    
                    %% Load and manipulate data
                    clear alltrials Vertices_Interest V1_L V1_R Vertices
                    scout_name = source_list{source};
                    
                    % Load vertices
                    close all
                    cd(pathoutput)
                    cd ..
                    vertices = load(sprintf('%s_vertices.mat',scout_name));
                    vertices = vertices.(scout_name);
                    V1_L = vertices(1).Vertices;
                    V1_R = vertices(2).Vertices;
                    Vertices_Interest = [V1_L, V1_R];
                    clear vertices V1_L V1_R
                    
                    % Load cortex vertices locations
                    load('cortex.mat');
                    
                    % Load values of source
                    cd(fullfile(db_path,sub_name,cond_name))
                    files           = dir('matrix_scout*.mat');
                    for i=1:length(files)
                        load(files(i).name,'Comment')
                        if contains(Comment,[scout_name,'All']) 
                            aux        = load(files(i).name);
                            aux.Value = aux.Value*10^12;
                        end
                    end
                    
                    % Check
                    if size(dataAverage.trial,1)*length(Vertices_Interest)~=length(aux.Description)
                        disp('PROBLEM')
                    end
               
                    % Rearange matrix of data
                    dataAverage.trial = [];
                    k=1;
                    for i=1:length(Vertices_Interest):length(aux.Value)-length(Vertices_Interest)+1
                        dataAverage.trial(k,:,:)= aux.Value(i:i+length(Vertices_Interest)-1,:);
                        k=k+1;
                    end
                            
                    % Period to analyse
                    begsample   = cfg_input.begwindow; 
                    endsample   = cfg_input.endwindow;
                    period = endsample-begsample+1;
                    
                    % -- Define baseline --
                    % redefine structure
                    data_baseline.trial     = [];
                    data_baseline.time      = [];
                    for kk = 1:size(dataAverage.trial,1)
                        data_baseline.trial(kk,:,:) = dataAverage.trial(kk,:,1:cfg_input.baseline);
                    end
                    % do mean along time of baseline
                    for trial = 1:size(data_baseline.trial,1) % baseline
                        for ch = 1:size(data_baseline.trial,2) % baseline
                            data_baseline.trial(trial,ch,:)=mean(data_baseline.trial(trial,ch,:))*ones(1,cfg_input.baseline);
                        end
                    end
                    % extend baseline to the same period as activation
                    % (that we want to analyse)
                    data_baseline_2.trial = zeros(size(data_baseline.trial,1),size(data_baseline.trial,2),period);
                    for trial = 1:size(data_baseline.trial,1) % baseline
                        for ch = 1:size(data_baseline.trial,2) % baseline
                            data_baseline_2.trial(trial,ch,:)=data_baseline.trial(trial,ch,1)*ones(1,period);
                        end
                    end
                    data_baseline.trial = data_baseline_2.trial;
                    
                    % -- Define post-stimulus period --
                    % redefine structure
                    data_activation.trial   = [];
                    data_activation.time    = [];
                    for kk = 1:size(dataAverage.trial,1)
                        data_activation.trial(kk,:,:) = dataAverage.trial(kk,:,begsample:endsample);
                    end
                    
                    % Set time to the same as activation
                    data_activation.time(:,:)   = dataAverage.time(begsample:endsample);
                    data_baseline.time(:,:)     = data_activation.time;
                                      
                    % Define variables
                    time_vec    = data_baseline.time;
                    ncond   = 2;
                    ntrial  = size(data_baseline.trial,1);
                    nchan   = size(data_baseline.trial,2);
                    ntf     = size(data_baseline.trial,3);
                    
                    % Electrodes
                    brain               = load(fullfile(db_path_anat,sub_name,'tess_cortex_pial_low.mat'));
                    channames           = Vertices_Interest';
                    elec                = [];
                    elec.chanpos        = brain.Vertices(Vertices_Interest,:);
                    for i=1:length(channames)
                        elec.label{i} = ['',num2str(channames(i))];
                        elec.unit{i} = 'm';
                    end
                    elec.label = reshape(elec.label ,[length(channames) 1]);
                    elec.unit = elec.unit';
                                        
                    clear file begsample endsample i kk A B ...
                        data_baseline_2 period files aux aux_L aux_R ...
                        dataAverage aux_0
                    
                    %% Make cell array with all data in separate ft structures
                    
                    % Make structures
                    template_data              = [];
                    template_data.fsample      = cfg_input.sr; %chan x time
                    template_data.avg          = nan(nchan,ntf);% 1 x tf
                    template_data.time         = time_vec; %1 x time
                    template_data.label        = elec.label ;%chan x 1;
                    template_data.dimord       = 'chan_time';
                    template_data.elec         = elec;
                    template_data.cfg          = [];
                    
                    % Make cell array with all data in separate ft structures
                    clear alltrials
                    alltrials       = cell(ntrial,ncond);
                    alltrials(:,:)  = {template_data};
                    clear data_avg chanlocs;
                    
                    for trial = 1:ntrial % baseline
                        alltrials{trial,1}.avg = squeeze(data_baseline.trial(trial,:,:));
                    end
                    for trial = 1:ntrial % activation
                        alltrials{trial,2}.avg = squeeze(data_activation.trial(trial,:,:));
                    end
                    
                    save(fullfile(pathoutput,sprintf('Trials_for_stats_Source_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'alltrials','-v7.3')
                    
                    clear data_activation data_baseline
                    
                    %% Create average 
                    % Average over trials of baseline and post-stimulus period
                    cfg                 = [];
                    cfg.channel         = 'all';
                    cfg.latency         = 'all';
                    cfg.parameter       = 'avg';
                    GA_baseline         = ft_timelockgrandaverage(cfg,alltrials{:,1});
                    GA_activation       = ft_timelockgrandaverage(cfg,alltrials{:,2});
                    GA_activation.avg = abs(GA_activation.avg);
                    GA_baseline.avg = abs(GA_baseline.avg);
                    save(fullfile(pathoutput,['GA_activation_',sub_name,'_',cond_name,'_',scout_name,'.mat']),'GA_activation')
                    save(fullfile(pathoutput,['GA_baseline_',sub_name,'_',cond_name,'_',scout_name,'.mat']),'GA_baseline')
                    % Make difference on the average
                    cfg             = [];
                    cfg.operation   = 'subtract';
                    cfg.parameter   = 'avg';
                    GA_diff         = ft_math(cfg,GA_activation,GA_baseline);
                    save(fullfile(pathoutput,['GA_diff_',sub_name,'_',cond_name,'_',scout_name,'.mat']),'GA_diff')                      
                    clearvars GA_diff GA_baseline GA_activation stat
                    
                
                    %% Make neighbours structure
                    for i=1:length(Vertices_Interest)
                        neighbours(i).label = num2str(Vertices_Interest(i));
                        neighbours(i).neighblabel = [];
                        aux = find(cortex.VertConn(Vertices_Interest(i),:)==1)';
                        for j=1:length(aux)
                            neighbours(i).neighblabel{end+1,1} = num2str(aux(j));
                        end
                    end
                    save(fullfile(pathoutput,sprintf('Neighbours_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'neighbours')
                    
                    
                end % END SOURCE
                
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
           case 'run'
               
                 % LOOP FOR SOURCES
                 for source=1:length(source_list)
                    
                    clear alltrials Vertices_Interest V1_L V1_R Vertices
                    scout_name = source_list{source};           
                    
                    % load data
                    load(fullfile(pathoutput,sprintf('Trials_for_stats_Source_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'alltrials')
                    ntrial  = size(alltrials,1);
                    load(fullfile(pathoutput,sprintf('Neighbours_%s_%s_%s.mat',sub_name,cond_name,scout_name)))
                 
                    %% Cluster permutation
                    % Define variables
                    cfg                     = [];
                    cfg.channel             = 'all';
                    cfg.latency             = 'all';
                    cfg.method              = 'montecarlo';
                    cfg.statistic           = 'depsamplesT';
                    cfg.correctm            = 'cluster';
                    cfg.alpha               = 0.025;
                    cfg.clusteralpha        = 0.001;
                    cfg.avgoverchan         = 'no';
                    cfg.avgovertime         = 'no';
                    cfg.correcttail         = 'prob';
                    cfg.tail                =   0;
                    cfg.numrandomization    = 1000;
                    cfg.uvar                = 1;
                    cfg.ivar                = 2;
                    cfg.neighbours          = neighbours;  % same as defined for the between-trials experiment
                    clearvars neighbours
                    
                    % Define design
                    design = zeros(2,2*ntrial);
                    for i = 1:ntrial
                        design(1,i) = i;
                    end
                    for i = 1:ntrial
                        design(1,ntrial+i) = i;
                    end
                    design(2,1:ntrial)          = 1;
                    design(2,ntrial+1:2*ntrial) = 2;
                    
                    cfg.design = design;
                    
                    % Run analysis
                    [stat] = ft_timelockstatistics(cfg, alltrials{:,2}, alltrials{:,1}); % Activation - Baseline (therefore we are interested in the
                    % positive clusters)
                    
                    save(fullfile(pathoutput,sprintf('Stats_Source_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'stat')
                    clearvars stat
                    
                   
                    
                 end % END SOURCES
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'get metrics'
                
                % LOOP FOR SOURCE
                for source=1:length(source_list)
                    
                    clearvars L R separate_com_L separate_com_R separate_L separate_R aux_L aux_R aux
                    scout_name = source_list{source};
                    
                    % Load vertices
                    close all
                    cd(pathoutput)
                    cd ..
                    vertices = load(sprintf('%s_vertices.mat',scout_name));
                    vertices = vertices.(scout_name);
                    V1_L = vertices(1).Vertices;
                    V1_R = vertices(2).Vertices;
                    Vertices_Interest = [V1_L, V1_R];
                    
                    % Load data
                    stat        = load(fullfile(pathoutput,sprintf('Stats_Source_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'stat');
                    stat        = stat.stat;                    
                    GA_diff = load(fullfile(pathoutput,['GA_diff_',sub_name,'_',cond_name,'_',scout_name,'.mat']),'GA_diff');
                    GA_diff = GA_diff.GA_diff;
                    time_vec = stat.time*10^3;
                    
                    % Get significant clusters
                    pos_cluster_pvals = [stat.posclusters(:).prob];
                    neg_cluster_pvals = [stat.negclusters(:).prob];
                    
                    pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
                    pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
                    
                    neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
                    neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
                    
                    % Consider all the clusters and mask the ones that were
                    % absolute positive
                    alltog = pos + neg;
                    mask_positive = sign(GA_diff.avg);
                    mask_positive(mask_positive~=1)=0;
                    result = alltog.*mask_positive;

                    % Get number of vertices active
                    pos_V1_L = result(1:length(V1_L),:);
                    n_V1_L = [];
                    for i=1:size(result,2)
                        n_V1_L = [n_V1_L sum(pos_V1_L(:,i))];
                    end

                    n_V1_L = n_V1_L/length(V1_L);
                    pos_V1_R = result(length(V1_L)+1:end,:);
                    n_V1_R = [];
                    for i=1:size(result,2)
                        n_V1_R = [n_V1_R sum(pos_V1_R(:,i))];
                    end
                    n_V1_R = n_V1_R/length(V1_R);
                     
                    path_delays  = fullfile(path_data,sub_name,'EEG','Delays');
                    load(fullfile(path_delays,sprintf('CD_%sVisualField_LeftBrainV1V2.mat',cond_name')))
                    load(fullfile(path_delays,sprintf('CD_%sVisualField_RightBrainV1V2.mat',cond_name')))
                    load(fullfile(path_group_results,'time_vec.mat'),'time_vec')

                    L = left_brain_hemisphere;
                    R = right_brain_hemisphere;
                    
                    % Get metrics
                    separate_complete_L = squeeze(mean(L(:,:,find(time_vec*10^-3==stat.time(1)):find(time_vec*10^-3==stat.time(end))),1));
                    separate_complete_R = squeeze(mean(R(:,:,find(time_vec*10^-3==stat.time(1)):find(time_vec*10^-3==stat.time(end))),1));
                    clear aux_R
                    separate_L = separate_complete_L.*result(1:length(V1_L),:);
                    separate_R = separate_complete_R.*result(length(V1_L)+1:end,:);
                    
                    % Doing the mean considering the zero and non-zero
                    % voxels
                    separate_complete_L_am = abs(mean(separate_complete_L,1));
                    separate_complete_R_am = abs(mean(separate_complete_R,1));
                    
                    separate_L_am = abs(mean(separate_L,1));
                    separate_R_am = abs(mean(separate_R,1));

                    save(fullfile(pathoutput,sprintf('separate_complete_L_am_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'separate_complete_L_am')
                    save(fullfile(pathoutput,sprintf('separate_complete_R_am_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'separate_complete_R_am')
                    save(fullfile(pathoutput,sprintf('separate_L_am_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'separate_L_am')
                    save(fullfile(pathoutput,sprintf('separate_R_am_%s_%s_%s.mat',sub_name,cond_name,scout_name)),'separate_R_am')
            
                end % END SOURCES
              
        end % END STEPS
        
    end % END CONDITIONS
    
end % END SUBJECTS

end % END FUNCTION
