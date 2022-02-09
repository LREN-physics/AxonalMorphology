function [data] = Artefact_Removal(cfg,data)
%function [data] = my_artefactremoval(cfg,data)
%
%   MY_ARTEFACTREMOVAL allows removal of EEG artefacts specific for blinks, 
%   muscle movement, bad channels, ICA, and visual inspection by
%   user-interaction. 
%   Inputs is a data structure with segmented EEG and EOG data
%
%   example configuration file here
%   cfg                = [];
%   cfg.eeg_chans      = [1:63];
%   cfg.eog_chans      = [64 65];
%   cfg.layout         = lay;
%   cfg.tasks           = [1 2 3 4 5]
%                           1=blinks
%                           2=muscle
%                           3=bad channels
%                           4=ica
%                           5=visual
%
% update 26 04 2017, do not remove the channels of no interest
                          
% update 25 04 2018, following an update of fieldtrip the function for
% cfgx.tasks==4 does not work anymore. The function ft_rejectvisual does
% not keep the channels in the data structure even with the option keep
% channels = yes, but include only the ones that are specified in cfg.channel

% update 2022 

    disp('RUNNING MY_ARTEFACTREMOVAL FUNCTION')

    %make general configuration structure 
    cfgx = cfg;
    cfg  = [];
    
    lay = cfgx.layout;
    
    %% ---------------- Task 1 (EOG) -----------------------------
    if find(cfgx.tasks==1)
        %find blinks
        disp('SELECT BLINK TRIALS')
%         cfg                     = [];
%         cfg.method              = 'summary';%'channel';
%         cfg.keepchannel         = 'yes';
%         cfg.metric              = 'absmax';%'var';
%         cfg.layout              = lay;       % this allows for plotting
%         cfg.channel             = cfgx.eog;   %EOG channels
%         cfg.preproc.bpfilter    = 'yes';
%         cfg.preproc.bpfilttype  = 'but';
%         cfg.preproc.bpfreq      = [1 15];
%         cfg.preproc.bpfiltord   = 4;
%         cfg.preproc.rectify     = 'yes';
%         dummy                   = ft_rejectvisual(cfg,data);
%         %remove blinks
%         disp('REMOVE BLINK TRIALS')
%         cfg                     = [];
%         cfg                     = dummy.cfg;
%         data                   = ft_rejectartifact(cfg,data);
%         clear dummy;
        dataforEOG = data;
        for trial=1:length(data.trial)
            dataforEOG.trial{trial}(129,:)=data.trial{trial}(130,:)-data.trial{trial}(129,:);
             
        end
         
        dataforEOG    = rmfield(dataforEOG,'sampleinfo');        
        cfg = [];
        cfg.channel             = {'FP2','CZ','EOG1'}; 
        cfg.continuous    = 'no';
        cfg.viewmode      = 'vertical';
        cfg.ylim          = [-10 10];
        cfg.verticalpadding = 0.1;
        ft_databrowser(cfg,dataforEOG) 
        
        % chose components to reject
        disp('REMOVE TRIALS WITH EYE MOVEMENTS')
        Prompt(1,:) = {'Enter all trials to be rejected (space-separated - e.g. 1 4 7 12)', 'comp',[]};
        name='Trial Rejection';
        formats(1,1).type = 'edit';
        formats(1,1).format = 'text';
        formats(1,1).limits = [0 1];
        DefAns = struct([]);
        DefAns(1).comp = 'None';
        [answer,canceled] = inputsdlg(Prompt,name,formats,DefAns);
        fprintf('removing trials: %s \n', answer.comp)
        trials_rej = unique(str2num(answer.comp));
        trials_not_rej = 1:length(data.trial);
        trials_not_rej(trials_rej) = [];
        
        % remove trials
        cfg = [];
        cfg.trials = trials_not_rej;
        data = ft_preprocessing(cfg, data);
        close all
        
        % save both data and trials rejected
        data = {data, trials_rej};
    end   

    %% ---------------- Task 2 (never used) -----------------------------
    if find(cfgx.tasks==2)
        %find muscle artifacts
        disp('SELECT MUSCLE ARTEFACT TRIALS')
        cfg                     = [];
        cfg.method              = 'summary';%'channel';
        cfg.keepchannel         = 'yes';
        cfg.metric              = 'absmax';%'var';
        cfg.layout              = lay;       % this allows for plotting
        cfg.channel             = cfgx.eeg;    %EEG channels
        cfg.preproc.bpfilter    = 'yes';
        cfg.preproc.bpfreq      = [110 140];
        cfg.preproc.bpfiltord   =  8;
        cfg.preproc.bpfilttype  = 'but';
        cfg.preproc.rectify     = 'yes';
        cfg.preproc.boxcar      = 0.2;
        dummy                   = ft_rejectvisual(cfg,data);
        %remove muscle artifacts
        disp('REMOVE MUSCL ARTEFACT TRIALS')
        cfg                     = [];
        cfg                     = dummy.cfg;
        data                    = ft_rejectartifact(cfg,data);
        clear dummy;
    end

    %% ---------------- Task 3 (voltage cleaning) -----------------------------    
    if find(cfgx.tasks==3)
        %find bad channels
        disp('(Step 1 of 2) SELECT AND REMOVE BAD TRIALS -->> CHANNELS ARE KEPT <<--')
        
        data_eeg=data; 
        %run artefact rejection on EEG data only
        cfg                     = [];
        cfg.method              = 'summary';%'channel';
        cfg.keepchannel         = 'yes';
        cfg.metric              = 'absmax';%'var';
        cfg.layout              = lay;     % this allows for plotting
        cfg.channel             = cfgx.eeg;   % only EEG channels
        [data_eeg] = ft_rejectvisual(cfg,data_eeg);
     
        
        clear data;
        
        disp('(Step 2 of 2) SELECT AND REMOVE BAD CHANNELS -->> WILL BE REMOVED <<--')
                
        %run artefact rejection on EEG data only
        cfg                     = [];
        cfg.method              = 'summary';%'channel';
        cfg.keepchannel         = 'no';
        cfg.metric              = 'absmax';%'var';
        cfg.layout              = lay;     % this allows for plotting
        cfg.channel             = 'all';   % only EEG channels
        [data_eeg]    = ft_rejectvisual(cfg,data_eeg);
        data=data_eeg;
        clear data_eeg data_extra;
        
    end

    %% ---------------- Task 4 (ICA) -----------------------------
    if find(cfgx.tasks==4)
       
        % run ica
        disp('RUN ICA') 
        cfg               = [];
        cfg.channel       = cfgx.eeg; 
        cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
        ic_data           = ft_componentanalysis(cfg,data);
        ic_data = rmfield(ic_data,'sampleinfo');
        if isfield(ic_data,'trialinfo')
            ic_data = rmfield(ic_data,'trialinfo');
        end

        % inspect and select components to reject
        disp('SELECT ICA COMPONENTS TO REMOVE')
        cfg               = [];
        cfg.viewmode      = 'component';
        cfg.continuous    = 'no';
        cfg.compscale     = 'local';
        cfg.layout        = lay;
        cfg.verticalpadding = 0.1;
        ft_databrowser(cfg,ic_data) 
        
        % chose components to reject
        disp('REMOVE ICA COMPONENTS')
        Prompt(1,:) = {'Enter all components to be rejected (space-separated - e.g. 1 4 7 12)', 'comp',[]};
        name='Component Rejection';
        formats(1,1).type = 'edit';
        formats(1,1).format = 'text';
        formats(1,1).limits = [0 1];
        DefAns = struct([]);
        DefAns(1).comp = 'None';
        [answer,canceled] = inputsdlg(Prompt,name,formats,DefAns);
        fprintf('removing components: %s \n', answer.comp)
        comp_rej = unique(str2num(answer.comp));
        close all
    
        % remove components chosen
        cfg = [];
        cfg.component = comp_rej;
        data_clean=ft_rejectcomponent(cfg,ic_data);
        
        % save both ica and data clean
        data = {data_clean, ic_data, comp_rej};
        
    end
    
    %% ---------------- Task 5 (Visual) -----------------------------
    if find(cfgx.tasks==5)
        
        % visual inspection of raw data and manual data rejection
        disp('VISUAL INSPECTION OF THE DATA')
        data_show = rmfield(data,'sampleinfo');
        cfg                 = [];
        cfg.viewmode        = 'butterfly';       
        cfg.ylim            = 'maxmin';       
        ft_databrowser(cfg,data_show);
       
        % chose trials to reject
        disp('REMOVE TRIALS ')
        Prompt(1,:) = {'Enter all trials to be rejected (space-separated - e.g. 1 4 7 12)', 'trial',[]};
        name='Trials Rejection';
        formats(1,1).type = 'edit';
        formats(1,1).format = 'text';
        formats(1,1).limits = [0 1];
        DefAns = struct([]);
        DefAns(1).trial = 'None';
        [answer,canceled] = inputsdlg(Prompt,name,formats,DefAns);
        fprintf('removing trials: %s \n', answer.trial)
        trials_rej = unique(str2num(answer.trial));
        trials_not_rej = 1:length(data.trial);
        trials_not_rej(trials_rej) = [];
        
        % remove trials
        cfg = [];
        cfg.trials = trials_not_rej;
        data = ft_preprocessing(cfg, data);
        
        % save both data and trials rejected
        data = {data, trials_rej};

        
    end

end

