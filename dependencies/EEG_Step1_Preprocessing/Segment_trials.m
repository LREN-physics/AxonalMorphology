function [dataVout]=Segment_trials(dataV,nchans,ind,cond, baseline,poststim,sr)
% Function to segment trials with a given baseline and post-stimulus period


%% Loop for the onset indices for each trial 
for kk=1:length(ind)    
    
    % define trial
    dataVout.trial{kk}=dataV.trial{1}(1:nchans,ind(kk)-baseline:ind(kk)+poststim);    
    
    % define time for that trial
    step=(poststim/sr+baseline/sr)/(poststim+baseline);
    dataVout.time{kk}=[(-baseline)/sr:step:poststim/sr];
end

% define other things for that structure
dataVout.hdr.Fs = sr;
dataVout.hdr.nChans = size(dataVout.trial{1},1);
dataVout.hdr.nSamples = size(dataVout.trial{1},2);
dataVout.hdr.nTrials = size(dataVout.trial,2);
dataVout.hdr.label = dataV.label;
dataVout.label=dataV.label;
dataVout.fsample=dataV.fsample;
dataVout.cfg=dataV.cfg;
dataVout.sampleinfo=[ind-baseline; ind+poststim]';
dataVout.trialinfo=[ind; ones(1,length(ind))*cond]'; %cond is a label that identify the condition 

end