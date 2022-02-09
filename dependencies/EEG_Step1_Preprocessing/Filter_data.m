function [outputdata] = Filter_data(inputdata, cnt, type)
%function [data_filt] = filter_data(data, low, high, sr, type)
%
%data: chan x sp
%low: low bound of the bandpass filter
%high: hig bound of the bandpass filter
%sr: sampling rate
%type: different functions to be used 1=eeglab 2=fieldtrip

data                    = inputdata.trial{1}; 
low                     = cnt.bandpass(1);
high                    = cnt.bandpass(2);
sr                      = inputdata.fsample;

[nchan,nsp] = size(data);

filtorder               = 3*fix(sr/low);
filtorder = 33792;

%% eeglab
if type == 1

    revfilt                 = []; % 
    usefft                  = 0;


    eeglab;
    
    EEGinput                     = eeg_emptyset;
    EEGinput.data                = data;
    EEGinput.times               = [1:nsp].*(1/sr);
    EEGinput.nbchans             = nchan;
    EEGinput.srate               = sr;
    EEGinput                     = eeg_checkset(EEGinput);
    EEGinput.setname             = 'test';
    EEGinput                     = eeg_checkset(EEGinput);

    clear neeg ieeg neog ieog necg iecg;

    EEGinput.chanlocs            = [];
    EEGinput.trials = 1;
    EEGinput.event = [];
    EEGinput.event(1).eventtype = 'boundary';
    EEGinput.event(1).eventlatency  = 1;
    EEGinput.event(1).eventduration = 0;
    EEGinput.pnts = nsp;
    
    %EEGinput = pop_eegfiltnew(EEGinput,  49, 51, filtorder, 1, usefft);
    EEGout = pop_eegfiltnew(EEGinput,  low, high, filtorder, revfilt, usefft); 
    

    
    data_filt = EEGout.data;
    clear EEGinput EEGout;

%% fieldtrip
elseif type == 2

    filtertype              = 'fir';
    dir                     = 'twopass';
    instabilityfix          = 'no';
    wintype                 = 'hamming';

    [data_filt]                  = ft_preproc_bandpassfilter(data, sr, [low high], filtorder, filtertype, dir, instabilityfix, wintype);

end

outputdata.trial{1}         = double(data_filt);
outputdata.label            = inputdata.label;
outputdata.fsample          = inputdata.fsample;
outputdata.cfg              = inputdata.cfg;
outputdata.hdr              = inputdata.hdr;
outputdata.time             = inputdata.time;

end
