function [data] = Interpolate_BadChannels(cfg,data)
%function [data] = my_interpolateBadChannels3(cfg,data)
%
%   MY_INTERPOLATEBADCHANNELS identifies missing channels in data structure
%   using a channel label template layout and interpolates these channels
%   useing nearest neighbor interpolation
%
%   example configuration structure:
%     cfg                      = [];
%     cfg.bad_channels         = {}; %% channel label
%     cfg.layout_channel_names = 'gtec_63chans.mat';
%     cfg.layout               = lay;
%
%   version 4 : MDL 11/07/2019; using two interpolation steps
%   version 2 : CP 14/3/2016; using custom neighbour structure for gtec
%   version '': CP 4/1/2016; update fix missing channels 2/2/2016

'RUNNING MY_INTERPOLATEBADCHANNELS FUNCTION'

cfgx = cfg;
clear cfg;

lay = cfgx.layout;

%find missing channels (removed by hand) in dataset based on list of
%desired number of template channels
missingchans = {};
channames=lay.label(1:end-2); 
la=channames;
ld=data.label;
for i=1:length(la)
    if ~strcmp(la{i},ld)
        missingchans=[missingchans,{la{i}}];
    end;
end;
clear la ld;

%update apriori bad channel structure (remove from this list the
%missing channels identified by hand)
badchans = cfgx.bad_channels;

for i = 1:length(badchans)
    if ~isstr(badchans{i})
        badchans(i) = channames(badchans{i});
    end;
end;

%if bad channels were specified
if ~isempty(badchans) && ~isempty(missingchans)
    %index vector initially keeping all bad channels as specified
    ind = ones(1,length(badchans));
    %loop over bad channels
    for i = 1:length(badchans)
        %if one of the channels matches the missing channels, exclude
        %it from the bad missing channel array
        if find(strcmp(missingchans,badchans(i)))
            %by setting the index vector to zero
            ind(i)=0;
        end;
    end;
    %apply index vector to update bad channel vector
    badchans=badchans(logical(ind));
    clear ind;
end;

%put all channels to be interpolated in one list
allchans = [badchans,missingchans];


if ~isempty(allchans)
    %find neighbours
    cfg=[];
    cfg.method         = 'triangulation';%, 'triangulation' or 'template'
    cfg.layout         = lay;%filename of the layout, see FT_PREPARE_LAYOUT
    cfg.channel        = channames;%channels for which neighbours should be found
    cfg.feedback       = 'yes';
    neighbours         = ft_prepare_neighbours(cfg);
    %insert missing channels
    cfg = [];
    cfg.method         = 'nearest';%, 'average', 'spline' or 'slap' (default='nearest')
    cfg.layout         = lay;
    if ~isempty(badchans)
        cfg.badchannel     = badchans;%cell-array, see FT_CHANNELSELECTION for details
    end;
    if ~isempty(missingchans)
        cfg.missingchannel = missingchans;%cell-array, see FT_CHANNELSELECTION for details
    end;
    cfg.neighbours     = neighbours;%neighbourhoodstructure, see also FT_PREPARE_NEIGHBOURS
    cfg.trials         = 'all';% or a selection given as a 1xN vector (default = 'all')
    data               = ft_channelrepair(cfg, data);
    clear neighbours;
end;

%% sort channels according to template
channames;

n = length(channames);
ind = nan(n,1);
for i = 1:n
    try
        ind(i)=find(strcmpi(data.label,channames(i)));
    end;
end;

check_chans_missing = sum(isnan(ind));

if check_chans_missing > 0
    %error('error: missing channels!') % changed on the 11th of July,
    %in case the interpolation did not succed for all channels because
    % all neighbours of one or more bad channels were all to be
    % inteporlated
    % in this new step, we repeat the interpolation by treating the
    % missing channels as channels to be interpolated on the
    % interpolated data
    
       cfg=[];
    cfg.method         = 'triangulation';%, 'triangulation' or 'template'
%     if strcmpi(cfgx.layout_channel_names,'gtec_63chans.mat')
%         cfg.template       = 'gtec63_neighbours.mat';%name of the template file, e.g. CTF275_neighb.mat
%         %'easycap128ch-avg_neighb.mat';
%     elseif strcmpi(cfgx.layout_channel_names,'gtec_gamma_62+3chans.mat')
%         cfg.template       = 'gtec62_neighbours.mat';%name of the template file, e.g. CTF275_neighb.mat
%         %'easycap128ch-avg_neighb.mat';
%     end;
    cfg.layout         = lay;%filename of the layout, see FT_PREPARE_LAYOUT
    cfg.channel        = channames(isnan(ind));%channels for which neighbours could not be found in the first step
    neighbours         = ft_prepare_neighbours(cfg);
    
    
    %%%%%%
    
    missingchans        = channames(isnan(ind));%channels for which neighbours could not be found in the first step
    
    %insert missing channels
   
    cfg = [];
    cfg.method         = 'nearest';%, 'average', 'spline' or 'slap' (default='nearest')
    cfg.layout         = lay;
   
    cfg.missingchannel = missingchans;%cell-array, see FT_CHANNELSELECTION for details
   
    cfg.neighbours     = neighbours;%neighbourhoodstructure, see also FT_PREPARE_NEIGHBOURS
    cfg.trials         = 'all';% or a selection given as a 1xN vector (default = 'all')
    data               = ft_channelrepair(cfg, data);
    clear neighbours;
        
   
else
    data.label=data.label(ind);
    n=length(data.trial);
    for i = 1:n
        data.trial(i) = {data.trial{i}(ind,:)};
    end;
end;

%% sort channels according to template
channames;

n = length(channames);
ind = nan(n,1);
for i = 1:n
    try
        ind(i)=find(strcmpi(data.label,channames(i)));
    end;
end;

check_chans_missing2 = sum(isnan(ind));

if check_chans_missing2 > 0
    error('error: missing channels!')
else
    data.label=data.label(ind);
    n=length(data.trial);
    for i = 1:n
        data.trial(i) = {data.trial{i}(ind,:)};
    end;
end;

end
