function data = RenameElectrodesMicromed(data)
% Function to rename the electrodes obtained with the Micromed system. This
% is because at CHUV we cannot rename the electrodes in the EEG system itself, 
% due to limit in space/caracters in the Micromed system.
%
%   INPUTS: 
%       data (structure)            - structure of the EEG data
%
%   OUPUTS:
%       data (structure)            - structure of the EEG data corrected
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

idx = find(strcmp(data.label, 'FTT10'));
if ~isempty(idx)
    data.label{idx} = 'FTT10H';
end

idx = find(strcmp(data.label, 'TPP10'));
if ~isempty(idx)
    data.label{idx} = 'TPP10H';
end

idx = find(strcmp(data.label, 'POO1H'));
if ~isempty(idx)
    data.label{idx} = 'POO10H';
end

idx = find(strcmp(data.label, 'PPO1H'));
if ~isempty(idx)
    data.label{idx} = 'PPO10H';
end

idx = find(strcmp(data.label, 'PPO10'));
if ~isempty(idx)
    data.label{idx} = 'PPO1';
end

idx = find(strcmp(data.label, '01'));
if ~isempty(idx)
    data.label{idx} = 'O1';
end

idx = find(strcmp(data.label, '02'));
if ~isempty(idx)
    data.label{idx} = 'O2';
end

idx = find(strcmp(data.label, '0Z'));
if ~isempty(idx)
    data.label{idx} = 'OZ';
end


end



