function [Data,TriggersTime,TriggersName]=Read_trcdata(pathdata,datafilename)
% Reads EEG files with extension .TRC
%
%   INPUTS: 
%       pathdata (string)          - path to the subjects EEG data folder
%       datafilename (string)       - name of the file we want to read
%
%   OUPUTS:
%       Data (matrix)               - data values 
%       TriggersTime (vector)       - list of original trigger times read from the data       
%       TriggersName (vector)       - list of original trigger names read from the data
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

% Read data
cd(pathdata)
cfg = [];
cfg.dataset     = datafilename;
Data       = ft_preprocessing(cfg);

% Read events
events = ft_read_event(datafilename);
temp = [events.sample];
TriggersTime =  temp(1,:);

temp = [events.value];
TriggersName =  temp(1,:);

end





