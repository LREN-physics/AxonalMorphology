function [TriggersName,TriggersTime] = RenameTriggers_FirstTriggerRemoval(TriggersName_orig,TriggersTime_orig,path_data,sub_name,run)
% Function that:
% Part 1: renames the numerical triggers to categorical triggers
% Part 2: removes the first unwanted triggers
% Gets called in Preprocessing_FieldTrip.m 
%
%   INPUTS: 
%       TriggersName_orig (vector)  - list of original trigger names read from the data
%       TriggersTime_orig (vector)  - list of original trigger times read from the data
%       path_data (string)          - path to the subjects EEG data folder
%       sub_name (string)           - name of the subject to analyse
%       run (number)                - number of the block/run to analyse
%
%   OUPUTS:
%       TriggersName (cell array)   - list of original trigger names read from the data
%       TriggersTime (vector)       - list of original trigger times read from the data
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

%% Part 1

% Because the first one is the start of the block
TriggersTime = TriggersTime_orig(2:end)';
TriggersName_temp = TriggersName_orig(2:end)';

% Rename according to the triggers sent
TriggersName = cell(length(TriggersName_temp),1);
for i=1:length(TriggersName_temp)
    TriggersName{i} = num2str(TriggersName_temp(i));
    if strcmp(TriggersName{i},'1')
        TriggersName{i} = 'Left';
    elseif strcmp(TriggersName{i},'2')
        TriggersName{i} = 'Right';
    elseif strcmp(TriggersName{i},'4')
        TriggersName{i} = 'Catch';
    end
end

%% Part 2

% Remove other stupid characters that could have been saved
if contains(TriggersName,'"')
    TriggersName = erase(TriggersName,'"');
end


% Load data
cd(path_data)
cd ..
cd('Matlab_Outputs')

load(sprintf('%s_timestamps_block_%d.mat',sub_name,run), 'timestamps_check')
load(sprintf('%s_input_configuration.mat',sub_name))
load(sprintf('%s_output_block_%d.mat',sub_name,run))

% Identify triggers that did not passed control quality of the visual
% presentation & remove triggers where subject answer did not meet criteria
% (i.e out.RESPONSE = -1)
for i=1:length(TriggersName)
    if ((strcmp(TriggersName{i},'Left') || strcmp(TriggersName{i},'Right')) && abs(timestamps_check.Time_Dur(i)-inputs.DURATION_STIM*10^3) >= 3)  
        TriggersName{i} = 'REMOVE';
    elseif ((strcmp(TriggersName{i},'Left') || strcmp(TriggersName{i},'Right')) && out.RESPONSE(i)==-1)
        TriggersName{i} = 'REMOVE2';
    end
    
end

% Make plot for quality control
figure;
scatter(1:length(TriggersName_temp),TriggersName_temp,'k','Filled')
hold on
aux = TriggersName_temp;
aux(strcmp(TriggersName,'REMOVE')) = 0;
scatter(1:length(TriggersName_temp),aux,'r');
xlim([0, length(TriggersName_temp)])
title('Triggers used')
ylabel('Value of trigger')
xlabel('Trigger')
legend('All triggers','Triggers used')
cd(path_data)
saveas(gcf,sprintf('Triggers_used_block_%d.png',run))
disp('Number of trials removed where participant response did not meet criteria:')
disp(sum(strcmp(TriggersName,'REMOVE2')))

% Remove those unwanted triggers
TriggersTime=TriggersTime(~strcmp(TriggersName,'REMOVE'));
TriggersName=TriggersName(~strcmp(TriggersName,'REMOVE'));
TriggersTime=TriggersTime(~strcmp(TriggersName,'REMOVE2'));
TriggersName=TriggersName(~strcmp(TriggersName,'REMOVE2'));


end





