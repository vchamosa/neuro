%% For use after go/no-go behavioural analysis
%
% This script plots behavioural parameters across go/no-go training
% sessions of mice selected by the user. Note that the data must have been
% previously extracted using GoNoGo_extract.m
% Author: Victor Chamosa Pino
%

%% Load previously analysed mice data

mice = selectanimals('task','GoNoGo','responsible_person','Victor');
%LEARNED
mice = {'C57-1024' 'C57-1025' 'C57-1026' 'C57-259' 'C57-483' 'C57-485' 'C57-494' 'C57-745' 'C57-748' 'C57-749' 'MCos1164' 'MCos1165' 'MCos1166' 'MCos1497' 'MCos1498' 'MCos1499'};

% exclude mice here
mice_unwanted = horzcat(selectanimals('strain', 'Thy1G'),selectanimals('strain', 'Thy1GCaMP'),selectanimals('strain', 'Thy1GCaMP6s'));
mice = setdiff(mice, mice_unwanted);
for iMouse = {'MCos1501' 'MCos170701_B' 'MCos170701_D' 'MCos170701_F' 'MCos170801_C' 'MCos170801_D' 'MCos170801_F' 'MCos402' 'MCos403' 'MCos216' 'MCos217' 'MCos863' 'MCos864' 'MCos904' 'MCos905' 'MCos906'}
    % 'MCos1522' 'MCos1523' 'MCos1524' 'MCos1525' 'MCos1500' 'MCos1164' 'MCos1165' 'MCos1166' 'MCos1167'
    mice = mice(~strcmp(mice, iMouse));
end
numMice = numel(mice);

% Load the events saved in the mat 
metadata = getmousemetadata(mice);
Data = cell(1, numMice);
for iMouse = 1 : numMice
    Data{iMouse} = load([mice{iMouse} '_analysed.mat']);
end


%% Plotting 
% from here onward you can analyse different behavioural parameters
% If adding new analysis modalities, keep in mind that some functions might
% exist already - check Functions/PlottingFunctions/


%% Hits

[numHits,numHits_norm] = calc_numHits(Data,metadata);


%% Correct rejections

[numCRs,numCRs_norm] = calc_numCRs(Data,metadata);


%% Hit reaction time

RT_hits = calc_hitRT(Data,metadata);


%% False alarm reaction time
% not a hugely interesting metric
%RT_FAs = calc_FART(Data,metadata);


%% d'

[hitRate,falseAlarmRate,dPrime] = calc_dPrime(Data,metadata,'responsible_person','Victor');


%% Percentage of correct trials

percentCorrect = calc_percor(Data,metadata,'responsible_person','Victor');


%% Optimal window
% mind that it currently calculates optimal window for correction trials too
[opt_rew,opt_time,opt_percor,learned,learninglength] = calc_optwin(Data,metadata,'responsible_person','Victor');


%% Other plots
% Includes total rewards, hit & CR percent, IGMs, total pushes, and weight progression (commented out)
plot_others(Data,metadata,numHits,numHits_norm,numCRs,numCRs_norm,hitRate,falseAlarmRate,percentCorrect,'responsible_person','Victor')


%% Comparison plots
% can only take the category you are comparing - other specifications 
% please enact above when loading data
%
%


%% Hits comparison

%numHits_cmp = calc_numHits_cmp(Data,metadata,'responsible_person',{'Michelle','Victor'});


%% CRs comparison

%numCRs_cmp = calc_numCRs_cmp(Data,metadata,'responsible_person',{'Michelle','Victor'});


%% d' comparison

[dPrime_cmp,hitRate_cmp,falseAlarmRate_cmp,p_dPrime] = calc_dPrime_cmp(Data,metadata,'responsible_person',{'Michelle','Victor'});


%% Percentage correct comparison
% please ensure inputs to compare are in the same order as in the d' plot
[percor_cmp,p_percor] = calc_percor_cmp(Data,metadata,hitRate_cmp,'responsible_person',{'Michelle','Victor'});


%% Optimal window comparison

[optRew_cmp,optTime_cmp,optPercor_cmp,learned_cmp,learninglength_cmp,p_optTime] = calc_optwin_cmp(Data,metadata,'responsible_person',{'Michelle','Victor'});


%% Save figures?

qtn = input('Would you like to save any figures? y/n: ','s');
if  strcmpi(qtn,'y')
    dir = input('Please declare a folder to save them in ','s');
    figs = get(0,'children');
    for s = 1:length(figs)
        if length(figs(s).Visible) == 2
            figure(figs(s))
            filename = input('If you would like to save this figure write its filename, else press enter ','s');
            if filename ~= 0
                saveas(figs(s),fullfile(dir,filename),'svg')
            else continue
            end
        end
    end
end

close all


