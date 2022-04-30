%% For use after GoNoGo behavioural analysis
%
% This script plots behavioural parameters across go/no-go training
% sessions of mice selected by the user. Note that the data must have been
% previously extracted using GoNoGo_extract.m
%

%% Load previously analysed mice data

%mice = selectanimals('task','GoNoGo','responsible_person','Victor');
mice = {'C57-483' 'C57-485'}; 

% exclude mice here
mice_unwanted = horzcat(selectanimals('strain', 'Thy1G'),selectanimals('strain', 'Thy1GCaMP'),selectanimals('strain', 'Thy1GCaMP6s'));
mice = setdiff(mice, mice_unwanted);
for iMouse = { 'MCos170701_B' 'MCos170701_D' 'MCos170701_F' 'MCos170801_C' 'MCos170801_D' 'MCos170801_F' 'MCos402' 'MCos403' 'MCos216' 'MCos217' 'MCos863' 'MCos864' 'MCos904' 'MCos905' 'MCos906'}
    % 'MCos1522' 'MCos1523' 'MCos1524' 'MCos1525' 'MCos1500' 'MCos1164' 'MCos1165' 'MCos1166' 'MCos1167'
    mice = mice(~strcmp(mice, iMouse));
end
numMice = numel(mice);

% Load data & metadata
metadata = getmousemetadata(mice);
Data = cell(1, numMice);
tMetadata = cell(1, numMice);
for iMouse = 1 : numMice
    Data{iMouse} = load([mice{iMouse} '_analysed.mat']);
    tMetadata{iMouse} = loadjson([mice{iMouse} '.json']);
end


%% Compute variables

cutoff = 1.5;
window = 600000; %length of the window = 10 min
step = 15000; % 15 sec

% Determine laser & control trials
s_controlHits = cell(numMice,1);
s_laserHits = cell(numMice,1);
s_controlMisses = cell(numMice,1);
s_laserMisses = cell(numMice,1);
s_controlCRs = cell(numMice,1);
s_laserCRs = cell(numMice,1);
s_controlFAs = cell(numMice,1);
s_laserFAs = cell(numMice,1);
d_controlHits = cell(numMice,1);
d_laserHits = cell(numMice,1);
d_controlMisses = cell(numMice,1);
d_laserMisses = cell(numMice,1);
d_controlCRs = cell(numMice,1);
d_laserCRs = cell(numMice,1);
d_controlFAs = cell(numMice,1);
d_laserFAs = cell(numMice,1);
s_windows_hit_control = cell(numMice,1);
s_windows_miss_control = cell(numMice,1);
s_windows_CR_control = cell(numMice,1);
s_windows_FA_control = cell(numMice,1);
s_windows_hit_laser = cell(numMice,1);
s_windows_miss_laser = cell(numMice,1);
s_windows_CR_laser = cell(numMice,1);
s_windows_FA_laser = cell(numMice,1);
d_windows_hit_control = cell(numMice,1);
d_windows_miss_control = cell(numMice,1);
d_windows_CR_control = cell(numMice,1);
d_windows_FA_control = cell(numMice,1);
d_windows_hit_laser = cell(numMice,1);
d_windows_miss_laser = cell(numMice,1);
d_windows_CR_laser = cell(numMice,1);
d_windows_FA_laser = cell(numMice,1);
s_optGO_proportions = cell(numMice,1);
s_optNOGO_proportions = cell(numMice,1);
d_optGO_proportions = cell(numMice,1);
d_optNOGO_proportions = cell(numMice,1);
for mus = 1:numMice
    s_controlHits{mus} = cell(length(Data{mus}.Events),1);
    s_laserHits{mus} = cell(length(Data{mus}.Events),1);
    s_controlMisses{mus} = cell(length(Data{mus}.Events),1);
    s_laserMisses{mus} = cell(length(Data{mus}.Events),1);
    s_controlCRs{mus} = cell(length(Data{mus}.Events),1);
    s_laserCRs{mus} = cell(length(Data{mus}.Events),1);
    s_controlFAs{mus} = cell(length(Data{mus}.Events),1);
    s_laserFAs{mus} = cell(length(Data{mus}.Events),1);
    d_controlHits{mus} = cell(length(Data{mus}.Events),1);
    d_laserHits{mus} = cell(length(Data{mus}.Events),1);
    d_controlMisses{mus} = cell(length(Data{mus}.Events),1);
    d_laserMisses{mus} = cell(length(Data{mus}.Events),1);
    d_controlCRs{mus} = cell(length(Data{mus}.Events),1);
    d_laserCRs{mus} = cell(length(Data{mus}.Events),1);
    d_controlFAs{mus} = cell(length(Data{mus}.Events),1);
    d_laserFAs{mus} = cell(length(Data{mus}.Events),1);
    s_windows_hit_control{mus} = cell(length(Data{mus}.Events),1);
    s_windows_miss_control{mus} = cell(length(Data{mus}.Events),1);
    s_windows_CR_control{mus} = cell(length(Data{mus}.Events),1);
    s_windows_FA_control{mus} = cell(length(Data{mus}.Events),1);
    s_windows_hit_laser{mus} = cell(length(Data{mus}.Events),1);
    s_windows_miss_laser{mus} = cell(length(Data{mus}.Events),1);
    s_windows_CR_laser{mus} = cell(length(Data{mus}.Events),1);
    s_windows_FA_laser{mus} = cell(length(Data{mus}.Events),1);
    d_windows_hit_control{mus} = cell(length(Data{mus}.Events),1);
    d_windows_miss_control{mus} = cell(length(Data{mus}.Events),1);
    d_windows_CR_control{mus} = cell(length(Data{mus}.Events),1);
    d_windows_FA_control{mus} = cell(length(Data{mus}.Events),1);
    d_windows_hit_laser{mus} = cell(length(Data{mus}.Events),1);
    d_windows_miss_laser{mus} = cell(length(Data{mus}.Events),1);
    d_windows_CR_laser{mus} = cell(length(Data{mus}.Events),1);
    d_windows_FA_laser{mus} = cell(length(Data{mus}.Events),1);
    s_optGO_proportions{mus} = cell(length(Data{mus}.Events),1);
    s_optNOGO_proportions{mus} = cell(length(Data{mus}.Events),1);
    d_optGO_proportions{mus} = cell(length(Data{mus}.Events),1);
    d_optNOGO_proportions{mus} = cell(length(Data{mus}.Events),1);
    for day = 1:numel(Data{mus}.Events)
        if ~isfield(tMetadata{mus}{day},'depth')
            continue
        else
            sessionLength = Data{mus}.Events{day}.dur;
            % Resolve unmatched/unfinished laser binary changes if they exist
            if isstring(Data{mus}.Events{day}) % account for missing sessions
                Data{mus}.Events{day} = [];
                Data{mus}.Events{day}.LaserOn = [];
                Data{mus}.Events{day}.LaserOff = [];
            end
            laserOn_TimesALL = Data{mus}.Events{day}.LaserOn;
            laserOff_TimesALL = Data{mus}.Events{day}.LaserOff;
            laserOn_Times = laserOn_TimesALL;
            laserOff_Times = laserOff_TimesALL;
            if length(laserOff_TimesALL) < length(laserOn_TimesALL)
                if laserOn_TimesALL(1) <= 99999 % eliminate initial weirdism
                    laserOn_Times = laserOn_TimesALL(laserOn_TimesALL > 99999);
                end
                if length(laserOff_Times) < length(laserOn_Times) % if still uneven, account for initial control trial
                    laserOff_Times = [0; laserOff_TimesALL];
                end
            elseif length(laserOn_TimesALL) < length(laserOff_TimesALL)
                if laserOff_TimesALL(1) <= 99999 % eliminate initial weirdism
                    laserOff_Times = laserOff_TimesALL(laserOff_TimesALL > 99999);
                end
                if length(laserOn_Times) < length(laserOff_Times) % if still uneven, eliminate unfinished trial at the end
                    laserOff_Times = laserOff_Times(1:(end - 1));
                end
            end
            if length(laserOn_Times) == length(laserOff_Times)
                if laserOff_Times(1) <= 99999 % eliminate initial weirdism
                    laserOff_Times(1) = 0; % account for initial control trial
                elseif laserOn_Times(1) <= 99999 % eliminate initial weirdism
                    laserOn_Times = laserOn_Times(laserOn_Times > 99999);
                    laserOff_Times = laserOff_Times(1:(end - 1)); % eliminate unfinished trial at the end
                end
                if laserOn_Times(1) < laserOff_Times(1)
                    laserOff_Times = [0; laserOff_Times]; % account for initial control trial
                    laserOff_Times = laserOff_Times(1:(end - 1)); % and for the unfinished trial at the end
                end
            end
            laserOn_fullTrials = [];
            laserOff_fullTrials = [];
            for lpc = 1: length(laserOn_Times)
                laser_length(lpc) = laserOn_Times(lpc) - laserOff_Times(lpc);
                if laser_length(lpc) > (4*sessionLength/(40*60))
                    laserOn_fullTrials = [laserOn_fullTrials laserOn_Times(lpc)];
                    laserOff_fullTrials = [laserOff_fullTrials laserOff_Times(lpc)];
                end
            end
            % Analyse dendritic sessions
            if tMetadata{mus}{day}.depth == "dendritic"
                % Separate control and laser trials
                if tMetadata{mus}{day}.trial == "go"
                    d_controlCRs{mus}{day} = Data{mus}.Events{day}.correctRejections';
                    d_controlFAs{mus}{day} = Data{mus}.Events{day}.falseAlarms';
                    d_controlHits{mus}{day} = NaN(length(Data{mus}.Events{day}.hits),1);
                    d_laserHits{mus}{day} = NaN(length(Data{mus}.Events{day}.hits),1);
                    for t = 1:length(Data{mus}.Events{day}.hits) % Hits
                        for L = 1:length(laserOff_fullTrials)
                            if Data{mus}.Events{day}.hits(t) <= laserOff_fullTrials(1)
                                d_controlHits{mus}{day}(t) = Data{mus}.Events{day}.hits(t);
                            elseif Data{mus}.Events{day}.hits(t) >= laserOff_fullTrials(L) && Data{mus}.Events{day}.hits(t) <= laserOn_fullTrials(L)
                                d_controlHits{mus}{day}(t) = Data{mus}.Events{day}.hits(t);
                            end
                        end
                        if d_controlHits{mus}{day}(t) == 0
                            d_laserHits{mus}{day}(t) = Data{mus}.Events{day}.hits(t);
                        end
                    end
                    d_controlMisses{mus}{day} = NaN(length(Data{mus}.Events{day}.misses),1);
                    d_laserMisses{mus}{day} = NaN(length(Data{mus}.Events{day}.misses),1);
                    for t = 1:length(Data{mus}.Events{day}.misses) % Misses
                        for L = 1:length(laserOff_fullTrials)
                            if Data{mus}.Events{day}.misses(t) <= laserOff_fullTrials(1)
                                d_controlMisses{mus}{day}(t) = Data{mus}.Events{day}.misses(t);
                            elseif Data{mus}.Events{day}.misses(t) >= laserOff_fullTrials(L) && Data{mus}.Events{day}.misses(t) <= laserOn_fullTrials(L)
                                d_controlMisses{mus}{day}(t) = Data{mus}.Events{day}.misses(t);
                            end
                        end
                        if d_controlMisses{mus}{day}(t) == 0
                            d_laserMisses{mus}{day}(t) = Data{mus}.Events{day}.misses(t);
                        end
                    end
                elseif  tMetadata{mus}{day}.trial == "no-go"
                    d_controlHits{mus}{day} = Data{mus}.Events{day}.hits';
                    d_controlMisses{mus}{day} = Data{mus}.Events{day}.misses';
                    d_controlCRs{mus}{day} = NaN(length(Data{mus}.Events{day}.correctRejections),1);
                    d_laserCRs{mus}{day} = NaN(length(Data{mus}.Events{day}.correctRejections),1);
                    for t = 1:length(Data{mus}.Events{day}.correctRejections) % CR
                        for L = 1:length(laserOff_fullTrials)
                            if Data{mus}.Events{day}.correctRejections(t) <= laserOff_fullTrials(1)
                                d_controlCRs{mus}{day}(t) = Data{mus}.Events{day}.correctRejections(t);
                            elseif Data{mus}.Events{day}.correctRejections(t) >= laserOff_fullTrials(L) && Data{mus}.Events{day}.correctRejections(t) <= laserOn_fullTrials(L)
                                d_controlCRs{mus}{day}(t) = Data{mus}.Events{day}.correctRejections(t);
                            end
                        end
                        if d_controlCRs{mus}{day}(t) == 0
                            d_laserCRs{mus}{day}(t) = Data{mus}.Events{day}.correctRejections(t);
                        end
                    end
                    d_controlFAs{mus}{day} = NaN(length(Data{mus}.Events{day}.falseAlarms),1);
                    d_laserFAs{mus}{day} = NaN(length(Data{mus}.Events{day}.falseAlarms),1);
                    for t = 1:length(Data{mus}.Events{day}.falseAlarms) % FA
                        for L = 1:length(laserOff_fullTrials)
                            if Data{mus}.Events{day}.falseAlarms(t) <= laserOff_fullTrials(1)
                                d_controlFAs{mus}{day}(t) = Data{mus}.Events{day}.falseAlarms(t);
                            elseif Data{mus}.Events{day}.falseAlarms(t) >= laserOff_fullTrials(L) && Data{mus}.Events{day}.falseAlarms(t) <= laserOn_fullTrials(L)
                                d_controlFAs{mus}{day}(t) = Data{mus}.Events{day}.falseAlarms(t);
                            end
                        end
                        if d_controlFAs{mus}{day}(t) == 0
                            d_laserFAs{mus}{day}(t) = Data{mus}.Events{day}.falseAlarms(t);
                        end
                    end
                end
                
                d_controlHits{mus}{day} = d_controlHits{mus}{day}(~isnan(d_controlHits{mus}{day}));
                d_laserHits{mus}{day} = d_laserHits{mus}{day}(~isnan(d_laserHits{mus}{day}));
                d_controlMisses{mus}{day} = d_controlMisses{mus}{day}(~isnan(d_controlMisses{mus}{day}));
                d_laserMisses{mus}{day} = d_laserMisses{mus}{day}(~isnan(d_laserMisses{mus}{day}));
                d_controlCRs{mus}{day} = d_controlCRs{mus}{day}(~isnan(d_controlCRs{mus}{day}));
                d_laserCRs{mus}{day} = d_laserCRs{mus}{day}(~isnan(d_laserCRs{mus}{day}));
                d_controlFAs{mus}{day} = d_controlFAs{mus}{day}(~isnan(d_controlFAs{mus}{day}));
                d_laserFAs{mus}{day} = d_laserFAs{mus}{day}(~isnan(d_laserFAs{mus}{day}));
            % Analyse somatic sessions
            elseif tMetadata{mus}{day}.depth == "somatic"
                % Separate control and laser trials
                if tMetadata{mus}{day}.trial == "go"
                    s_controlCRs{mus}{day} = Data{mus}.Events{day}.correctRejections';
                    s_controlFAs{mus}{day} = Data{mus}.Events{day}.falseAlarms';
                    s_controlHits{mus}{day} = NaN(length(Data{mus}.Events{day}.hits),1);
                    s_laserHits{mus}{day} = NaN(length(Data{mus}.Events{day}.hits),1);
                    for t = 1:length(Data{mus}.Events{day}.hits) % Hits
                        for L = 1:length(laserOff_fullTrials)
                            if Data{mus}.Events{day}.hits(t) <= laserOff_fullTrials(1)
                                s_controlHits{mus}{day}(t) = Data{mus}.Events{day}.hits(t);
                            elseif Data{mus}.Events{day}.hits(t) >= laserOff_fullTrials(L) && Data{mus}.Events{day}.hits(t) <= laserOn_fullTrials(L)
                                s_controlHits{mus}{day}(t) = Data{mus}.Events{day}.hits(t);
                            end
                        end
                        if s_controlHits{mus}{day}(t) == 0
                            s_laserHits{mus}{day}(t) = Data{mus}.Events{day}.hits(t);
                        end
                    end
                    s_controlMisses{mus}{day} = NaN(length(Data{mus}.Events{day}.misses),1);
                    s_laserMisses{mus}{day} = NaN(length(Data{mus}.Events{day}.misses),1);
                    for t = 1:length(Data{mus}.Events{day}.misses) % Misses
                        for L = 1:length(laserOff_fullTrials)
                            if Data{mus}.Events{day}.misses(t) <= laserOff_fullTrials(1)
                                s_controlMisses{mus}{day}(t) = Data{mus}.Events{day}.misses(t);
                            elseif Data{mus}.Events{day}.misses(t) >= laserOff_fullTrials(L) && Data{mus}.Events{day}.misses(t) <= laserOn_fullTrials(L)
                                s_controlMisses{mus}{day}(t) = Data{mus}.Events{day}.misses(t);
                            end
                        end
                        if s_controlMisses{mus}{day}(t) == 0
                            s_laserMisses{mus}{day}(t) = Data{mus}.Events{day}.misses(t);
                        end
                    end
                elseif  tMetadata{mus}{day}.trial == "no-go"
                    s_controlHits{mus}{day} = Data{mus}.Events{day}.hits';
                    s_controlMisses{mus}{day} = Data{mus}.Events{day}.misses';
                    s_controlCRs{mus}{day} = NaN(length(Data{mus}.Events{day}.correctRejections),1);
                    s_laserCRs{mus}{day} = NaN(length(Data{mus}.Events{day}.correctRejections),1);
                    for t = 1:length(Data{mus}.Events{day}.correctRejections) % CR
                        for L = 1:length(laserOff_fullTrials)
                            if Data{mus}.Events{day}.correctRejections(t) <= laserOff_fullTrials(1)
                                s_controlCRs{mus}{day}(t) = Data{mus}.Events{day}.correctRejections(t);
                            elseif Data{mus}.Events{day}.correctRejections(t) >= laserOff_fullTrials(L) && Data{mus}.Events{day}.correctRejections(t) <= laserOn_fullTrials(L)
                                s_controlCRs{mus}{day}(t) = Data{mus}.Events{day}.correctRejections(t);
                            end
                        end
                        if s_controlCRs{mus}{day}(t) == 0
                            s_laserCRs{mus}{day}(t) = Data{mus}.Events{day}.correctRejections(t);
                        end
                    end
                    s_controlFAs{mus}{day} = NaN(length(Data{mus}.Events{day}.falseAlarms),1);
                    s_laserFAs{mus}{day} = NaN(length(Data{mus}.Events{day}.falseAlarms),1);
                    for t = 1:length(Data{mus}.Events{day}.falseAlarms) % FA
                        for L = 1:length(laserOff_fullTrials)
                            if Data{mus}.Events{day}.falseAlarms(t) <= laserOff_fullTrials(1)
                                s_controlFAs{mus}{day}(t) = Data{mus}.Events{day}.falseAlarms(t);
                            elseif Data{mus}.Events{day}.falseAlarms(t) >= laserOff_fullTrials(L) && Data{mus}.Events{day}.falseAlarms(t) <= laserOn_fullTrials(L)
                                s_controlFAs{mus}{day}(t) = Data{mus}.Events{day}.falseAlarms(t);
                            end
                        end
                        if s_controlFAs{mus}{day}(t) == 0
                            s_laserFAs{mus}{day}(t) = Data{mus}.Events{day}.falseAlarms(t);
                        end
                    end
                end
                s_controlHits{mus}{day} = s_controlHits{mus}{day}(~isnan(s_controlHits{mus}{day}));
                s_laserHits{mus}{day} = s_laserHits{mus}{day}(~isnan(s_laserHits{mus}{day}));
                s_controlMisses{mus}{day} = s_controlMisses{mus}{day}(~isnan(s_controlMisses{mus}{day}));
                s_laserMisses{mus}{day} = s_laserMisses{mus}{day}(~isnan(s_laserMisses{mus}{day}));
                s_controlCRs{mus}{day} = s_controlCRs{mus}{day}(~isnan(s_controlCRs{mus}{day}));
                s_laserCRs{mus}{day} = s_laserCRs{mus}{day}(~isnan(s_laserCRs{mus}{day}));
                s_controlFAs{mus}{day} = s_controlFAs{mus}{day}(~isnan(s_controlFAs{mus}{day}));
                s_laserFAs{mus}{day} = s_laserFAs{mus}{day}(~isnan(s_laserFAs{mus}{day}));
                % Compute optimal window metrics
                windowNumber = length(0:step:sessionLength-window);
                window_current = 1;
                for w = 0:step:sessionLength-window
                    s_windows_hit_control{mus}{day}(window_current) = length(s_controlHits{mus}{day}(s_controlHits{mus}{day} > w & s_controlHits{mus}{day} < w+window));
                    wincalc_hit_control(window_current) = length(s_controlHits{mus}{day}(s_controlHits{mus}{day} > w & s_controlHits{mus}{day} < w+window)) + 1;
                    s_windows_miss_control{mus}{day}(window_current) = length(s_controlMisses{mus}{day}(s_controlMisses{mus}{day} > w & s_controlMisses{mus}{day} < w+window));
                    wincalc_miss_control(window_current) = length(s_controlMisses{mus}{day}(s_controlMisses{mus}{day} > w & s_controlMisses{mus}{day} < w+window)) + 1;
                    s_windows_CR_control{mus}{day}(window_current) = length(s_controlCRs{mus}{day}(s_controlCRs{mus}{day} > w & s_controlCRs{mus}{day} < w+window));
                    wincalc_CR_control(window_current) = length(s_controlCRs{mus}{day}(s_controlCRs{mus}{day} > w & s_controlCRs{mus}{day} < w+window)) + 1;
                    s_windows_FA_control{mus}{day}(window_current) = length(s_controlFAs{mus}{day}(s_controlFAs{mus}{day} > w & s_controlFAs{mus}{day} < w+window));
                    wincalc_FA_control(window_current) = length(s_controlFAs{mus}{day}(s_controlFAs{mus}{day} > w & s_controlFAs{mus}{day} < w+window)) + 1;
                    s_windows_hit_laser{mus}{day}(window_current) = length(s_laserHits{mus}{day}(s_laserHits{mus}{day} > w & s_laserHits{mus}{day} < w+window));
                    s_windows_miss_laser{mus}{day}(window_current) = length(s_laserMisses{mus}{day}(s_laserMisses{mus}{day} > w & s_laserMisses{mus}{day} < w+window));
                    s_windows_CR_laser{mus}{day}(window_current) = length(s_laserCRs{mus}{day}(s_laserCRs{mus}{day} > w & s_laserCRs{mus}{day} < w+window));
                    s_windows_FA_laser{mus}{day}(window_current) = length(s_laserFAs{mus}{day}(s_laserFAs{mus}{day} > w & s_laserFAs{mus}{day} < w+window));
                    window_current = window_current + 1;
                end
                windows_hitRate_control = s_windows_hit_control{mus}{day} ./ (s_windows_hit_control{mus}{day} + s_windows_miss_control{mus}{day});
                windows_falseAlarmRate_control = s_windows_FA_control{mus}{day} ./ (s_windows_FA_control{mus}{day} + s_windows_CR_control{mus}{day});
                %windows_dPrime = norminv(windows_hitRate) - norminv(windows_falseAlarmRate);
                windows_dPrime_control = (norminv(wincalc_hit_control ./ (wincalc_hit_control + wincalc_miss_control)) - norminv(wincalc_FA_control ./ (wincalc_FA_control + wincalc_CR_control)))';
                optimalWindow = find(windows_dPrime_control > cutoff & windows_hitRate_control > 0.50 & windows_falseAlarmRate_control < 0.50); % control trials only and removing red tail
                optimalWindow_d = cell(1);
                optimalWindow_starts = [];
                optimalWindow_ends =  [];
                for ow = 1:length(optimalWindow)
                    if ow == 1
                        optimalWindow_d{end}(ow) = optimalWindow(ow);
                        optimalWindow_starts(end+1) = (optimalWindow(ow)-1)*step;
                    elseif ow == length(optimalWindow)
                        optimalWindow_d{end}(ow) = optimalWindow(ow);
                        optimalWindow_ends(end+1) = (optimalWindow(ow)-1)*step+600000;
                    else
                        if (optimalWindow(ow) + 1) < optimalWindow(ow+1)
                            optimalWindow_d{end}(ow) = optimalWindow(ow);
                            optimalWindow_ends(end+1) = (optimalWindow(ow)-1)*step+600000;
                            optimalWindow_starts(end+1) = (optimalWindow(ow+1)-1)*step;
                            optimalWindow_d{end + 1} = [];
                        else
                            optimalWindow_d{end}(ow) = optimalWindow(ow);
                        end
                    end
                end
                optimalWindow_d = cellfun(@nonzeros, optimalWindow_d, 'UniformOutput', false);
                if ~isempty(optimalWindow)
                    optimal_hits_control = 0;
                    optimal_misses_control = 0;
                    optimal_CR_control = 0;
                    optimal_FA_control = 0;
                    optimal_hits_laser = 0;
                    optimal_misses_laser = 0;
                    optimal_CR_laser = 0;
                    optimal_FA_laser = 0;
                    opt_hits_control = [];
                    opt_misses_control = [];
                    opt_CR_control = [];
                    opt_FA_control = [];
                    opt_hits_laser = [];
                    opt_misses_laser = [];
                    opt_CR_laser = [];
                    opt_FA_laser = [];
                    for opwin = 1:length(optimalWindow_d)
                        optimal_hits_control = optimal_hits_control + length(find(s_controlHits{mus}{day} > optimalWindow_starts(opwin) & s_controlHits{mus}{day} < optimalWindow_ends(opwin)));
                        opt_hits_control = vertcat(opt_hits_control, s_controlHits{mus}{day}(find(s_controlHits{mus}{day} > optimalWindow_starts(opwin) & s_controlHits{mus}{day} < optimalWindow_ends(opwin))));
                        optimal_misses_control = optimal_misses_control + length(find(s_controlMisses{mus}{day} > optimalWindow_starts(opwin) & s_controlMisses{mus}{day} < optimalWindow_ends(opwin)));
                        opt_misses_control = vertcat(opt_misses_control, s_controlMisses{mus}{day}(find(s_controlMisses{mus}{day} > optimalWindow_starts(opwin) & s_controlMisses{mus}{day} < optimalWindow_ends(opwin))));
                        optimal_CR_control = optimal_CR_control + length(find(s_controlCRs{mus}{day} > optimalWindow_starts(opwin) & s_controlCRs{mus}{day} < optimalWindow_ends(opwin)));
                        opt_CR_control = vertcat(opt_CR_control, s_controlCRs{mus}{day}(find(s_controlCRs{mus}{day} > optimalWindow_starts(opwin) & s_controlCRs{mus}{day} < optimalWindow_ends(opwin))));
                        optimal_FA_control = optimal_FA_control + length(find(s_controlFAs{mus}{day} > optimalWindow_starts(opwin) & s_controlFAs{mus}{day} < optimalWindow_ends(opwin)));
                        opt_FA_control = vertcat(opt_FA_control, s_controlFAs{mus}{day}(find(s_controlFAs{mus}{day} > optimalWindow_starts(opwin) & s_controlFAs{mus}{day} < optimalWindow_ends(opwin))));
                        optimal_hits_laser = optimal_hits_laser + length(find(s_laserHits{mus}{day} > optimalWindow_starts(opwin) & s_laserHits{mus}{day} < optimalWindow_ends(opwin)));
                        opt_hits_laser = vertcat(opt_hits_laser, s_laserHits{mus}{day}(find(s_laserHits{mus}{day} > optimalWindow_starts(opwin) & s_laserHits{mus}{day} < optimalWindow_ends(opwin))));
                        optimal_misses_laser = optimal_misses_laser + length(find(s_laserMisses{mus}{day} > optimalWindow_starts(opwin) & s_laserMisses{mus}{day} < optimalWindow_ends(opwin)));
                        opt_misses_laser = vertcat(opt_misses_laser, s_laserMisses{mus}{day}(find(s_laserMisses{mus}{day} > optimalWindow_starts(opwin) & s_laserMisses{mus}{day} < optimalWindow_ends(opwin))));
                        optimal_CR_laser = optimal_CR_laser + length(find(s_laserCRs{mus}{day} > optimalWindow_starts(opwin) & s_laserCRs{mus}{day} < optimalWindow_ends(opwin)));
                        opt_CR_laser = vertcat(opt_CR_laser, s_laserCRs{mus}{day}(find(s_laserCRs{mus}{day} > optimalWindow_starts(opwin) & s_laserCRs{mus}{day} < optimalWindow_ends(opwin))));
                        optimal_FA_laser = optimal_FA_laser + length(find(s_laserFAs{mus}{day} > optimalWindow_starts(opwin) & s_laserFAs{mus}{day} < optimalWindow_ends(opwin)));
                        opt_FA_laser = vertcat(opt_FA_laser, s_laserFAs{mus}{day}(find(s_laserFAs{mus}{day} > optimalWindow_starts(opwin) & s_laserFAs{mus}{day} < optimalWindow_ends(opwin))));
                    end
                    optimalGO_control = optimal_hits_control + optimal_misses_control;
                    optimalNOGO_control = optimal_CR_control + optimal_FA_control;
                    optimalGO_laser = optimal_hits_laser + optimal_misses_laser;
                    optimalNOGO_laser = optimal_CR_laser + optimal_FA_laser;
                    optimal_length = (((length(optimalWindow)*step)/600000)/3)*4;
                    % Half pushes
                    opt_control_missIncPush = []; % optimal window control misses
                    for miss = 1:length(opt_misses_control)
                        trialIndex = find(Data{mus}.Events{day}.FL_off == opt_misses_control(miss));
                        if isempty(find(sensor_backVector(Data{mus}.Events{day}.FL_off(trialIndex):Data{mus}.Events{day}.BL_on(trialIndex)))) % if the back sensor is OFF during the tone, no push was started
                        else
                            opt_control_missIncPush = [opt_control_missIncPush opt_misses_control(miss)];
                        end
                    end
                    opt_laser_missIncPush = []; % optimal window laser misses
                    for miss = 1:length(opt_misses_laser)
                        trialIndex = find(Data{mus}.Events{day}.FL_off == opt_misses_laser(miss));
                        if isempty(find(sensor_backVector(Data{mus}.Events{day}.FL_off(trialIndex):Data{mus}.Events{day}.BL_on(trialIndex)))) % if the back sensor is OFF during the tone, no push was started
                        else
                            opt_laser_missIncPush = [opt_laser_missIncPush opt_misses_laser(miss)];
                        end
                    end
                    opt_control_FAIncPush = []; % optimal window control false alarms
                    for fa = 1:length(opt_FA_control)
                        trialIndex = find(Data{mus}.Events{day}.FL_off == opt_FA_control(fa));
                        if isempty(find(sensor_backVector(Data{mus}.Events{day}.FL_off(trialIndex):Data{mus}.Events{day}.BL_on(trialIndex)))) % if the back sensor is OFF during the tone, no push was started
                        else
                            opt_control_FAIncPush = [opt_control_FAIncPush opt_FA_control(fa)];
                        end
                    end
                    opt_laser_FAIncPush = []; % optimal window laser false alarms
                    for fa = 1:length(opt_FA_laser)
                        trialIndex = find(Data{mus}.Events{day}.FL_off == opt_FA_laser(fa));
                        if isempty(find(sensor_backVector(Data{mus}.Events{day}.FL_off(trialIndex):Data{mus}.Events{day}.BL_on(trialIndex)))) % if the back sensor is OFF during the tone, no push was started
                        else
                            opt_laser_FAIncPush = [opt_laser_FAIncPush opt_FA_laser(fa)];
                        end
                    end
                    % Calculate proportions
                    s_optGO_control_proportions = 0;
                    s_optGO_laser_proportions = 0;
                    s_optNOGO_control_proportions = 0;
                    s_optNOGO_laser_proportions = 0;
                    if tMetadata{mus}{day}.trial == "go"
                        s_optGO_control_proportions = [(optimal_hits_control/optimalGO_control) ((optimal_misses_control-length(opt_control_missIncPush))/optimalGO_control) (length(opt_control_missIncPush)/optimalGO_control)];
                        s_optGO_laser_proportions = [(optimal_hits_laser/optimalGO_laser) ((optimal_misses_laser-length(opt_laser_missIncPush))/optimalGO_laser) (length(opt_laser_missIncPush)/optimalGO_laser)];
                        s_optGO_proportions{mus}{day} = vertcat(s_optGO_control_proportions,s_optGO_laser_proportions);
                    elseif  tMetadata{mus}{day}.trial == "no-go"
                        s_optNOGO_control_proportions = [(optimal_CR_control/optimalNOGO_control) ((optimal_FA_control-length(opt_control_FAIncPush))/optimalNOGO_control) (length(opt_control_FAIncPush)/optimalNOGO_control)];
                        s_optNOGO_laser_proportions = [(optimal_CR_laser/optimalNOGO_laser) (optimal_FA_laser-length(opt_laser_FAIncPush)/optimalNOGO_laser) (length(opt_laser_FAIncPush)/optimalNOGO_laser)];
                        s_optNOGO_proportions{mus}{day} = vertcat(s_optNOGO_control_proportions,s_optNOGO_laser_proportions);
                    end
                end
            end
        end
    end
end


% figure
% plot()


% Optimal window


% Push duration



%% Plots

% Proportions


% Reaction time


% Push duration



%% Save figures?

% qtn = input('Would you like to save any figures? y/n: ','s');
% if  strcmpi(qtn,'y')
%     dir = input('Please declare a folder to save them in ','s');
%     figs = get(0,'children');
%     for s = 1:length(figs)
%         if length(figs(s).Visible) == 2
%             figure(figs(s))
%             filename = input('If you would like to save this figure write its filename, else press enter ','s');
%             if filename ~= 0
%                 saveas(figs(s),fullfile(dir,filename),'svg')
%             else continue
%             end
%         end
%     end
% end
% 
% close all


