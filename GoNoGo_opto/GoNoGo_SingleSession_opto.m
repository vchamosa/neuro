%% Analysis of a single go/no-go opto session
%
% This script analyses and plots behavioural parameters of an opto go/no-go session.
% Author: Victor Chamosa Pino, some sections heavily modified 
% from a script by Michelle Sanchez Rivera
%

clear all;
%LOAD THE FILE AND READ IT%
file = '21717000.abf';
imaging = true;

[filepath,name,ext] = fileparts(file);
if ext == '.daq'
    abf = daqread(file);
elseif ext == '.abf'
    abf = abfLoad2(file);
elseif ext == '.bin'
    fid2=fopen(file,'r');
    [abf,count] = fread(fid2,[7,inf],'double');
    fclose(fid2); clear count;
    abf = abf';
elseif ext == '.mat'
    abf = load(file);
    abf = abf.chanData;
end


%% Load behavioural channels

% Load the channels as logical vectors
if ext == '.daq' | ext == '.abf'
    %     sensor_front       = abf(:,3);
    %     sensor_back        = abf(:,2);
    %     signal_reward      = abf(:,4);
    %     signal_lock        = abf(:,5);
    %     signal_tone        = abf(:,6);
    sensor_front       = abf(:,2);
    sensor_back        = abf(:,1);
    signal_reward      = abf(:,3);
    signal_lock        = abf(:,4);
    signal_tone        = abf(:,5);
    signal_laser       = abf(:,7);
elseif ext == '.mat'
    %if imaging
    sensor_front       = abf(:,3); %2
    sensor_back        = abf(:,2); %1
    signal_reward      = abf(:,5); %4
    signal_lock        = abf(:,4); %3
    signal_tone        = abf(:,6); %5
    %else
    sensor_front       = abf(:,2); %2
    sensor_back        = abf(:,1); %1
    signal_reward      = abf(:,4); %4
    signal_lock        = abf(:,3); %3
    signal_tone        = abf(:,5); %5
    %end
end

if ext == '.mat'
    sensorTh = 4;   % Threshold for sensor state change
    toneTh   = 4;   % Threshold for detecting tone
else
    sensorTh   = 2.5;      % Threshold for sensor state change
    toneTh     = 0.1;      % Threshold for detecting tone
end

sessionLength = length(abf);                       %length of training session
sessionMins = sessionLength/60000;                 % Length of session in minutes
sensor_frontVector  = (sensor_front > sensorTh);   % Time points at which front sensor is HIGH
sensor_backVector   = (sensor_back > sensorTh);    % Time points at which back sensor is HIGH
signal_rewardVector = (signal_reward > sensorTh);  % Time points at which reward signal is HIGH
signal_lockVector   = (signal_lock > sensorTh);    % Time points at which the lock signal is HIGH
tone_Vector         = (signal_tone > toneTh);      % Time points at which tone is ON
signal_laserVector  = (signal_laser > sensorTh);   % Time points at which laser is OFF


%% Initial calculations

% Initialise variables
tonesGO                  = 0;
tonesNOGO                = 0;
hits                     = 0;
misses                   = 0;
correctRejections        = 0;
falseAlarms              = 0;
incompletePushesCued     = 0;
uncuedPushes             = 0;
uncuedPushes_complete    = 0;
uncuedPushes_incomplete  = 0;

% Define binary change timepoints
toneTimesALL = find(tone_Vector(1:end-1)==0 & tone_Vector(2:end)==1); % finds the times in which the tone signal switches from OFF to ON
toneEndsALL = find(tone_Vector(1:end-1)==1 & tone_Vector(2:end)==0); % finds the times in which the tone signal switches from ON to OFF
laserOff_TimesALL = find(signal_laserVector(1:end-1)==1 & signal_laserVector(2:end)==0);
laserOn_TimesALL = find(signal_laserVector(1:end-1)==0 & signal_laserVector(2:end)==1);

if ext == '.bin' | ext == '.mat'
    for t = 2:length(toneTimesALL)
        if toneTimesALL(t) - toneTimesALL(t-1) < 80 % defective tones in binary files
            toneTimesALL(t-1) = 0; toneTimesALL(t) = 0;
            toneEndsALL(t-1) = 0; toneEndsALL(t) = 0;
        end
    end
    toneTimesALL     = toneTimesALL(find(toneTimesALL));
    toneEndsALL      = toneEndsALL(find(toneEndsALL));
end

% Remove no-go second 'bump'
for t = 2:length(toneTimesALL)
    if toneTimesALL(t) - toneTimesALL(t-1) < 1700 % if two tones start within one second, it is NoGo tone
        toneTimesALL(t) = 0;% so we remove the "second" tone's start index
        toneEndsALL(t) = 0;
    end
end

% Make a vector that lists all trials in order differentiating them by go or no-go
trialTypesALL = toneTimesALL;
for t = 1 : length(trialTypesALL)-1
    if toneTimesALL(t) ~= 0 && toneTimesALL(t+1) == 0
        trialTypesALL(t) = 2; % no-go tone
    elseif toneTimesALL(t) ~= 0 && toneTimesALL(t+1) ~= 0
        trialTypesALL(t) = 1; % go tone
    end
end

% Clean them up
toneTimesALL = toneTimesALL(find(toneTimesALL)); %so we remove the "second" tone's start index
toneEndsALL = toneEndsALL(find(toneEndsALL));
trialTypesALL = trialTypesALL(find(trialTypesALL));
trialTypesALL(length(trialTypesALL)) = 1; %YES?%
trialOutcomesALL = zeros(length(trialTypesALL),1);
%trialEndsALL = zeros(length(trialTypesALL),1);

% Resolve unmatched/unfinished laser binary changes if they exist
laserOn_Times = laserOn_TimesALL;
laserOff_Times = laserOff_TimesALL;
if length(laserOff_TimesALL) < length(laserOn_TimesALL)
    laserOff_Times = [0; laserOff_TimesALL]; % account for initial control trial(s)
elseif length(laserOn_TimesALL) < length(laserOff_TimesALL)
    laserOff_Times = laserOff_Times(1:(end - 1)); % eliminate unfinished trial at the end
elseif length(laserOn_TimesALL) == length(laserOff_TimesALL)
    if laserOn_Times(1) < laserOff_Times(1)
        laserOff_Times = [0; laserOff_Times]; % account for initial control trial
        laserOff_Times = laserOff_Times(1:(end - 1)); % and for the unfinished trial at the end
    end
end
laserOn_fullTrials = [];
laserOff_fullTrials = [];
for lpc = 1: length(laserOn_Times)
    laser_length(lpc) = laserOn_Times(lpc) - laserOff_Times(lpc);
    if laser_length(lpc) > (4*sessionLength/(sessionMins*60))
        laserOn_fullTrials = [laserOn_fullTrials laserOn_Times(lpc)];
        laserOff_fullTrials = [laserOff_fullTrials laserOff_Times(lpc)];
    end
end

% Separate control and laser trials
controlTrials = zeros(length(toneTimesALL),1);
for t = 1:length(toneTimesALL)
    for L = 1:length(laserOff_fullTrials)
        if toneTimesALL(t) <= laserOff_fullTrials(1)
            controlTrials(t) = 1;
        elseif toneTimesALL(t) >= laserOff_fullTrials(L) && toneTimesALL(t) <= laserOn_fullTrials(L)
            controlTrials(t) = 1;
        end
    end
    if controlTrials(t) == 0 && trialTypesALL(t) == 1
        trialTypesALL(t) = 3; % Laser go trial
    elseif controlTrials(t) == 0 && trialTypesALL(t) == 2
        trialTypesALL(t) = 4; % Laser no-go trial
    end
end
laserTrialtimes = toneTimesALL(~any(controlTrials,2));

% Define trial outcomes
trials_hit   = [];
trials_CR    = [];
trials_miss  = [];
trials_FA    = [];
control_hit  = [];
control_CR   = [];
control_miss = [];
control_FA   = [];
laser_hit    = [];
laser_CR     = [];
laser_miss   = [];
laser_FA     = [];
for n = 1:length(trialOutcomesALL)
    % Define the beginning and end of each trial
    nthTrial_start = toneTimesALL(n);
    if n<length(toneTimesALL)
        nthTrial_end = toneTimesALL(n+1)-1;
    else
        nthTrial_end = length(signal_tone);
    end
    %trialEndsALL(n) = nthTrial_end;
    if isempty(find(signal_rewardVector(nthTrial_start:nthTrial_end))) % if the reward signal is OFF during the whole trial
        trialOutcomesALL(n) = 0; %non-rewarded
    else
        trialOutcomesALL(n) = 1; % rewarded
    end
end

for i = 1:length(trialOutcomesALL)
    if trialTypesALL(i) == 1 && trialOutcomesALL(i) == 1
        trials_hit = [trials_hit toneTimesALL(i)];
        control_hit = [control_hit toneTimesALL(i)];
        trialtest(i) = 5;
    elseif trialTypesALL(i) == 1 && trialOutcomesALL(i) == 0
        trials_miss = [trials_miss toneTimesALL(i)];
        control_miss = [control_miss toneTimesALL(i)];
        trialtest(i) = 7;
    elseif trialTypesALL(i) == 2 && trialOutcomesALL(i) == 1
        trials_CR = [trials_CR toneTimesALL(i)];
        control_CR = [control_CR toneTimesALL(i)];
                trialtest(i) = 1;
    elseif trialTypesALL(i) == 2 && trialOutcomesALL(i) == 0
        trials_FA = [trials_FA toneTimesALL(i)];
        control_FA = [control_FA toneTimesALL(i)];
                trialtest(i) = 3;
    elseif trialTypesALL(i) == 3 && trialOutcomesALL(i) == 1
        trials_hit = [trials_hit toneTimesALL(i)];
        laser_hit = [laser_hit toneTimesALL(i)];
                trialtest(i) = 6;
    elseif trialTypesALL(i) == 3 && trialOutcomesALL(i) == 0
        trials_miss = [trials_miss toneTimesALL(i)];
        laser_miss = [laser_miss toneTimesALL(i)];
                trialtest(i) = 8;
    elseif trialTypesALL(i) == 4 && trialOutcomesALL(i) == 1
        trials_CR = [trials_CR toneTimesALL(i)];
        laser_CR = [laser_CR toneTimesALL(i)];
                trialtest(i) = 2;
    elseif trialTypesALL(i) == 4 && trialOutcomesALL(i) == 0
        trials_FA = [trials_FA toneTimesALL(i)];
        laser_FA = [laser_FA toneTimesALL(i)];
                trialtest(i) = 4;
    end
end
trialtest = trialtest';

tonesGO                   = length(find(trialTypesALL == 1)) + length(find(trialTypesALL == 3));
controlGO                 = length(find(trialTypesALL == 1));
laserGO                   = length(find(trialTypesALL == 3));
tonesNOGO                 = length(find(trialTypesALL == 2)) + length(find(trialTypesALL == 4));
controlNOGO               = length(find(trialTypesALL == 2));
laserNOGO                 = length(find(trialTypesALL == 4));
hits                      = length(trials_hit);
hits_control              = length(control_hit);
hits_laser                = length(laser_hit);
misses                    = length(trials_miss);
misses_control            = length(control_miss);
misses_laser              = length(laser_miss);
correctRejections         = length(trials_CR);
correctRejections_control = length(control_CR);
correctRejections_laser   = length(laser_CR);
falseAlarms               = length(trials_FA);
falseAlarms_control       = length(control_FA);
falseAlarms_laser         = length(laser_FA);

% Discrimination metrics
hitRate = hits/tonesGO; % all trials
falseAlarmRate = falseAlarms/tonesNOGO;
dPrime = norminv((hits+1)/(tonesGO+2)) - norminv((falseAlarms+1)/(tonesNOGO+2));

hitRate_control = hits_control/controlGO; % control trials only
falseAlarmRate_control = falseAlarms_control/controlNOGO;
dPrime_control = norminv((hits_control+1)/(controlGO+2)) - norminv((falseAlarms_control+1)/(controlNOGO+2));


%% Create motivation plot if the behaviour is go/no-go

if tonesNOGO == 0
    disp('This session is not go/no-go');
    pause
end

% Initialise values
window = floor(sessionLength/4); % length of the window = 10 min
step = floor(sessionLength/160000); % step = 15 s
windowNumber = length(0:step:sessionLength-window);
windows_hit          = zeros(windowNumber,1);
windows_miss         = zeros(windowNumber,1);
windows_CR           = zeros(windowNumber,1);
windows_FA           = zeros(windowNumber,1);
windows_hit_control  = zeros(windowNumber,1);
windows_miss_control = zeros(windowNumber,1);
windows_CR_control   = zeros(windowNumber,1);
windows_FA_control   = zeros(windowNumber,1);
windows_hit_laser    = zeros(windowNumber,1);
windows_miss_laser   = zeros(windowNumber,1);
windows_CR_laser     = zeros(windowNumber,1);
windows_FA_laser     = zeros(windowNumber,1);

% Calculate window trial outcomes
window_current = 1;
for w = 0:step:sessionLength-window
    windows_hit(window_current) = length(trials_hit(trials_hit > w & trials_hit < w+window));
    wincalc_hit(window_current) = length(trials_hit(trials_hit > w & trials_hit < w+window)) + 1;
    windows_miss(window_current) = length(trials_miss(trials_miss > w & trials_miss < w+window));
    wincalc_miss(window_current) = length(trials_miss(trials_miss > w & trials_miss < w+window)) + 1;
    windows_CR(window_current) = length(trials_CR(trials_CR > w & trials_CR < w+window));
    wincalc_CR(window_current) = length(trials_CR(trials_CR > w & trials_CR < w+window)) + 1;
    windows_FA(window_current) = length(trials_FA(trials_FA > w & trials_FA < w+window));
    wincalc_FA(window_current) = length(trials_FA(trials_FA > w & trials_FA < w+window)) + 1;
    windows_hit_control(window_current) = length(control_hit(control_hit > w & control_hit < w+window));
    wincalc_hit_control(window_current) = length(control_hit(control_hit > w & control_hit < w+window)) + 1;
    windows_miss_control(window_current) = length(trials_miss(control_miss > w & control_miss < w+window));
    wincalc_miss_control(window_current) = length(trials_miss(control_miss > w & control_miss < w+window)) + 1;
    windows_CR_control(window_current) = length(trials_CR(control_CR > w & control_CR < w+window));
    wincalc_CR_control(window_current) = length(trials_CR(control_CR > w & control_CR < w+window)) + 1;
    windows_FA_control(window_current) = length(trials_FA(control_FA > w & control_FA < w+window));
    wincalc_FA_control(window_current) = length(trials_FA(control_FA > w & control_FA < w+window)) + 1;
    windows_hit_laser(window_current) = length(laser_hit(laser_hit > w & laser_hit < w+window));
    windows_miss_laser(window_current) = length(laser_miss(laser_miss > w & laser_miss < w+window));
    windows_CR_laser(window_current) = length(laser_CR(laser_CR > w & laser_CR < w+window));
    windows_FA_laser(window_current) = length(laser_FA(laser_FA > w & laser_FA < w+window));
    window_current = window_current + 1;
end

% Calculate Hit rate and False Alarm rate for all trials
windows_hitRate = windows_hit ./ (windows_hit + windows_miss); % Log-linear method
windows_falseAlarmRate = windows_FA ./ (windows_FA + windows_CR);
windows_dPrime = (norminv(wincalc_hit ./ (wincalc_hit + wincalc_miss)) - norminv(wincalc_FA ./ (wincalc_FA + wincalc_CR)))';

% Calculate Hit rate and False Alarm rate for control trials
windows_hitRate_control = windows_hit_control ./ (windows_hit_control + windows_miss_control); % Log-linear method
windows_falseAlarmRate_control = windows_FA_control ./ (windows_FA_control + windows_CR_control);
windows_dPrime_control = (norminv(wincalc_hit_control ./ (wincalc_hit_control + wincalc_miss_control)) - norminv(wincalc_FA_control ./ (wincalc_FA_control + wincalc_CR_control)))';

% Calculate Hit rate and False Alarm rate for laser trials
windows_hitRate_laser = windows_hit_laser ./ (windows_hit_laser + windows_miss_laser);
windows_falseAlarmRate_laser = windows_FA_laser ./ (windows_FA_laser + windows_CR_laser);

%Plot all trials
ROC = figure('Renderer', 'painters'); % comment out the parentheses if figure generation is problematic
axis([0, 1, 0, 1]);
patch([0 0.5 0], [0 0.5 0.5], 'blue', 'EdgeAlpha',0, 'FaceAlpha', 0.2);
patch([0 0.5 0.5 0], [0.5 0.5 1 1], 'green', 'EdgeAlpha',0, 'FaceAlpha', 0.2);
patch([0.5 1 0.5], [0.5 1 1], 'red', 'EdgeAlpha',0, 'FaceAlpha', 0.2);
patch([0 1 1], [0 0 1], 'black', 'EdgeAlpha',0, 'FaceAlpha', 0.2);
hold on
%Plots the curve for the desired cutoff d'
cutoff = 1.5;
X = [0.01:0.01:0.99];
Y = [0.01:0.01:0.99];
[X,Y] = meshgrid(X,Y);
dPrimeContour = norminv(Y)-norminv(X); % Y = Hit Rate, X = False Alarm Rate
V1 = [cutoff,cutoff]; % set the desired cutoff value here
V2 = [dPrime, dPrime];
contour(X,Y,dPrimeContour,V1,'LineColor','k','LineStyle','--','LineWidth',1);
contour(X,Y,dPrimeContour,V2,'LineColor','k','LineStyle','-','LineWidth',2);
colour = linspace(0,3,windowNumber); plotsize = 100;
scatter(windows_falseAlarmRate,windows_hitRate,plotsize,colour,'filled','MarkerEdgeColor','black');
colormap(gray);
xlabel('False Alarm rate'), ylabel('Hit Rate');
title('All trials','FontSize',14,'FontWeight','bold');

% Plot control trials
ROC_control = figure('Renderer', 'painters'); % comment out the parentheses if figure generation is problematic
axis([0, 1, 0, 1]);
patch([0 0.5 0], [0 0.5 0.5], 'blue', 'EdgeAlpha',0, 'FaceAlpha', 0.2);
patch([0 0.5 0.5 0], [0.5 0.5 1 1], 'green', 'EdgeAlpha',0, 'FaceAlpha', 0.2);
patch([0.5 1 0.5], [0.5 1 1], 'red', 'EdgeAlpha',0, 'FaceAlpha', 0.2);
patch([0 1 1], [0 0 1], 'black', 'EdgeAlpha',0, 'FaceAlpha', 0.2);
hold on
%Plots the curve for the desired cutoff d'
V2 = [dPrime_control, dPrime_control];
contour(X,Y,dPrimeContour,V1,'LineColor','k','LineStyle','--','LineWidth',1);
contour(X,Y,dPrimeContour,V2,'LineColor','k','LineStyle','-','LineWidth',2);
colour = linspace(0,3,windowNumber); plotsize = 100;
scatter(windows_falseAlarmRate_control,windows_hitRate_control,plotsize,colour,'filled','MarkerEdgeColor','black');
colormap(gray);
xlabel('False Alarm rate'), ylabel('Hit Rate');
title('Control trials','FontSize',14,'FontWeight','bold');


 %% Optimal window metrics
% 
 optimalWindow = find(windows_dPrime_control > cutoff & windows_hitRate_control > 0.50 & windows_falseAlarmRate_control < 0.50); % control trials only and removing red tail
% optimalWindow_allTrials = find(windows_dPrime > cutoff & windows_hitRate > 0.50 & windows_falseAlarmRate < 0.50); % removing red tail
% 
% optimalWindow_d = cell(1);
% optimalWindow_starts = [];
% optimalWindow_ends =  [];
% for ow = 1:length(optimalWindow) % isolate optimal window time blocks
%     if ow == 1
%         optimalWindow_d{end}(ow) = optimalWindow(ow);
%         optimalWindow_starts(end+1) = (optimalWindow(ow)-1)*step;
%     elseif ow == length(optimalWindow)
%         optimalWindow_d{end}(ow) = optimalWindow(ow);
%         if optimalWindow(ow) == windowNumber
%             optimalWindow_ends(end+1) = ((optimalWindow(ow)-1)*step)+window;
%         else
%             optimalWindow_ends(end+1) = (optimalWindow(ow)-1)*step+window;
%         end
%     else
%         if (optimalWindow(ow) + 1) < optimalWindow(ow+1)
%             optimalWindow_d{end}(ow) = optimalWindow(ow);
%             optimalWindow_ends(end+1) = (optimalWindow(ow)-1)*step+window;
%             optimalWindow_starts(end+1) = (optimalWindow(ow+1)-1)*step;
%             optimalWindow_d{end + 1} = [];
%         else
%             optimalWindow_d{end}(ow) = optimalWindow(ow);
%         end
%     end
% end
% optimalWindow_d = cellfun(@nonzeros, optimalWindow_d, 'UniformOutput', false);
% 
% if isempty(optimalWindow)
%     disp('optimal control behaviour not achieved in this session :(')
% else
%     opt_hits_control = [];
%     opt_misses_control = [];
%     opt_CR_control = [];
%     opt_FA_control = [];
%     opt_hits_laser = [];
%     opt_misses_laser = [];
%     opt_CR_laser = [];
%     opt_FA_laser = [];
%     for opwin = 1:length(optimalWindow_d) % isolate optimal trials by type and stim
%         opt_hits_control = [opt_hits_control control_hit(find(control_hit > optimalWindow_starts(opwin) & control_hit < optimalWindow_ends(opwin)))];
%         opt_misses_control = [opt_misses_control control_miss(find(control_miss > optimalWindow_starts(opwin) & control_miss < optimalWindow_ends(opwin)))];
%         opt_CR_control = [opt_CR_control control_CR(find(control_CR > optimalWindow_starts(opwin) & control_CR < optimalWindow_ends(opwin)))];
%         opt_FA_control = [opt_FA_control control_FA(find(control_FA > optimalWindow_starts(opwin) & control_FA < optimalWindow_ends(opwin)))];
%         opt_hits_laser = [opt_hits_laser laser_hit(find(laser_hit > optimalWindow_starts(opwin) & laser_hit < optimalWindow_ends(opwin)))];
%         opt_misses_laser = [opt_misses_laser laser_miss(find(laser_miss > optimalWindow_starts(opwin) & laser_miss < optimalWindow_ends(opwin)))];
%         opt_CR_laser = [opt_CR_laser laser_CR(find(laser_CR > optimalWindow_starts(opwin) & laser_CR < optimalWindow_ends(opwin)))];
%         opt_FA_laser = [opt_FA_laser laser_FA(find(laser_FA > optimalWindow_starts(opwin) & laser_FA < optimalWindow_ends(opwin)))];
%     end
%     opt_hits_control = unique(opt_hits_control); % remove duplicates
%     opt_hits_laser = unique(opt_hits_laser);
%     opt_misses_control = unique(opt_misses_control);
%     opt_misses_laser = unique(opt_misses_laser);
%     opt_CR_control = unique(opt_CR_control);
%     opt_CR_laser = unique(opt_CR_laser);
%     opt_FA_control = unique(opt_FA_control);
%     opt_FA_laser = unique(opt_FA_laser);
%     optimal_CR_control = length(opt_CR_control); % calculate number of trials above in different configurations
%     optimal_FA_control = length(opt_FA_control);
%     optimal_hits_control = length(opt_hits_control);
%     optimal_misses_control = length(opt_misses_control);
%     optimal_CR_laser = length(opt_CR_laser);
%     optimal_FA_laser = length(opt_FA_laser);
%     optimal_hits_laser = length(opt_hits_laser);
%     optimal_misses_laser = length(opt_misses_laser);
%     optimalGO_control = optimal_hits_control + optimal_misses_control;
%     optimalNOGO_control = optimal_CR_control + optimal_FA_control;
%     optimalGO_laser = optimal_hits_laser + optimal_misses_laser;
%     optimalNOGO_laser = optimal_CR_laser + optimal_FA_laser;
%     optimalGO = optimalGO_control + optimalGO_laser;
%     optimalNOGO = optimalNOGO_control + optimalNOGO_laser;
%     optimal_length = (sum(optimalWindow_ends - optimalWindow_starts)/sessionLength)*(sessionLength/600000); % optimal window duration
% end


%% Binned view

% Initialise bin values
bin = floor(sessionLength/8); % length of the window = 5 min
binNumber = length(0:bin:sessionLength-bin);
bins_hit          = zeros(binNumber,1);
bins_miss         = zeros(binNumber,1);
bins_CR           = zeros(binNumber,1);
bins_FA           = zeros(binNumber,1);
bins_hit_control  = zeros(binNumber,1);
bins_miss_control = zeros(binNumber,1);
bins_CR_control   = zeros(binNumber,1);
bins_FA_control   = zeros(binNumber,1);
bins_hit_laser    = zeros(binNumber,1);
bins_miss_laser   = zeros(binNumber,1);
bins_CR_laser     = zeros(binNumber,1);
bins_FA_laser     = zeros(binNumber,1);

% Calculate window trial outcomes
bin_current = 1;
for w = 0:bin:sessionLength-bin
    bins_hit(bin_current) = length(trials_hit(trials_hit > w & trials_hit < w+bin));
    bins_miss(bin_current) = length(trials_miss(trials_miss > w & trials_miss < w+bin));
    bins_CR(bin_current) = length(trials_CR(trials_CR > w & trials_CR < w+bin));
    bins_FA(bin_current) = length(trials_FA(trials_FA > w & trials_FA < w+bin));
    bins_hit_control(bin_current) = length(control_hit(control_hit > w & control_hit < w+bin));
    bincalc_hit_control(bin_current) = length(control_hit(control_hit > w & control_hit < w+bin)) + 1;
    bins_miss_control(bin_current) = length(trials_miss(control_miss > w & control_miss < w+bin));
    bincalc_miss_control(bin_current) = length(trials_miss(control_miss > w & control_miss < w+bin)) + 1;
    bins_CR_control(bin_current) = length(trials_CR(control_CR > w & control_CR < w+bin));
    bincalc_CR_control(bin_current) = length(trials_CR(control_CR > w & control_CR < w+bin)) + 1;
    bins_FA_control(bin_current) = length(trials_FA(control_FA > w & control_FA < w+bin));
    bincalc_FA_control(bin_current) = length(trials_FA(control_FA > w & control_FA < w+bin)) + 1;
    bins_hit_laser(bin_current) = length(laser_hit(laser_hit > w & laser_hit < w+bin));
    bins_miss_laser(bin_current) = length(laser_miss(laser_miss > w & laser_miss < w+bin));
    bins_CR_laser(bin_current) = length(laser_CR(laser_CR > w & laser_CR < w+bin));
    bins_FA_laser(bin_current) = length(laser_FA(laser_FA > w & laser_FA < w+bin));
    bin_current = bin_current + 1;
end

% In-session laser & control success rates
bins_perHit_laser = bins_hit_laser./(bins_hit_laser + bins_miss_laser);
bins_perCR_laser = bins_CR_laser./(bins_CR_laser + bins_FA_laser);
bins_perHit_control = bins_hit_control./(bins_hit_control + bins_miss_control);
bins_perCR_control = bins_CR_control./(bins_CR_control + bins_FA_control);
bins_perFA_control = bins_FA_control./(bins_CR_control + bins_FA_control);
bins_dPrime_control = (norminv(bincalc_hit_control ./ (bincalc_hit_control + bincalc_miss_control)) - norminv(bincalc_FA_control ./ (bincalc_FA_control + bincalc_CR_control)))';

bin_optimalWindow = find(bins_dPrime_control > cutoff & bins_perHit_control > 0.50 & bins_perFA_control < 0.50); % control trials only and removing red tail

optimalWindow_win = zeros(binNumber,1);

% Optimal window
if isempty(optimalWindow)
    disp('optimal control behaviour not achieved in this session :(')
else
    optB_hits_control = [];
    optB_misses_control = [];
    optB_CR_control = [];
    optB_FA_control = [];
    optB_hits_laser = [];
    optB_misses_laser = [];
    optB_CR_laser = [];
    optB_FA_laser = [];
    Bopt_hits_control = [];
    Bopt_misses_control = [];
    Bopt_CR_control = [];
    Bopt_FA_control = [];
    Bopt_hits_laser = [];
    Bopt_misses_laser = [];
    Bopt_CR_laser = [];
    Bopt_FA_laser = [];
    BoptimalWindow_starts = [];
    BoptimalWindow_ends =  [];
    bOW_blocks = 1;
    for ow = 1:length(bin_optimalWindow)
        optimalWindow_win(bin_optimalWindow(ow)) = 1;
        optB_hits_control(ow) = bins_hit_control(bin_optimalWindow(ow));
        optB_misses_control(ow) =  bins_miss_control(bin_optimalWindow(ow));
        optB_CR_control(ow) = bins_CR_control(bin_optimalWindow(ow));
        optB_FA_control(ow) = bins_FA_control(bin_optimalWindow(ow));
        optB_hits_laser(ow) = bins_hit_laser(bin_optimalWindow(ow));
        optB_misses_laser(ow) = bins_miss_laser(bin_optimalWindow(ow));
        optB_CR_laser(ow) = bins_CR_laser(bin_optimalWindow(ow));
        optB_FA_laser(ow) = bins_FA_laser(bin_optimalWindow(ow));
    end
    optimalB_GO_control = optB_hits_control + optB_misses_control;
    optimalB_NOGO_control = optB_CR_control + optB_FA_control;
    optimalB_GO_laser = optB_hits_laser + optB_misses_laser;
    optimalB_NOGO_laser = optB_CR_laser + optB_FA_laser;
    optimalB_GO = optimalB_GO_control + optimalB_GO_laser;
    optimalB_NOGO = optimalB_NOGO_control + optimalB_NOGO_laser;
    for ow = 1:length(bin_optimalWindow) % isolate optimal window time blocks
        if ow == 1
            BoptimalWindow_starts(end+1) = (bin_optimalWindow(ow)-1)*bin;
        elseif bin_optimalWindow(ow) > (bin_optimalWindow(ow-1) + 1)
            BoptimalWindow_ends(end+1) = bin_optimalWindow(ow-1)*bin;
            BoptimalWindow_starts(end+1) = (bin_optimalWindow(ow)-1)*bin;
            bOW_blocks = bOW_blocks + 1;
        end
        if ow == length(bin_optimalWindow)
            BoptimalWindow_ends(end+1) = bin_optimalWindow(ow)*bin;
        end
    end
    for ow = 1:bOW_blocks
        Bopt_hits_control = [Bopt_hits_control control_hit(find(control_hit > BoptimalWindow_starts(ow) & control_hit < BoptimalWindow_ends(ow)))];
        Bopt_misses_control = [Bopt_misses_control control_miss(find(control_miss > BoptimalWindow_starts(ow) & control_miss < BoptimalWindow_ends(ow)))];
        Bopt_CR_control = [Bopt_CR_control control_CR(find(control_CR > BoptimalWindow_starts(ow) & control_CR < BoptimalWindow_ends(ow)))];
        Bopt_FA_control = [Bopt_FA_control control_FA(find(control_FA > BoptimalWindow_starts(ow) & control_FA < BoptimalWindow_ends(ow)))];
        Bopt_hits_laser = [Bopt_hits_laser laser_hit(find(laser_hit > BoptimalWindow_starts(ow) & laser_hit < BoptimalWindow_ends(ow)))];
        Bopt_misses_laser = [Bopt_misses_laser laser_miss(find(laser_miss > BoptimalWindow_starts(ow) & laser_miss < BoptimalWindow_ends(ow)))];
        Bopt_CR_laser = [Bopt_CR_laser laser_CR(find(laser_CR > BoptimalWindow_starts(ow) & laser_CR < BoptimalWindow_ends(ow)))];
        Bopt_FA_laser = [Bopt_FA_laser laser_FA(find(laser_FA > BoptimalWindow_starts(ow) & laser_FA < BoptimalWindow_ends(ow)))];
    end
end


%% Other metrics

% Uncued pushes UNFINISHED
leverLockALL = find(signal_lockVector(1:end-1)==0 & signal_lockVector(2:end)==1);% finds the times in which the lever locks
toneReset = []; %a list of lever locks that occur after a tone (as part of a trial)
for r = 1:length(toneTimesALL)
    % find the first reset that happens after the rth tone
    test = leverLockALL -  toneTimesALL(r);
    resetIndex = find(test>0,1);
    toneReset = [toneReset leverLockALL(resetIndex)];
end
uncuedPushes = setdiff(leverLockALL,toneReset);
% for igm = 1:length(uncuedPushes) % separate between full and partial pushes
%     trialIndex = find(toneTimesALL == uncuedPushes(igm));
%     if isempty(find(sensor_backVector(toneTimesALL(trialIndex):toneEndsALL(trialIndex)))) % if the back sensor is OFF during the tone, no push was started
%     else
%        trials_missIncPush = [trials_missIncPush uncuedPushes(igm)];
%     end
% end


% Differentiate between full and incomplete pushes
if ~isempty(optimalWindow)
    optB_control_missIncPush = []; % optimal window control misses
    for miss = 1:length(Bopt_misses_control)
        trialIndex = find(toneTimesALL == Bopt_misses_control(miss));
        if isempty(find(sensor_backVector(toneTimesALL(trialIndex):toneEndsALL(trialIndex)))) % if the back sensor is OFF during the tone, no push was started
        else
            optB_control_missIncPush = [optB_control_missIncPush Bopt_misses_control(miss)];
        end
    end
    optB_laser_missIncPush = []; % optimal window laser misses
    for miss = 1:length(Bopt_misses_laser)
        trialIndex = find(toneTimesALL == Bopt_misses_laser(miss));
        if isempty(find(sensor_backVector(toneTimesALL(trialIndex):toneEndsALL(trialIndex)))) % if the back sensor is OFF during the tone, no push was started
        else
            optB_laser_missIncPush = [optB_laser_missIncPush Bopt_misses_laser(miss)];
        end
    end
    optB_control_FAIncPush = []; % optimal window control false alarms
    for fa = 1:length(Bopt_FA_control)
        trialIndex = find(toneTimesALL == Bopt_FA_control(fa));
        if isempty(find(sensor_backVector(toneTimesALL(trialIndex):toneEndsALL(trialIndex)))) % if the back sensor is OFF during the tone, no push was started
        else
            optB_control_FAIncPush = [optB_control_FAIncPush Bopt_FA_control(fa)];
        end
    end
    optB_laser_FAIncPush = []; % optimal window laser false alarms
    for fa = 1:length(Bopt_FA_laser)
        trialIndex = find(toneTimesALL == Bopt_FA_laser(fa));
        if isempty(find(sensor_backVector(toneTimesALL(trialIndex):toneEndsALL(trialIndex)))) % if the back sensor is OFF during the tone, no push was started
        else
            optB_laser_FAIncPush = [optB_laser_FAIncPush Bopt_FA_laser(fa)];
        end
    end
    
    
% Reaction time & push duration
    pushStartALL = find(sensor_backVector(1:end-1)==0 & sensor_backVector(2:end)==1);
    pushEndALL = find(sensor_frontVector(1:end-1)==0 & sensor_frontVector(2:end)==1);
    optRT_hit = [];
    optRT_hit_laser = [];
    opt_pushdur_hit = [];
    opt_pushdur_hit_laser = [];
    for h = 1:length(Bopt_hits_control)%control
        %find the push that happens after the hth hit trial begins
        start_diff = pushStartALL - Bopt_hits_control(h);
        end_diff = pushEndALL - Bopt_hits_control(h);
        start_pushIndex = find(start_diff>0,1);
        end_pushIndex = find(end_diff>0,1);
        optRT_hit = [optRT_hit start_diff(start_pushIndex)/10000];
        opt_pushdur_hit = [opt_pushdur_hit ((end_diff(end_pushIndex) - start_diff(start_pushIndex))/10000)];
    end
    for h = 1:length(Bopt_hits_laser)%laser
        %find the push that happens after the hth hit trial begins
        start_diff = pushStartALL - Bopt_hits_laser(h);
        end_diff = pushEndALL - Bopt_hits_laser(h);
        start_pushIndex = find(start_diff>0,1);
        end_pushIndex = find(end_diff>0,1);
        optRT_hit_laser = [optRT_hit_laser start_diff(start_pushIndex)/10000];
        opt_pushdur_hit_laser = [opt_pushdur_hit_laser ((end_diff(end_pushIndex) - start_diff(start_pushIndex))/10000)];
    end
    
    optRT_FA = [];
    optRT_FA_laser = [];
    opt_pushdur_FA = [];
    opt_pushdur_FA_laser = [];
    for fa = 1:length(Bopt_FA_control)%control
        %find the push that happens after the fath FA trial begins
        start_diff = pushStartALL - Bopt_FA_control(fa);
        end_diff = pushEndALL - Bopt_FA_control(fa);
        start_pushIndex = find(start_diff>0,1);
        end_pushIndex = find(end_diff>0,1);
        optRT_FA = [optRT_FA start_diff(start_pushIndex)/10000];
        opt_pushdur_FA = [opt_pushdur_FA ((end_diff(end_pushIndex) - start_diff(start_pushIndex))/10000)];
    end
    for fa = 1:length(Bopt_FA_laser)%laser
        %find the push that happens after the fath FA trial begins
        start_diff = pushStartALL - Bopt_FA_laser(fa);
        end_diff = pushEndALL - Bopt_FA_laser(fa);
        start_pushIndex = find(start_diff>0,1);
        end_pushIndex = find(end_diff>0,1);
        optRT_FA_laser = [optRT_FA_laser start_diff(start_pushIndex)/10000];
        opt_pushdur_FA_laser = [opt_pushdur_FA_laser ((end_diff(end_pushIndex) - start_diff(start_pushIndex))/10000)];
    end
    
    optRT_missinc = [];
    optRT_missinc_laser = [];
    for h = 1:length(optB_control_missIncPush)%control
        %find the push that happens after the hth hit trial begins
        test = pushStartALL - optB_control_missIncPush(h);
        pushIndex = find(test>0,1);
        optRT_missinc = [optRT_missinc test(pushIndex)/10000];
    end
    for h = 1:length(optB_laser_missIncPush)%laser
        %find the push that happens after the hth hit trial begins
        test = pushStartALL - optB_laser_missIncPush(h);
        pushIndex = find(test>0,1);
        optRT_missinc_laser = [optRT_missinc_laser test(pushIndex)/10000];
    end
    
    optRT_FAinc = [];
    optRT_FAinc_laser = [];
    for fa = 1:length(optB_control_FAIncPush)%control
        %find the push that happens after the fath FA trial begins
        test = pushStartALL - optB_control_FAIncPush(fa);
        pushIndex = find(test>0,1);
        optRT_FAinc = [optRT_FAinc test(pushIndex)/10000];
    end
    for fa = 1:length(optB_laser_FAIncPush)%laser
        %find the push that happens after the fath FA trial begins
        test = pushStartALL - optB_laser_FAIncPush(fa);
        pushIndex = find(test>0,1);
        optRT_FAinc_laser = [optRT_FAinc_laser test(pushIndex)/10000];
    end
end


%% Plots

% Binned session trajectories
figure % total
hold on
yyaxis right
bar(optimalWindow_win,1,'FaceColor','green','FaceAlpha',0.4,'EdgeColor','none');
plot(bins_perCR_control,'-','LineWidth',2,'Marker','o','MarkerSize',7,'Color','b')
plot(bins_perHit_control,'-','LineWidth',2,'Marker','o','MarkerSize',7,'Color','r')
title('Control trials','FontSize',14,'FontWeight','bold');
ylim([0 1]);
ylabel("Percentage correct",'FontSize',13,'FontWeight','bold');
yyaxis left
plot(bins_dPrime_control,'-','LineWidth',2,'Marker','o','MarkerSize',7,'Color','k')
ylim([-1.5 3.5]);
ylabel("d'",'FontSize',13,'FontWeight','bold');
hold off
xlim([0 binNumber+1]);
xticklabels({'','0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40'})
xlabel('Time','FontSize',13,'FontWeight','bold');

figure % control
hold on
bar(optimalWindow_win,1,'FaceColor','green','FaceAlpha',0.4,'EdgeColor','none');
if laserGO < 10
    plot(bins_perCR_control,'-','LineWidth',2,'Marker','o','MarkerSize',7,'Color',[0.8500 0.3250 0.0980])
    plot(bins_perCR_laser,'-','LineWidth',2,'Marker','o','MarkerSize',7,'Color',[0 0.4470 0.7410])
    title('Laser vs control no-go trials','FontSize',14,'FontWeight','bold');
elseif laserNOGO < 3
    plot(bins_perHit_control,'-','LineWidth',2,'Marker','o','MarkerSize',7,'Color',[0.8500 0.3250 0.0980])
    plot(bins_perHit_laser,'-','LineWidth',2,'Marker','o','MarkerSize',7,'Color',[0 0.4470 0.7410])
    title('Laser vs control go trials','FontSize',14,'FontWeight','bold');
end
hold off
ylim([0 1]);
ylabel("Percentage correct",'FontSize',13,'FontWeight','bold');
xlim([0 binNumber+1]);
xticklabels({'','0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40'})
xlabel('Time','FontSize',13,'FontWeight','bold');


% Proportions (sum)
trial_categories = categorical({'Control','Laser'});
if ~isempty(bin_optimalWindow) % optimal window trials only
    if laserGO < 10
        optNOGO_control_sumprop = [(sum(optB_CR_control)/sum(optimalB_NOGO_control)) ((sum(optB_FA_control) - length(optB_control_FAIncPush))/sum(optimalB_NOGO_control)) (length(optB_control_FAIncPush)/sum(optimalB_NOGO_control))];
        optNOGO_laser_sumprop = [(sum(optB_CR_laser)/sum(optimalB_NOGO_laser)) ((sum(optB_FA_laser) - length(optB_laser_FAIncPush))/sum(optimalB_NOGO_laser)) (length(optB_laser_FAIncPush)/sum(optimalB_NOGO_laser))];
        optNOGO_sumprop = vertcat(optNOGO_control_sumprop,optNOGO_laser_sumprop);
        figure
        b = bar(trial_categories,optNOGO_sumprop,'stacked');
        b(1).FaceColor = [0 0.68 0];
        ylim([0 1]);
        title('Optimal no-go summed outcome proportions','FontSize',14,'FontWeight','bold');
    elseif laserNOGO < 3
        optGO_control_sumprop = [(sum(optB_hits_control)/sum(optimalB_GO_control)) ((sum(optB_misses_control) - length(optB_control_missIncPush))/sum(optimalB_GO_control)) (length(optB_control_missIncPush)/sum(optimalB_GO_control))];
        optGO_laser_sumprop = [(sum(optB_hits_laser)/sum(optimalB_GO_laser)) ((sum(optB_misses_laser) - length(optB_laser_missIncPush))/sum(optimalB_GO_laser)) (length(optB_laser_missIncPush)/sum(optimalB_GO_laser))];
        optGO_sumprop = vertcat(optGO_control_sumprop,optGO_laser_sumprop);
        figure
        b = bar(trial_categories,optGO_sumprop,'stacked');
        b(1).FaceColor = [0 0.68 0];
        ylim([0 1]);
        title('Optimal go summed outcome proportions','FontSize',14,'FontWeight','bold');
    end
end


% Proportions (median)
if ~isempty(bin_optimalWindow) % optimal window trials only
    if laserGO < 10
        optNOGO_control_medprop = [(nanmedian(optB_CR_control./optimalB_NOGO_control)) (nanmedian(optB_FA_control./optimalB_NOGO_control))];
        optNOGO_laser_medprop = [(nanmedian(optB_CR_laser./optimalB_NOGO_laser)) (nanmedian(optB_FA_laser./optimalB_NOGO_laser))];
        optNOGO_medprop = vertcat(optNOGO_control_medprop,optNOGO_laser_medprop);
        figure
        b = bar(trial_categories,optNOGO_medprop,'stacked');
        b(1).FaceColor = [0 0.68 0];
        ylim([0 1]);
        title('Optimal no-go median outcome proportions','FontSize',14,'FontWeight','bold');
    elseif laserNOGO < 3
        optGO_control_medprop = [(nanmedian(optB_hits_control./optimalB_GO_control)) (nanmedian(optB_misses_control./optimalB_GO_control))];
        optGO_laser_medprop = [(nanmedian(optB_hits_laser./optimalB_GO_laser)) (nanmedian(optB_misses_laser./optimalB_GO_laser))];
        optGO_medprop = vertcat(optGO_control_medprop,optGO_laser_medprop);
        figure
        b = bar(trial_categories,optGO_medprop,'stacked');
        b(1).FaceColor = [0 0.68 0];
        ylim([0 1]);
        title('Optimal go median outcome proportions','FontSize',14,'FontWeight','bold');
    end
end


% Data saves
%RT_21501002 = {optRT_FA  optRT_FA_laser  optRT_hit};
%save('D:\\ALM_Ardbeg\\Behavioural_data\\reaction_times.mat','RT_21501002','-append','-nocompression')

%prop_21501002 = {bins_perCR_laser bins_perCR_control};
%save('D:\\ALM_Ardbeg\\Behavioural_data\\proportions.mat','prop_21501002','-append','-nocompression')

%push_21501002 = {opt_pushdur_FA opt_pushdur_FA_laser opt_pushdur_hit};
%save('D:\\ALM_Ardbeg\\Behavioural_data\\push_durations.mat','push_21501002','-append','-nocompression')

% Event times files for video analysis
%save('D:\\ALM_Ardbeg\\Behavioural_data\\eventTimes_21705003.mat','Bopt_hits_control','Bopt_misses_control','Bopt_CR_control','Bopt_FA_control','Bopt_hits_laser','Bopt_misses_laser','Bopt_CR_laser','Bopt_FA_laser','uncuedPushes','sessionLength');
%save('D:\\ALM_Ardbeg\\Behavioural_data\\behavMetrics_21705003.mat','optRT_hit','optRT_hit_laser','opt_pushdur_hit','opt_pushdur_hit_laser','optRT_FA','optRT_FA_laser','opt_pushdur_FA','opt_pushdur_FA_laser');


