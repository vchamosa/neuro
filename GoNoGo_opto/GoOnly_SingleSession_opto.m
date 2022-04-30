%% Analysis of a single go-only opto session
%
% This script analyses and plots behavioural parameters of an opto go-only session.
% Author: Victor Chamosa Pino, some sections heavily modified 
% from a script by Michelle Sanchez Rivera
%

clear all;
%LOAD THE FILE AND READ IT%
file = '21314001.abf';
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
    signal_shutter     = abf(:,6);
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
    sensorTh   = 4;      % Threshold for sensor state change
    toneTh  = 4;   % Threshold for detecting tone
else
    sensorTh   = 2.5;      % Threshold for sensor state change
    toneTh     = 0.1;      % Threshold for detecting tone
end

sessionLength = length(abf);                       % Length of training session
sessionMins = sessionLength/60000;                 % Length of session in minutes
sensor_frontVector  = (sensor_front > sensorTh);   % Time points at which front sensor is HIGH
sensor_backVector   = (sensor_back > sensorTh);    % Time points at which back sensor is HIGH
signal_rewardVector = (signal_reward > sensorTh);  % Time points at which reward signal is HIGH
signal_lockVector   = (signal_lock > sensorTh);    % Time points at which  the lock signal is HIGH
tone_Vector         = (signal_tone > toneTh);      % Time points at which tone is ON
signal_laserVector  = (signal_laser > sensorTh);   % Time points at which laser is OFF
signal_shutterVector  = (signal_shutter > sensorTh);   % Time points at which shutter is OFF


%% Initial calculations

% Initialise variables
hits                     = 0;
misses                   = 0;
incompletePushesCued     = 0;
uncuedPushes             = 0;
uncuedPushes_complete    = 0; 
uncuedPushes_incomplete  = 0; 

% Define binary change timepoints
toneTimesALL = find(tone_Vector(1:end-1)==0 & tone_Vector(2:end)==1); % finds the times in which the tone signal switches from OFF to ON
toneEndsALL = find(tone_Vector(1:end-1)==1 & tone_Vector(2:end)==0); % finds the times in which the tone signal switches from ON to OFF
laserOff_TimesALL = find(signal_laserVector(1:end-1)==1 & signal_laserVector(2:end)==0);
laserOn_TimesALL = find(signal_laserVector(1:end-1)==0 & signal_laserVector(2:end)==1);
shutterOff_TimesALL = find(signal_shutterVector(1:end-1)==1 & signal_shutterVector(2:end)==0);
shutterOn_TimesALL = find(signal_shutterVector(1:end-1)==0 & signal_shutterVector(2:end)==1);
if ext == '.bin' | ext == '.mat'
    for t = 2:length(toneTimesALL)
        if toneTimesALL(t) - toneTimesALL(t-1) < 80 % defective tones in binary files
            toneTimesALL(t-1) = 0; toneTimesALL(t) = 0;
            toneEndsALL(t-1) = 0; toneEndsALL(t) = 0;
        end
    end
    toneTimesALL = toneTimesALL(find(toneTimesALL));
    toneEndsALL = toneEndsALL(find(toneEndsALL));
end


% Remove unnecessary shutter values
shutterOn_Times = zeros(length(shutterOn_TimesALL),1);
shutterOff_Times = zeros(length(shutterOff_TimesALL),1);
for shut = 1:length(shutterOn_TimesALL)
    if shut == 1
        shutterOn_Times(shut) = shutterOn_TimesALL(shut); %start of first shutter block
    elseif shut == length(shutterOn_TimesALL)
        shutterOff_Times(shut) = shutterOff_TimesALL(shut); %end of last shutter block
    else
        shutter_dif(shut) = shutterOn_TimesALL(shut+1) - shutterOff_TimesALL(shut);
        if shutter_dif(shut) > (2400+(shutterOn_TimesALL(2) - shutterOff_TimesALL(1)))
            shutterOn_Times(shut) = shutterOn_TimesALL(shut+1); %start of shutter block
            shutterOff_Times(shut) = shutterOff_TimesALL(shut); %end of shutter block
        end
    end
end
shutterOn_Times = shutterOn_Times(find(shutterOn_Times));
shutterOff_Times = shutterOff_Times(find(shutterOff_Times));

% Make a trials vector
trialTypesALL = ones(length(toneTimesALL),1);
trialOutcomesALL = zeros(length(trialTypesALL),1); % Clean them up

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
laserOn_earlyPush = [];
laserOff_earlyPush = [];
for lpc = 1: length(laserOn_Times)
    laser_length(lpc) = laserOn_Times(lpc) - laserOff_Times(lpc);
    if laser_length(lpc) > (4*sessionLength/(sessionMins*60))
        laserOn_fullTrials = [laserOn_fullTrials laserOn_Times(lpc)];
        laserOff_fullTrials = [laserOff_fullTrials laserOff_Times(lpc)];
    else
        laserOn_earlyPush = [laserOn_earlyPush laserOn_Times(lpc)];
        laserOff_earlyPush = [laserOff_earlyPush laserOff_Times(lpc)];
    end
end


% Separate control and laser trials (both pre- and post-cue)
controlTrials = zeros(length(toneTimesALL),1);
earlyTrials = zeros(length(shutterOn_Times),1);
for sht = 1:length(shutterOn_Times)
    for tne = 1:length(toneTimesALL)
        if toneTimesALL(tne) >= shutterOn_Times(sht) && toneTimesALL(tne) <= shutterOff_Times(sht)
            earlyTrials(sht) = 1; % full laser trial
        end
        for L = 1:length(laserOff_fullTrials)
            if toneTimesALL(tne) <= laserOff_fullTrials(1)
                trialTypesALL(tne) = 2; % full control trial
            elseif toneTimesALL(tne) >= laserOff_fullTrials(L) && toneTimesALL(tne) <= laserOn_fullTrials(L)
                trialTypesALL(tne) = 2; % full control trial
            end
        end
    end
    for L = 1:length(laserOff_Times)
        if shutterOn_Times(sht) >= laserOff_Times(L) && shutterOff_Times(sht) <= laserOn_Times(L)
            if earlyTrials(sht) == 1
                earlyTrials(sht) = 3; % full control trial
            else
                earlyTrials(sht) = 2; % early control trial
            end
        end
    end
end

%laserTrialtimes = toneTimesALL(~any(controlTrials,2));


% Define trial outcomes
trials_hit   = [];
trials_miss  = [];
control_hit  = [];
control_miss = [];
laser_hit    = [];
laser_miss   = [];
control_early = [];
laser_early   = [];
for n = 1:length(trialOutcomesALL)
    % Define the beginning and end of each trial
    nthTrial_start = toneTimesALL(n);
    if n<length(toneTimesALL)
       nthTrial_end = toneTimesALL(n+1)-1;
    else
        nthTrial_end = length(signal_tone);
    end
    if isempty(find(signal_rewardVector(nthTrial_start:nthTrial_end))) % if the reward signal is OFF during the whole trial
        trialOutcomesALL(n) = 0; %non-rewarded
    else
        trialOutcomesALL(n) = 1; % rewarded
    end
end

for i = 1:length(trialOutcomesALL)
    if trialTypesALL(i) == 2 && trialOutcomesALL(i) == 1
        trials_hit = [trials_hit toneTimesALL(i)];
        control_hit = [control_hit toneTimesALL(i)];
    elseif trialTypesALL(i) == 2 && trialOutcomesALL(i) == 0
        trials_miss = [trials_miss toneTimesALL(i)];
        control_miss = [control_miss toneTimesALL(i)];
    elseif trialTypesALL(i) == 1 && trialOutcomesALL(i) == 1
        trials_hit = [trials_hit toneTimesALL(i)];
        laser_hit = [laser_hit toneTimesALL(i)];
    elseif trialTypesALL(i) == 1 && trialOutcomesALL(i) == 0
        trials_miss = [trials_miss toneTimesALL(i)];
        laser_miss = [laser_miss toneTimesALL(i)];
    end
end
for e = 1:length(shutterOn_Times)
    if earlyTrials(e) == 2
        control_early = [control_early shutterOn_Times(e)];
    elseif earlyTrials(e) == 0
        laser_early = [laser_early shutterOn_Times(e)];
    end
end

controlGO                 = length(find(trialTypesALL == 2));
laserGO                   = length(find(trialTypesALL == 1));
hits                      = length(trials_hit);
hits_control              = length(control_hit);
hits_laser                = length(laser_hit);
misses                    = length(trials_miss);
misses_control            = length(control_miss);
misses_laser              = length(laser_miss);
early_control             = length(control_early);
early_laser               = length(laser_early);


%% Other metrics

% Initialise values
window = floor(sessionLength/8); %length of the window = 5 min%
windowNumber = length(0:window:sessionLength-window);
windows_hit           = zeros(windowNumber,1);
windows_miss          = zeros(windowNumber,1);
windows_hit_control   = zeros(windowNumber,1);
windows_miss_control  = zeros(windowNumber,1);
windows_hit_laser     = zeros(windowNumber,1);
windows_miss_laser    = zeros(windowNumber,1);
windows_early_control = zeros(windowNumber,1);
windows_early_laser   = zeros(windowNumber,1);

% Calculate window trial outcomes
window_current = 1;
for w = 0:window:sessionLength-window
    windows_hit(window_current) = length(trials_hit(trials_hit > w & trials_hit < w+window));
    windows_miss(window_current) = length(trials_miss(trials_miss > w & trials_miss < w+window));
    windows_hit_control(window_current) = length(control_hit(control_hit > w & control_hit < w+window));
    windows_miss_control(window_current) = length(trials_miss(control_miss > w & control_miss < w+window));
    windows_hit_laser(window_current) = length(laser_hit(laser_hit > w & laser_hit < w+window));
    windows_miss_laser(window_current) = length(laser_miss(laser_miss > w & laser_miss < w+window));
    windows_early_control(window_current) = length(control_early(control_early > w & control_early < w+window));
    windows_early_laser(window_current) = length(laser_early(laser_early > w & laser_early < w+window));
    window_current = window_current + 1;
end


% Reaction time
pushStartALL = find(sensor_backVector(1:end-1)==0 & sensor_backVector(2:end)==1);
RT_hit = [];
RT_hit_laser = [];
for h = 1:length(control_hit)%control
    % find the push that happens after the hth hit trial begins
    test = pushStartALL - trials_hit(h);
    pushIndex = find(test>0,1); 
    RT_hit = [RT_hit test(pushIndex)];
end
for h = 1:length(laser_hit)%laser
    % find the push that happens after the hth hit trial begins
    test = pushStartALL - trials_hit(h);
    pushIndex = find(test>0,1); 
    RT_hit_laser = [RT_hit test(pushIndex)];
end

% Push duration


% Uncued pushes
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
control_missIncPush = []; % control misses
for miss = 1:length(control_miss)
    trialIndex = find(toneTimesALL == control_miss(miss));
    if isempty(find(sensor_backVector(toneTimesALL(trialIndex):toneEndsALL(trialIndex)))) % if the back sensor is OFF during the tone, no push was started
    else
        control_missIncPush = [control_missIncPush control_miss(miss)];
    end
end
laser_missIncPush = []; % laser misses
for miss = 1:length(laser_miss)
    trialIndex = find(toneTimesALL == laser_miss(miss));
    if isempty(find(sensor_backVector(toneTimesALL(trialIndex):toneEndsALL(trialIndex)))) % if the back sensor is OFF during the tone, no push was started
    else
        laser_missIncPush = [laser_missIncPush laser_miss(miss)];
    end
end


%% Plots

% Summed proportion of trial outcomes
trial_categories = categorical({'Control','Laser'});
GO_control_proportions = [(hits_control/controlGO) ((misses_control-length(control_missIncPush))/controlGO) (length(control_missIncPush)/controlGO)];
GO_laser_proportions = [(hits_laser/laserGO) ((misses_laser-length(laser_missIncPush))/laserGO) (length(laser_missIncPush)/laserGO)];
GO_proportions = vertcat(GO_control_proportions,GO_laser_proportions);
figure
b = bar(trial_categories,GO_proportions,'stacked');
b(1).FaceColor = [0 0.68 0];
title('Summed trial outcomes','FontSize',14,'FontWeight','bold');

% Median proportion of trial outcomes
GO_control_medprop = [(nanmedian(windows_hit_control./(windows_hit_control+windows_miss_control))) (nanmedian(windows_miss_control./(windows_hit_control+windows_miss_control)))];
GO_laser_medprop = [(nanmedian(windows_hit_laser./(windows_hit_laser+windows_miss_laser))) (nanmedian(windows_miss_laser./(windows_hit_laser+windows_miss_laser)))];
GO_medprop = vertcat(GO_control_medprop,GO_laser_medprop);
figure
b = bar(trial_categories,GO_medprop,'stacked');
b(1).FaceColor = [0 0.68 0];
title('Median trial outcomes','FontSize',14,'FontWeight','bold');

% Summed proportion of early pushes
early_control_proportions = [(hits_control/(early_control+controlGO)) ((misses_control-length(control_missIncPush))/(early_control+controlGO)) (length(control_missIncPush)/(early_control+controlGO)) (early_control/(early_control+controlGO))];
early_laser_proportions = [(hits_laser/(early_laser+laserGO)) ((misses_laser-length(laser_missIncPush))/(early_laser+laserGO)) (length(laser_missIncPush)/(early_laser+laserGO)) (early_laser/(early_laser+laserGO))];
early_proportions = vertcat(early_control_proportions,early_laser_proportions);
figure
b = bar(trial_categories,early_proportions,'stacked');
b(1).FaceColor = [0 0.68 0];
title('Summed early pushes','FontSize',14,'FontWeight','bold');

% Median proportion of early pushes
early_control_medprop = [(nanmedian(windows_hit_control./(windows_hit_control+windows_miss_control+windows_early_control))) (nanmedian(windows_miss_control./(windows_hit_control+windows_miss_control+windows_early_control))) (nanmedian(windows_early_control./(windows_hit_control+windows_miss_control+windows_early_control)))];
early_laser_medprop = [(nanmedian(windows_hit_laser./(windows_hit_laser+windows_miss_laser+windows_early_laser))) (nanmedian(windows_miss_laser./(windows_hit_laser+windows_miss_laser+windows_early_laser))) (nanmedian(windows_early_laser./(windows_hit_laser+windows_miss_laser+windows_early_laser)))];
early_medprop = vertcat(early_control_medprop,early_laser_medprop);
figure
b = bar(trial_categories,early_medprop,'stacked');
b(1).FaceColor = [0 0.68 0];
b(3).FaceColor = [0.4940 0.1840 0.5560];
title('Median early pushes','FontSize',14,'FontWeight','bold');

% In-session laser & control success rates
windows_perHit = windows_hit./(windows_hit + windows_miss);
windows_perHit_control = windows_hit_control./(windows_hit_control + windows_miss_control);
windows_perHit_laser = windows_hit_laser./(windows_hit_laser + windows_miss_laser);
windows_perearlyHit_control = windows_hit_control./(windows_hit_control + windows_miss_control + windows_early_control);
windows_perearlyHit_laser = windows_hit_laser./(windows_hit_laser + windows_miss_laser + windows_early_laser);
windows_perearly_control = windows_early_control./(windows_hit_control + windows_miss_control + windows_early_control);
windows_perearly_laser = windows_early_laser./(windows_hit_laser + windows_miss_laser + windows_early_laser);

% Actual trials
figure
hold on
plot(windows_perHit_control,'-','LineWidth',2,'Marker','o','MarkerSize',7,'Color',[0.8500 0.3250 0.0980])
plot(windows_perHit_laser,'-','LineWidth',2,'Marker','o','MarkerSize',7,'Color',[0 0.4470 0.7410])
hold off
ylim([0 1]);
ylabel("Percentage correct",'FontSize',13,'FontWeight','bold');
xlim([0 windowNumber+1]);
xticklabels({'','0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40'})
xlabel('Time','FontSize',13,'FontWeight','bold');
title('Laser vs control go trials','FontSize',14,'FontWeight','bold');

% Including early pushes
figure
hold on
plot(windows_perearlyHit_control,'-','LineWidth',2,'Marker','o','MarkerSize',7,'Color',[0.8500 0.3250 0.0980])
plot(windows_perearlyHit_laser,'-','LineWidth',2,'Marker','o','MarkerSize',7,'Color',[0 0.4470 0.7410])
plot(windows_perearly_control,'-','LineWidth',2,'Marker','o','MarkerSize',7,'Color',[0.8 0.2 0.38])
plot(windows_perearly_laser,'-','LineWidth',2,'Marker','o','MarkerSize',7,'Color',[0.2 0.2 0.56])
hold off
ylim([0 1]);
ylabel("Percentage correct",'FontSize',13,'FontWeight','bold');
xlim([0 windowNumber+1]);
xticklabels({'','0-5','5-10','10-15','15-20','20-25','25-30','30-35','35-40'})
xlabel('Time','FontSize',13,'FontWeight','bold');
title('Laser vs control go trials','FontSize',14,'FontWeight','bold');



