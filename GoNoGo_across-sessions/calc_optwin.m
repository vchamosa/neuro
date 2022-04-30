function [opt_rew,opt_time,opt_percor,learned,learninglength] = calc_optwin(Data,metadata,varargin)
%% Calculates optimal window parameters across sessions of go/no-go behaviour
%
%  INPUTS
% Must include Data and metadata in order to compute optimal window stats.
% If no further specification is given it will do so to all mice provided
% by both inputs. If further specification is required, include name/value 
% pairs of categories as found in metadata and animals.json. If filters are
% used in combination the anaylsis will involve mice that satisfy both
% criteria.
%
% For example, calc_optwin(Data,metadata,'responsible_person','Isolde', 
% 'strain','MCos') will plot the correct rejections of MCos mice trained by 
% Isolde.
% Author: Victor Chamosa Pino
%


% Stop if not enough inputs supplied
if round(nargin/2) ~= nargin/2
    error('Please ensure you are providing sufficient input');
    return
end


% Parse optional inputs
if nargin > 3
    % Find animals that meet filter criteria
    input = reshape(varargin,2,[]);
    miceNames = cellfun(@(x) getfield(x,'mouseID'),metadata,'UniformOutput',false);
    selectedMiceLog = zeros(size(miceNames));
    for mouse = 1:numel(miceNames)
        selected = 1;
        for pair = input
            name = pair{1};
            value = pair{2};
            if strcmpi(name,'strain')
                if contains(value,'thy','IgnoreCase',true)
                    value = 'thy';
                elseif contains(value,'vgat','IgnoreCase',true)
                    value = 'vga';
                end
            end
            % if any name does not have desired value, deselect mouse
            if ~contains(metadata{mouse}.(name),value,'IgnoreCase',true)
                selected = 0;
            end
        end
        % Generate values for analysis
        selectedMiceLog(mouse) = selected;
        micedata = Data(logical(selectedMiceLog));
        numMice = numel(micedata);
    end
    miceNames = miceNames(logical(selectedMiceLog));
%      tMetadata = cell(1, numMice);
%      for iMouse = 1 : numMice
%          tMetadata{iMouse} = loadjson([miceNames{iMouse} '.json']);
%      end
    
else
    % Quick adjustments if no further specification has been made
    numMice = numel(metadata);
    micedata = Data;
    tMetadata = cell(1, numMice);
    miceNames = cell(1, numMice);
    for iMouse = 1 : numMice
        miceNames{iMouse} = metadata{1,iMouse}.mouseID;
%         tMetadata{iMouse} = loadjson([miceNames{iMouse} '.json']);
    end
end


%% logicals for random mode - only necessary if analysing alternating task structure, which has been deprecated
% randomSes = cell(1, numMice);
% for iMouse = 1 : numMice
%     numDays = numel(micedata{iMouse}.Events);
%     randoms = zeros(1, numDays);
%     for day = 1 : numDays
%         if isfield(tMetadata{iMouse}{day},'mode') && strcmp(tMetadata{iMouse}{day}.mode, 'random')
%             randoms(day) = 1;
%         end
%     end
%     randomSes{iMouse} = randoms;
% end 


%% Optimal window calculations and plots

cutoff = 1.5;
window = 600000; %length of the window = 10 min
step = 15000; %step = 15s
opt_wnd = cell(numMice,60);
opt_time = NaN(numMice,60);
optnorm_time = NaN(numMice,10);
opt_hits = zeros(numMice,60);
opt_CRs = zeros(numMice,60);
opt_FAs = zeros(numMice,60);
opt_misses = zeros(numMice,60);
for i = 1:numMice
    for day = 1:numel(micedata{i}.Events)
        if isstring(micedata{i}.Events{day}) %|| randomSes{i}(day) == 0 % add for random sessions only
            opt_wnd{i,day} = NaN;
            %opt_time(i,day) = NaN;
        else
            sessionLength = micedata{i}.Events{day}.dur;
            windowNumber = length(0:step:sessionLength-window);
            windows_hit  = zeros(windowNumber,1);
            windows_miss = zeros(windowNumber,1);
            windows_CR   = zeros(windowNumber,1);
            windows_FA   = zeros(windowNumber,1);

            window_current = 1;
            for w = 0:step:sessionLength-window
                windows_hit(window_current) = numel(micedata{i}.Events{day}.hits(micedata{i}.Events{day}.hits > w & micedata{i}.Events{day}.hits < w+window));
                windows_miss(window_current) = numel(micedata{i}.Events{day}.misses(micedata{i}.Events{day}.misses > w & micedata{i}.Events{day}.misses < w+window));
                windows_CR(window_current) = numel(micedata{i}.Events{day}.correctRejections(micedata{i}.Events{day}.correctRejections > w & micedata{i}.Events{day}.correctRejections < w+window));
                windows_FA(window_current) = numel(micedata{i}.Events{day}.falseAlarms(micedata{i}.Events{day}.falseAlarms > w & micedata{i}.Events{day}.falseAlarms < w+window));
                window_current = window_current + 1;
            end

            windows_hitRate = windows_hit ./ (windows_hit + windows_miss);
            windows_falseAlarmRate = windows_FA ./ (windows_FA + windows_CR);
            windows_hitRate(windows_hitRate == 0) = 0.0001; windows_hitRate(windows_hitRate == 1)  = 0.999;
            windows_falseAlarmRate(windows_falseAlarmRate == 0) = 0.0001; windows_falseAlarmRate(windows_falseAlarmRate == 1) = 0.999;
            windows_dPrime = norminv(windows_hitRate) - norminv(windows_falseAlarmRate);

            optimalWindow = find(windows_dPrime > cutoff & windows_hitRate > 0.5);% & windows_falseAlarmRate < 0.5);
            if ~isempty(optimalWindow)
                optimalWindow_start = (optimalWindow(1)-1)*step;
                optimalWindow_end = (optimalWindow(end)-1)*step+window;
                opt_wnd{i,day} = [optimalWindow_start optimalWindow_end];
                opt_time(i,day) = ((optimalWindow_end-optimalWindow_start)/micedata{i}.Events{day}.sampleRate)/60;
                opt_hits(i,day) = numel(find(micedata{i}.Events{day}.hits > optimalWindow_start & micedata{i}.Events{day}.hits < optimalWindow_end));
                opt_misses(i,day) = numel(find(micedata{i}.Events{day}.misses > optimalWindow_start & micedata{i}.Events{day}.misses < optimalWindow_end));
                opt_CRs(i,day) = numel(find(micedata{i}.Events{day}.correctRejections > optimalWindow_start & micedata{i}.Events{day}.correctRejections < optimalWindow_end));
                opt_FAs(i,day) = numel(find(micedata{i}.Events{day}.falseAlarms > optimalWindow_start & micedata{i}.Events{day}.falseAlarms < optimalWindow_end));
            else
                opt_time(i,day) = 0;
            end
        end
    end
    
    % For normalised training length
    sesh_norm(i) = round(numel(micedata{i}.Events)/10); % in 10% jumps
    for sn = 1:10
        session_unit = sesh_norm(i)*sn;
        if session_unit > numel(micedata{i}.Events)
            session_unit = numel(micedata{i}.Events);
        end
        optnorm_time(i,sn) = opt_time(i,session_unit);
    end
end
opt_rew = opt_hits + opt_CRs;
opt_tones = opt_rew + opt_FAs + opt_misses;
opt_percor = (opt_rew./opt_tones)*100;

% decrease in optimal window-free sessions
% for day = 1:50
%     non_opt(day) = length(find(opt_time(:,day) == 0));
%     norm_opt(day) = length(find(opt_time(:,day) >= 20))/length(nonzeros(~isnan(opt_time(:,day))));
% end
% 
% figure
% bar(non_opt)
% xlim([0 50])
% ylim([0 numMice])
% ylabel('Number of mice','FontSize',12,'FontWeight','bold')
% xlabel('Session','FontSize',12,'FontWeight','bold')
% title('Non-optimal sessions','FontSize',14,'FontWeight','bold')
% 
% figure
% bar(norm_opt)
% xlim([0 50])
% ylim([0 1])
% ylabel('Normalised no. of mice','FontSize',12,'FontWeight','bold')
% xlabel('Session','FontSize',12,'FontWeight','bold')
% title('Optimal sessions','FontSize',14,'FontWeight','bold')

% remove values of 0 to plot OW length when there is an optimal window
opt_rew(opt_rew == 0) = NaN;
opt_time(opt_time == 0) = NaN; %note that this removes sessions that were random but had no optimal window

% time spent in optimal window
figure
hold on
plotSpread(opt_time,'distributionColors',[0 0.6 1])
%plot(opt_time','Marker','o','LineWidth',0.5,'Color',[0 0.6 1]);
plot(nanmedian(opt_time,1),'Marker','o','LineWidth',1.5,'Color','b');
hold off
yline(20,':','LineWidth',2,'Color','k');
ylabel("Optimal window length (min)",'FontSize',12,'FontWeight','bold');
xlabel('Session','FontSize',12,'FontWeight','bold');
xlim([0 50]);
xticks([0:5:50]);
ylim([0 41]);
title('Optimal window duration','FontSize',14,'FontWeight','bold');

% time spent in optimal window in normalised training length
ci = 0.95;
alpha = 1 - ci;
T_multiplier = tinv((1-alpha/2), (numMice-1));
ci95 = T_multiplier*std(optnorm_time,0,1,'omitnan')/sqrt(numMice);

figure
hold on
plot(nanmedian(optnorm_time,1),'Marker','o','LineWidth',1.5,'Color','b');
errorbar(nanmedian(optnorm_time,1),ci95,'Marker','o','LineWidth',1.5,'Color','b');
hold off
yline(20,':','LineWidth',2,'Color','k');
ylabel("Optimal window length (min)",'FontSize',12,'FontWeight','bold');
xlabel('Training progression (%)','FontSize',12,'FontWeight','bold');
xlim([0 11]);
xticks([0:1:10]);
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
ylim([0 41]);
title('Optimal window duration','FontSize',14,'FontWeight','bold');

% all windows, one bin
% opt_time_all = reshape(opt_time,[],1);
% figure
% plotSpread(opt_time_all,'distributionColors',[0 0.6 1])
% violin(opt_time_all,'facecolor',[0 0.8 1],'edgecolor',[],'facealpha',0.3,'mc',[],'medc',[],'bw',0.5);
% errorbar(nanmedian(opt_time_all,1),std(opt_time_all,0,1,'omitnan'),'Marker','o','LineWidth',1.5,'Color','b');
% yline(20,':','LineWidth',2,'Color','k');
% ylabel("Optimal window length (min)",'FontSize',12,'FontWeight','bold');
% xlabel('Sessions (all)','FontSize',12,'FontWeight','bold');
% ylim([0 45]);
% xticks([]);
% title('Binned optimal window duration','FontSize',14,'FontWeight','bold');

%percentage correct within optimal window
% figure
% plot(opt_percor);
% errorbar(nanmean(opt_percor,1),std(opt_percor,0,1,'omitnan'),'Marker','o','LineWidth',1);
% ylabel("Percent correct (%)",'FontSize',12,'FontWeight','bold');
% xlabel('Session','FontSize',12,'FontWeight','bold');
% ylim([0 100]);
% title('Percent correct','FontSize',14,'FontWeight','bold');


%% Learned mice
learned_allmice = zeros(1,numMice);
num_days = learned_allmice;
learned_dur = cell(numMice,1);
minOptDur = 20; %the minimum duration of the optimal window for it to count as a good session
for m = 1:numMice
    goodSes = opt_time(m,:) > minOptDur;
    [labeledVector, numRegions] = bwlabel(goodSes);
    measurements = regionprops(labeledVector, goodSes, 'Area', 'PixelList');
    %test{m} = regionprops(labeledVector, goodSes, 'Area', 'PixelList');
    for r = 1:numRegions
        if measurements(r).Area >= 3
            num_days(m) = measurements(r).PixelList(1,1) + 2;
            learned_allmice(m) = 1;
            learned_dur{m} = opt_time(m,(num_days(m)+1):end); % OW durations after becoming expert
            break
        end
    end
end
learnedIndexes = find(learned_allmice);
learned = miceNames(learnedIndexes); % Mice that learned
learninglength = num_days(learnedIndexes); % Time taken to expert level



