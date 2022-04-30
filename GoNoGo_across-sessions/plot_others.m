function plot_others(Data,metadata,numHits,numHits_norm,numCRs,numCRs_norm,hitRate,falseAlarmRate,percentCorrect,varargin)
%% Experimental, miscellaneous, lesser and WIP plots
%
%  INPUTS
% This function requires various inputs, please keep an eye out for where
% they're needed. It's recommended to comment out uninteresting plots
% and use this function as if it only offers one or two plots that are
% relevant to the moment. As said above, they are mostly experimental.
% Author: Victor Chamosa Pino
%


% Stop if not enough inputs supplied
if round(nargin/2) == nargin/2 || nargin == 1
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
     tMetadata = cell(1, numMice);
     for iMouse = 1 : numMice
         tMetadata{iMouse} = loadjson([miceNames{iMouse} '.json']);
     end
    
else
    % Quick adjustments if no further specification has been made
    numMice = numel(metadata);
    micedata = Data;
    tMetadata = cell(1, numMice);
    miceNames = cell(1, numMice);
    for iMouse = 1 : numMice
        miceNames{iMouse} = metadata{1,iMouse}.mouseID;
        tMetadata{iMouse} = loadjson([miceNames{iMouse} '.json']);
    end
end


%% Total rewards

numRew = numHits + numCRs;
numRew_norm = numHits_norm + numCRs_norm;

%calculate 95% CI
ci = 0.95;
alpha = 1 - ci;
T_multiplier = tinv((1-alpha/2), (numMice-1));
ci95 = T_multiplier*std(numRew,0,1,'omitnan')/sqrt(numMice);
ci95_hit = T_multiplier*std(numHits,0,1,'omitnan')/sqrt(numMice);
ci95_CR = T_multiplier*std(numCRs,0,1,'omitnan')/sqrt(numMice);
ci95_norm = T_multiplier*std(numRew_norm,0,1,'omitnan')/sqrt(numMice);

%plot
figure
hold on
plot(numRew','Color',[0.85 0.85 0.85]);
%plot(nanmedian(numRew),'Marker','o','LineWidth',1.5,'Color','k')
errorbar(nanmedian(numRew,1),ci95,'Marker','o','LineWidth',1.5,'Color','k');
hold off
ylabel('No. rewards','FontSize',12,'FontWeight','bold');
xlabel('Session','FontSize',12,'FontWeight','bold');
title('Total rewards','FontSize',14,'FontWeight','bold');

%separate hit & CR traces
figure
hold on
plot(numHits','Color',[1 0.6 0.6]);
plot(numCRs','Color',[0.6 0.6 1]);
errorbar(nanmedian(numHits,1),ci95_hit,'Marker','o','LineWidth',1.5,'Color','r');
errorbar(nanmedian(numCRs,1),ci95_CR,'Marker','o','LineWidth',1.5,'Color','b');
hold off
ylim([0 90]);
ylabel('No. rewards','FontSize',12,'FontWeight','bold');
xlabel('Session','FontSize',12,'FontWeight','bold');
title('Total rewards separated','FontSize',14,'FontWeight','bold');

% plot across normalised training time
figure
hold on
%plot(numRew_norm','Color',[0.85 0.85 0.85]);
errorbar(nanmedian(numRew_norm,1),ci95_norm,'Marker','o','LineWidth',1.5,'Color','k');
hold off
xlim([0 11]);
xticks([0:1:10]);
xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
ylabel('No. rewards','FontSize',12,'FontWeight','bold');
xlabel('Session','FontSize',12,'FontWeight','bold');
title('Total rewards','FontSize',14,'FontWeight','bold');


%% Hit percent

percentHits = NaN(numMice,60);
for i = 1:numMice
    for day = 1 : numel(micedata{i}.Events)
        if ~isnan(numHits(i,day))
            percentHits(i, day) = numHits(i,day)/(numHits(i,day)+numel(micedata{i}.Events{day}.misses));
        end
    end
end
percentHits(percentHits == 0) = 0.0001;
percentHits = percentHits*100;

%calculate 95% CI
ci95 = T_multiplier*std(percentHits,0,1,'omitnan')/sqrt(numMice);

%plot
figure;
hold on
plot(percentHits','Color',[0.85 0.85 0.85]);
%plot(nanmedian(percentCorrect),'Marker','o','LineWidth',1.5,'Color','k')
errorbar(nanmedian(percentHits,1),ci95,'Marker','o','LineWidth',1.5,'Color','k');
hold off
ylabel("Correct trials (%)",'FontSize',12,'FontWeight','bold');
xlabel('Session','FontSize',12,'FontWeight','bold');
ylim([0 100]);
title('Percent hits','FontSize',14,'FontWeight','bold');


%% CR percent

ratioCRs = NaN(numMice,60);
for i = 1:numMice
    for day = 1 : numel(micedata{i}.Events)
        if ~isnan(numCRs(i,day))
            ratioCRs(i,day) = numCRs(i,day)/(numCRs(i,day)+numel(micedata{i}.Events{day}.falseAlarms));
        end
    end
end
ratioCRs(ratioCRs == 0) = 0.1;
percentCRs = ratioCRs*100;

%calculate 95% CI
ci95 = T_multiplier*std(percentCRs,0,1,'omitnan')/sqrt(numMice);

%plot
figure;
hold on
plot(percentCRs','Color',[0.85 0.85 0.85]);
%plot(nanmedian(percentCorrect),'Marker','o','LineWidth',1.5,'Color','k')
errorbar(nanmedian(percentCRs,1),ci95,'Marker','o','LineWidth',1.5,'Color','k');
hold off
ylabel("Correct trials (%)",'FontSize',12,'FontWeight','bold');
xlabel('Session','FontSize',12,'FontWeight','bold');
ylim([0 100]);
title('Percent correct rejections','FontSize',14,'FontWeight','bold');


%% Performance segmented by tone length

tone4 = cell(numMice,1);
tone2 = cell(numMice,1);
tone250 = cell(numMice,1);
for i = 1:numMice
    for day = 1:numel(micedata{i}.Events)
        if day == 1
            tone4{i}.day(1) = day;
            tone4{i}.numHits(1) = numHits(i,day);
            tone4{i}.numCRs(1) = numCRs(i,day);
            tone4{i}.numRew(1) = numRew(i,day);
            tone4{i}.percentCorrect(1) = percentCorrect(i,day);
        elseif tMetadata{i}{day}.tone == 2 && tMetadata{i}{day-1}.tone == 4
            tone4{i}.day(2) = day-1;
            tone4{i}.numCRs(2) = numCRs(i,day-1);
            tone4{i}.numHits(2) = numHits(i,day-1);
            tone4{i}.numRew(2) = numRew(i,day-1);
            tone4{i}.percentCorrect(2) = percentCorrect(i,day-1);
            tone2{i}.day(1) = day;
            tone2{i}.numCRs(1) = numCRs(i,day);
            tone2{i}.numHits(1) = numHits(i,day);
            tone2{i}.numRew(1) = numRew(i,day);
            tone2{i}.percentCorrect(1) = percentCorrect(i,day);
        elseif (tMetadata{i}{day}.tone == 250 || tMetadata{i}{day}.tone == 240) && tMetadata{i}{day-1}.tone == 2
            tone2{i}.day(2) = day-1;
            tone2{i}.numCRs(2) = numCRs(i,day-1);
            tone2{i}.numHits(2) = numHits(i,day-1);
            tone2{i}.numRew(2) = numRew(i,day-1);
            tone2{i}.percentCorrect(2) = percentCorrect(i,day-1);
            tone250{i}.day(1) = day;
            tone250{i}.numCRs(1) = numCRs(i,day);
            tone250{i}.numHits(1) = numHits(i,day);
            tone250{i}.numRew(1) = numRew(i,day);
            tone250{i}.percentCorrect(1) = percentCorrect(i,day);
        elseif day == numel(micedata{i}.Events)
            tone250{i}.day(2) = day;
            tone250{i}.numCRs(2) = numCRs(i,day);
            tone250{i}.numHits(2) = numHits(i,day);
            tone250{i}.numRew(2) = numRew(i,day);
            tone250{i}.percentCorrect(2) = percentCorrect(i,day);
        end
    end
end

tones = NaN(numMice,8);
tone_day = NaN(numMice,2);
for i = 1:numMice
    tones(i,1) = tone4{i}.numRew(1);
    tones(i,2) = tone4{i}.numRew(2);
    tones(i,4) = tone2{i}.numRew(1);
    tones(i,5) = tone2{i}.numRew(2);
    tones(i,7) = tone250{i}.numRew(1);
    tones(i,8) = tone250{i}.numRew(2);
    tone_day(i,1) = tone4{i}.day(2);
    tone_day(i,2) = tone2{i}.day(2);
end

%calculate 95% CI
ci95 = T_multiplier*std(tones,0,1,'omitnan')/sqrt(numMice);

figure
hold on
for i = 1:numMice
    plot(tones(i,:),'Color',[0.85 0.85 0.85])
end
errorbar(nanmedian(tones),ci95,'LineWidth',2,'Color','k','Marker','o')
hold off
ylabel('No. rewards','FontSize',12,'FontWeight','bold')
xlabel('Session','FontSize',12,'FontWeight','bold')
xticks([1 2 4 5 7 8])
xticklabels({'First','Last','First','Last','First','Last'})
ylim([0 150])
title('Number of rewards across tone lengths','FontSize',14,'FontWeight','bold')

a = 2;

%% Internally-generated movements

% igm = nan(numMice,60);
% for i = 1:numMice
%     for day = 1:numel(micedata{i}.Events)
%         if isstring(micedata{i}.Events{day})
%             igm(i,day) = NaN;
%         elseif isempty(micedata{i}.Events{day})
%             igm(i,day) = NaN;
%         else
%             igm(i,day) = numel(micedata{i}.Events{day}.preCues);
%         end
%     end
% end
% 
% %calculate 95% CI
% ci95 = T_multiplier*std(igm,0,1,'omitnan')/sqrt(numMice);
% 
% %plot
% figure
% hold on
% plot(igm');
% errorbar(nanmedian(igm,1),ci95,'Marker','o','LineWidth',1.5,'Color','k');
% hold off
% ylim([0 200]);
% ylabel('Internally-generated movements','FontSize',12,'FontWeight','bold');
% xlabel('Session','FontSize',12,'FontWeight','bold');
% title('Internally-generated movements','FontSize',14,'FontWeight','bold');


%% Total pushes

% falseAlarms = nan(numMice,60);
% for i = 1:numMice
%     for day = 1:numel(micedata{i}.Events)
%         if isstring(micedata{i}.Events{day})
%             falseAlarms(i,day) = NaN;
%         elseif isempty(micedata{i}.Events{day})
%             falseAlarms(i,day) = NaN;
%         else
%             falseAlarms(i,day) = numel(micedata{i}.Events{day}.falseAlarms);
%         end
%     end
% end
% 
% numPush = numHits + igm + falseAlarms;
% 
% %calculate 95% CI
% ci95 = T_multiplier*std(numPush,0,1,'omitnan')/sqrt(numMice);
% 
% %plot
% figure;
% hold on
% plot(numPush');
% errorbar(nanmedian(numPush,1),ci95,'Marker','o','LineWidth',1.5,'Color','k');
% hold off
% ylim([0 250]);
% ylabel('No. pushes','FontSize',12,'FontWeight','bold');
% xlabel('Session','FontSize',12,'FontWeight','bold');
% title('Total pushes','FontSize',14,'FontWeight','bold');


%% Weight

% weights = [17.76 18.51 17.41 17.65 24.4 24.1 22.8 22.1 22.1 24.35 22.8 21.76 22.47] % Initial weights for 13 of Victor's mice. Clunky, I know.
% weights = weights'
% 
% Weight = NaN(numMice,60);
% for i = 1:numMice
%     for day = 1 : numel(tMetadata{i})
%         if ~isstring(tMetadata{i}{day}.weight)
%             if ~isempty(tMetadata{i}{day}.weight)
%                 Weight(i, day) = tMetadata{i}{day}.weight;
%             end
%         end
%     end
% end
% perWeight = 100*(Weight./weights)
% 
% %calculate 95% CI
% ci95 = T_multiplier*std(perWeight,0,1,'omitnan')/sqrt(numMice);
% 
% %plot
% figure;
% hold on
% plot(perWeight');
% %plot(nanmedian(perWeight),'Marker','o','LineWidth',1.5,'Color','k')
% errorbar(nanmedian(perWeight,1),ci95,'Marker','o','LineWidth',1.5,'Color','k');
% hold off
% yline(85,':','LineWidth',2,'Color','r');
% ylabel("Bodyweight (%)",'FontSize',12,'FontWeight','bold');
% xlabel('Session','FontSize',12,'FontWeight','bold');
% ylim([70 135]);
% title('Weight progression','FontSize',14,'FontWeight','bold');


