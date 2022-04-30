function [RT_hits] = calc_hitRT(Data,metadata,varargin)
%% Calculates and plots hit reaction time across sessions of go/no-go behaviour
%
%  INPUTS
% Must include Data and metadata in order to compute reaction time. If no
% further specification is given it will do so to all mice provided by both
% inputs. If further specification is required, include name/value pairs of
% categories as found in metadata and animals.json. If filters are used in 
% combination the anaylsis will involve mice that satisfy both criteria.
%
% For example, calc_hitRT(Data,metadata,'responsible_person','Dorothy', 
% 'strain','MCos') will plot the correct rejections of MCos mice trained by 
% Dorothy.
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
else
    % Quick adjustments if no further specification has been made
    numMice = numel(metadata);
    micedata = Data;
end


% Compute hit reaction time per session
RT_hits = NaN(numMice,60);
RT_hits2 = NaN(numMice,60);
for i = 1:numMice
    for day = 1 : numel(micedata{i}.Events)
        if ~isstring(micedata{i}.Events{day})
            if ~isempty(micedata{i}.Events{day})
                %dayRT_hit = zeros(1, numel(micedata{i}.Events{day}.hits));
                dayRT_hits = zeros(1, numel(micedata{i}.Events{day}.hits));
                samplerate = round(micedata{i}.Events{day}.sampleRate);
                samplerate_test(i,day) = samplerate;
                for h = 1:numel(micedata{i}.Events{day}.hits)
                    % find the push that happens after the hth hit trial begins
                    test = micedata{i}.Events{day}.BL_on - micedata{i}.Events{day}.hits(h);
                    pushIndex = find(test>0,1);
                    if isempty(pushIndex)
                        %dayRT_hit(h) = NaN; %fix thiiiisss
                        dayRT_hits(h) = NaN;
                    else
                        %dayRT_hit(h) = test(pushIndex);
                        dayRT_hits(h) = test(pushIndex)/samplerate;
                    end
                end
                %RT_hits(i,day) = nanmedian(dayRT_hit);
                RT_hits2(i,day) = nanmedian(dayRT_hits);
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
        RT_hits_norm(i,sn) = RT_hits2(i,session_unit);
    end
end

 
% %plot
% figure
% hold on
% plot(RT_hits2');
% plot(nanmedian(RT_hits2),'Marker','o','LineWidth',1.5,'Color','k')
% hold off
% %errorbar(nanmean(RT_hits2,1),std(RT_hits,0,1,'omitnan'),'Marker','o','LineWidth',1);
% %yline(0,'--','LineWidth',1.5,'Color','k');
% ylabel('Reaction time','FontSize',12,'FontWeight','bold');
% xlabel('Session','FontSize',12,'FontWeight','bold');
% title('Reaction time (hits)','FontSize',14,'FontWeight','bold');


%calculate 95% CI
ci = 0.95;
alpha = 1 - ci;
T_multiplier = tinv((1-alpha/2), (numMice-1));
ci95 = T_multiplier*std(RT_hits_norm,0,1,'omitnan')/sqrt(numMice);

%plot normalised to training length
figure
hold on
plot(RT_hits_norm');
plot(nanmedian(RT_hits_norm),'Marker','o','LineWidth',1.5,'Color','k')
errorbar(nanmedian(RT_hits_norm,1),ci95,'Marker','o','LineWidth',1.5,'Color','k');
hold off
%errorbar(nanmean(RT_hits,1),std(RT_hits,0,1,'omitnan'),'Marker','o','LineWidth',1);
%yline(0,'--','LineWidth',1.5,'Color','k');
ylim([0 1.5])
ylabel('Reaction time (s)','FontSize',12,'FontWeight','bold');
xlabel('Session','FontSize',12,'FontWeight','bold');
title('Reaction time (hits)','FontSize',14,'FontWeight','bold');



