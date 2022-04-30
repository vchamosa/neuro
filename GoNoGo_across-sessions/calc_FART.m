function [RT_FAs] = calc_FART(Data,metadata,varargin)
%% Calculates and plots false alarm reaction time across sessions of go/no-go behaviour
%
%  INPUTS
% Must include Data and metadata in order to compute reaction time. If no
% further specification is given it will do so to all mice provided by both
% inputs. If further specification is required, include name/value pairs of
% categories as found in metadata and animals.json. If filters are used in 
% combination the anaylsis will involve mice that satisfy both criteria.
%
% For example, calc_FART(Data,metadata,'responsible_person','Jeremiah', 
% 'strain','MCos') will plot the correct rejections of MCos mice trained by 
% Jeremiah.
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


% Compute false alarm reaction time per session
RT_FAs = NaN(numMice, 60);
for i = 1:numMice
    for day = 1 : numel(micedata{i}.Events)
        if ~isstring(micedata{i}.Events{day})
            if ~isempty(micedata{i}.Events{day})
                dayRT_FA = zeros(1, numel(micedata{i}.Events{day}.falseAlarms));
                for h = 1:numel(micedata{i}.Events{day}.falseAlarms)
                    % find the push that happens after the hth false alarm trial begins
                    test = micedata{i}.Events{day}.BL_on - micedata{i}.Events{day}.falseAlarms(h);
                    pushIndex = find(test>0,1);
                    if isempty(pushIndex) %|| test(pushIndex)>4000 % FALSE ALARMS NEEDS FIXING CAUSE THIS SHOULD NEVER HAPPEN
                        dayRT_FA(h) = NaN;
                    else
                        dayRT_FA(h) = test(pushIndex);
                    end
                end
                RT_FAs(i,day) = nanmedian(dayRT_FA);
            end
        end
    end
end


%plot
figure;
errorbar(nanmean(RT_FAs,1),std(RT_FAs,0,1,'omitnan'),'Marker','o','LineWidth',1);
yline(0,'--','LineWidth',1.5,'Color','k');
ylabel('Reaction time','FontSize',12,'FontWeight','bold');
xlabel('Session','FontSize',12,'FontWeight','bold');
title('Reaction time (false alarms)','FontSize',14,'FontWeight','bold');

% figure;
% plot(RT_FAs');


