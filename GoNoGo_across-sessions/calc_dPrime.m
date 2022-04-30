function [hitRate,falseAlarmRate,dPrime] = calc_dPrime(Data,metadata,varargin)
%% Calculates and plots d' across sessions go/no-go behaviour
%
%  INPUTS
% Must include Data and metadata in order to compute discrimination index.
% If no further specification is given it will do so to all mice provided 
% by both inputs. If further specification is required, include name/value 
% pairs of categories as found in metadata and animals.json. If filters are
% used in combination the anaylsis will involve mice that satisfy both
% criteria.
%
% For example, calc_dPrime(Data,metadata,'responsible_person','Henrietta', 
% 'strain','MCos') will plot the correct rejections of MCos mice trained by 
% Henrietta.
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


% Compute hit rate, false alarm rate and d' per session
hitRate = NaN(numMice, 60);
falseAlarmRate = NaN(numMice, 60);
hitRate_calc = NaN(numMice, 60);
falseAlarmRate_calc = NaN(numMice, 60);
for i = 1:numMice
    for day = 1 : numel(micedata{i}.Events)
        if ~isempty(micedata{i}.Events{day})
             if ~isstring(micedata{i}.Events{day}) 
                hitRate_calc(i,day) = (numel(micedata{i}.Events{day}.hits)+1)/((numel(micedata{i}.Events{day}.hits)+1)+(numel(micedata{i}.Events{day}.misses)+1)); % Constantinos method
                falseAlarmRate_calc(i,day) = (numel(micedata{i}.Events{day}.falseAlarms)+1)/((numel(micedata{i}.Events{day}.falseAlarms)+1)+(numel(micedata{i}.Events{day}.correctRejections)+1)); % Constantinos method
                 hitRate(i,day) = numel(micedata{i}.Events{day}.hits)/(numel(micedata{i}.Events{day}.hits)+numel(micedata{i}.Events{day}.misses)); % direct calculation
                 falseAlarmRate(i,day) = numel(micedata{i}.Events{day}.falseAlarms)/(numel(micedata{i}.Events{day}.falseAlarms)+numel(micedata{i}.Events{day}.correctRejections)); % direct calculation
%                 if numel(micedata{i}.Events{day}.hits) == 0 % Kuchibhotla method
%                     hitRate(i,day) = 1/(2 * (numel(micedata{i}.Events{day}.hits)+numel(micedata{i}.Events{day}.misses)));
%                 elseif numel(micedata{i}.Events{day}.misses) == 0
%                     hitRate(i,day) = 1 - (1/(2 * (numel(micedata{i}.Events{day}.hits)+numel(micedata{i}.Events{day}.misses))));
%                 else
%                     hitRate(i,day) = numel(micedata{i}.Events{day}.hits)/(numel(micedata{i}.Events{day}.hits)+numel(micedata{i}.Events{day}.misses));
%                 end
%                 if    numel(micedata{i}.Events{day}.correctRejections) == 0
%                     falseAlarmRate(i,day) = 1/(2 * (numel(micedata{i}.Events{day}.falseAlarms)+numel(micedata{i}.Events{day}.correctRejections)));
%                 elseif numel(micedata{i}.Events{day}.falseAlarms) == 0
%                     falseAlarmRate(i,day) = 1 - (1/(2 * (numel(micedata{i}.Events{day}.falseAlarms)+numel(micedata{i}.Events{day}.correctRejections))));
%                 else
%                     falseAlarmRate(i,day) = numel(micedata{i}.Events{day}.falseAlarms)/(numel(micedata{i}.Events{day}.falseAlarms)+numel(micedata{i}.Events{day}.correctRejections));
%                 end % end Kuchibhotla method
%                 hitRate(i,day) = numel(micedata{i}.Events{day}.hits)/(numel(micedata{i}.Events{day}.hits)+numel(micedata{i}.Events{day}.misses)); % direct calculation
%                 falseAlarmRate(i,day) = numel(micedata{i}.Events{day}.falseAlarms)/(numel(micedata{i}.Events{day}.falseAlarms)+numel(micedata{i}.Events{day}.correctRejections)); % direct calculation
             end
        end
    end
end
% hitRate(hitRate == 0) = 0.0001; % old method
% hitRate(hitRate == 1) = 0.999;
% falseAlarmRate(falseAlarmRate == 0) = 0.0001;
% falseAlarmRate(falseAlarmRate == 1) = 0.999; % end old method
 dPrime = norminv(hitRate_calc) - norminv(falseAlarmRate_calc);
% dPrime(dPrime == Inf) = NaN; % this ignores 1 and 0 values


%calculate 95% CI
ci = 0.95;
alpha = 1 - ci;
T_multiplier = tinv((1-alpha/2), (numMice-1));
ci95 = T_multiplier*std(dPrime,0,1,'omitnan')/sqrt(numMice);


%plot
figure;
hold on
plot(dPrime','Color',[0.85 0.85 0.85]);
%plot(nanmedian(dPrime),'Marker','o','LineWidth',1.5,'Color','k')
errorbar(nanmedian(dPrime,1),ci95,'Marker','o','LineWidth',1.5,'Color','k');
hold off
yline(0,':','LineWidth',2,'Color','k');
ylabel("d'",'FontSize',13,'FontWeight','bold');
xlabel('Session','FontSize',12,'FontWeight','bold');
ylim([-1.5 3]);
title("Discriminability index (d')",'FontSize',14,'FontWeight','bold');


