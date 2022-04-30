function [numCRs,numCRs_norm] = calc_numCRs(Data,metadata,varargin)
%% Calculates and plots number of correct rejections across sessions of go/no-go behaviour
%
%  INPUTS
% Must include Data and metadata in order to compute number of correct 
% rejections. If no further specification is given it will do so to all 
% mice provided by both inputs. If further specification is required, 
% include name/value pairs of categories as found in metadata and 
% animals.json. If filters are used in combination the anaylsis will 
% involve mice that satisfy both criteria.
%
% For example, calc_numCRs(Data,metadata,'responsible_person','Woodrow', 
% 'strain','MCos') will plot the correct rejections of MCos mice trained by 
% Woodrow.
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


% Compute number of correct rejections per session
numCRs = NaN(numMice,60);
for i = 1:numMice
    for day = 1 : numel(micedata{i}.Events)
        if ~isstring(micedata{i}.Events{day})
            if ~isempty(micedata{i}.Events{day})
                numCRs(i,day) = numel(micedata{i}.Events{day}.correctRejections);
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
        numCRs_norm(i,sn) = numCRs(i,session_unit);
    end
end


%calculate 95% CI
% ci = 0.95 ;
% alpha = 1 - ci;
% T_multiplier = tinv((1-alpha/2), (numMice-1));
% ci95 = T_multiplier*std(numCRs,0,1,'omitnan')/sqrt(numMice);
% 
% %plot
% figure
% hold on
% plot(numCRs','Color',[0.85 0.85 0.85]);
% %plot(nanmedian(numCRs),'Marker','o','LineWidth',1.5,'Color','k')
% errorbar(nanmedian(numCRs,1),ci95,'Marker','o','LineWidth',1.5,'Color','k');
% hold off
% ylabel('No. rewards','FontSize',12,'FontWeight','bold');
% xlabel('Session','FontSize',12,'FontWeight','bold');
% ylim([0 130]);
% title('Correct rejections','FontSize',14,'FontWeight','bold');
% 
% 
% % plot CRs across normalised training time
% ci = 0.95;
% alpha = 1 - ci;
% T_multiplier = tinv((1-alpha/2), (numMice-1));
% ci95 = T_multiplier*std(numCRs_norm,0,1,'omitnan')/sqrt(numMice);
% 
% figure
% hold on
% plot(numCRs_norm','Color',[0.85 0.85 0.85]);
% errorbar(nanmedian(numCRs_norm,1),ci95,'Marker','o','LineWidth',1.5,'Color','k');
% hold off
% ylabel('No. rewards','FontSize',12,'FontWeight','bold');
% xlabel('Session','FontSize',12,'FontWeight','bold');
% ylim([0 130]);
% xlim([0 11]);
% xticks([0:1:10]);
% xticklabels({'0','10','20','30','40','50','60','70','80','90','100'})
%title('Correct rejections','FontSize',14,'FontWeight','bold');



