function [numHits_cmp] = calc_numHits_cmp(Data,metadata,varargin)
%% Compares hits across sessions of go/no-go behaviour
%
%  INPUTS
% Must include Data and metadata in order to compute number of hits. It
% also requires further specification, but can only take in a single
% name/value pair as per animals.json. However, one can input as many
% values as desired for comparison, within curly brackets.
%
% For example, calc_numHits_cmp(Data,metadata,'strain',{'MCos','VGAT','Thy1'})
% will plot a single figure containing the hits across sessions of each
% strain in the list.
% Author: Victor Chamosa Pino
%


% Stop if not enough inputs supplied
if round(nargin/2) ~= nargin/2
    error('Please ensure you are providing sufficient input');
    return
end


% Parse optional inputs
if nargin > 3
    input = reshape(varargin,2,[]);
    for var = 1:length(input(1,:))
        inputs = cell(numel(input{2,var})+1,length(input(1,:)));
    end
    for var = 1:length(input(1,:))
        if ~iscell(input{2,var})
            inputs(:,var) = input(2,var)
        else
            for val = 1:numel(input{2,var})
                inputs{val+1,var} = input{2,var}(val);
            end
        end
    end
    %call function
    inputs(1,:) = input(1,:);
    set(0,'DefaultFigureVisible','off')
    for cmp = 1:(length(inputs(:,1))-1)
        numHits = calc_numHits(Data,metadata,inputs{1,1},inputs{cmp+1,1});
        numHits_cmp{cmp,1} = numHits;
    end
    %plot figure
    set(0,'DefaultFigureVisible','on')
    figure;
    hold on
    for plot = 1:length(numHits_cmp(:,1))
        errorbar(nanmean(numHits_cmp{plot,1},1),std(numHits_cmp{plot,1},0,1,'omitnan'),'Marker','o','LineWidth',1,'Color',[rand,rand,rand]);
    end
    ylabel('Hits','FontSize',12,'FontWeight','bold');
    xlabel('Session','FontSize',12,'FontWeight','bold');
    title('Number of hits','FontSize',14,'FontWeight','bold');
    hold off
    legend(input{2,1});
end


    
