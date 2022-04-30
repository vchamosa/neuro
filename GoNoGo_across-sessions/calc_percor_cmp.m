function [percor_cmp,p_percor] = calc_percor_cmp(Data,metadata,hitRate_cmp,varargin);
%% Compares percentage of correct trials across sessions of go/no-go behaviour
%
%  INPUTS
% Must include Data, metadata and hitRate_cmp in order to compute 
% percentage correct trials. It also requires further specification,
% but can only take in a single name/value pair as per animals.json. 
% However, one can input as many values as desired for comparison, within 
% curly brackets. 
% Importantly, the inputted values must be in the same order as in
% calc_dPrime_cmp, so that the hitRate values are in the same order too.
% This script currently does not account for the value itself, but its
% position in the list you provide.
%
% For example, calc_percor_cmp(Data,metadata,hitRate_cmp,'strain',{'MCos',
% 'VGAT','Thy1'})will plot a single figure containing the percentage of 
% correct trials across sessions of each strain in the list.
% Author: Victor Chamosa Pino
%


% Stop if not enough inputs supplied
if round(nargin/2) == nargin/2 || nargin == 1
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
        hitRate = hitRate_cmp{cmp};
        percentCorrect = calc_percor(Data,metadata,hitRate,inputs{1,1},inputs{cmp+1,1});
        percor_cmp{cmp,1} = percentCorrect;
    end
    
    %plot figure
    set(0,'DefaultFigureVisible','on')
    Colour = distinguishable_colors(cmp);
    figure;
    hold on
    for plt = 1:length(percor_cmp(:,1))
        sem = std(percor_cmp{plt,1},0,1,'omitnan')/sqrt(length(percor_cmp{plt,1})); %%%
        errorbar(nanmedian(percor_cmp{plt,1},1),sem,'Marker','o','LineWidth',1,'Color',Colour(plt,:)); %%%
    end
    ylabel("Correct trials (%)",'FontSize',12,'FontWeight','bold');
    xlabel('Session','FontSize',12,'FontWeight','bold');
    ylim([0 100]);
    title('Percent correct','FontSize',14,'FontWeight','bold');
    hold off
    legend(input{2,1});
    
    %stats
    for stat = 1:length(percor_cmp(:,1))
        percor_cmp_st(:,stat) = nanmean(percor_cmp{stat}(:,:));
    end
    p_percor = friedman(percor_cmp_st,[],'off');
end


    
