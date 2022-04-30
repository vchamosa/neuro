function [dPrime_cmp,hitRate_cmp,falseAlarmRate_cmp,p_dPrime] = calc_dPrime_cmp(Data,metadata,varargin)
%% Compares discrimination index values across sessions of go/no-go behaviour
%
%  INPUTS
% Must include Data and metadata in order to compute discrimination index.
% It also requires further specification, but can only take in a single
% name/value pair as per animals.json. However, one can input as many
% values as desired for comparison, within curly brackets.
%
% For example, calc_dPrime_cmp(Data,metadata,'strain',{'MCos','VGAT','Thy1'})
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
        [hitRate,falseAlarmRate,dPrime] = calc_dPrime(Data,metadata,inputs{1,1},inputs{cmp+1,1});
        dPrime_cmp{cmp,1} = dPrime;
        hitRate_cmp{cmp,1} = hitRate;
        falseAlarmRate_cmp{cmp,1} = falseAlarmRate;
    end
    
    %plot figure
    set(0,'DefaultFigureVisible','on')
    Colour = distinguishable_colors(cmp);
    figure;
    hold on
    for plt = 1:length(dPrime_cmp(:,1))
        sem = std(dPrime_cmp{plt,1},0,1,'omitnan')/sqrt(length(dPrime_cmp{plt,1})); %%%
        errorbar(nanmedian(dPrime_cmp{plt,1},1),sem,'Marker','o','LineWidth',1,'Color',Colour(plt,:));
    end
    yline(0,':','LineWidth',0.75,'Color','k');
    ylabel("d'",'FontSize',13,'FontWeight','bold');
    xlabel('Session','FontSize',12,'FontWeight','bold');
    ylim([-1.5 2.5]);
    title("Discriminability index (d')",'FontSize',14,'FontWeight','bold');
    hold off
    legend(input{2,1});
    
    %stats
    for stat = 1:length(dPrime_cmp(:,1))
        dPrime_cmp_st(:,stat) = nanmean(dPrime_cmp{stat}(:,:));
        for t = 1: length(dPrime_cmp_st(:,stat))
        if isnan(dPrime_cmp_st(t,stat))
            dPrime_cmp_st(t,stat) = 0;
        end
        end
    end
    p_dPrime = friedman(dPrime_cmp_st,2,'off');
    
end
    
