function [optRew_cmp,optTime_cmp,optPercor_cmp,learned_cmp,learninglength_cmp,p_optTime] = calc_optwin_cmp(Data,metadata,varargin)
%% Compares optimal window parameters across sessions of go/no-go behaviour
%
%  INPUTS
% Must include Data and metadata in order to compute optimal window stats.
% It also requires further specification, but can only take in a single
% name/value pair as per animals.json. However, one can input as many
% values as desired for comparison, within curly brackets.
%
% For example, calc_optwin_cmp(Data,metadata,'strain',{'MCos','VGAT','Thy1'})
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
        [opt_rew,opt_time,opt_percor,learned,learninglength] = calc_optwin(Data,metadata,inputs{1,1},inputs{cmp+1,1});
        optRew_cmp{cmp,1} = opt_rew;
        optTime_cmp{cmp,1} = opt_time;
        optPercor_cmp{cmp,1} = opt_percor;
        learned_cmp{cmp,1} = learned;
        learninglength_cmp{cmp,1} = learninglength;
    end

    
    %plots
    set(0,'DefaultFigureVisible','on')
    Colour = distinguishable_colors(cmp);

    %time spent in optimal windows
    figure;
    hold on
    for idx = 1:length(optTime_cmp(:,1))
        plotSpread(optTime_cmp{idx,1},'distributionColors',Colour(idx,:))
        plot(nanmedian(optTime_cmp{idx,1},1),'Marker','o','LineWidth',1.5,'Color',Colour(idx,:)/2);
    end
    hold off
    yline(20,':','LineWidth',2,'Color','k');
    ylabel("Optimal window length (min)",'FontSize',12,'FontWeight','bold');
    xlabel('Session','FontSize',12,'FontWeight','bold');
    xticks([0:10:60]);
    xlim([0 50]);
    ylim([0 41]);
    title('Optimal window duration','FontSize',14,'FontWeight','bold');
    p = findobj(gca,'Type','line','Marker','o');
    legend([p],flip(input{2,1}),'Location','southeast');  
    
    
    %all windows, one bin
%     figure;
%     hold on 
%     for idx = 1:length(optTime_cmp(:,1))
%         optTime_cmp_all = reshape(optTime_cmp{idx,1},[],1);
%         plotSpread(optTime_cmp_all,'distributionColors',Colour(idx,:))
%         violin(optTime_cmp_all,'facecolor',Colour(idx,:),'edgecolor',Colour(idx,:)/2,'facealpha',0.2,'mc',[],'medc',[],'bw',1);
%         errorbar(nanmedian(optTime_cmp_all,1),std(optTime_cmp_all,0,1,'omitnan'),'Marker','o','LineWidth',1.5,'Color',Colour(idx,:)/2);
%     end
%     hold off
%     yline(20,':','LineWidth',2,'Color','k');
%     ylabel("Optimal window length (min)",'FontSize',12,'FontWeight','bold');
%     xlabel('Sessions (all)','FontSize',12,'FontWeight','bold');
%     ylim([0 45]);
%     xticks([]);
%     title('Binned optimal window duration','FontSize',14,'FontWeight','bold');
%     p = findobj(gca,'Type','line');
%     legend([p],flip(input{2,1}),'Location','southeast');
    
    
    %percentage correct within optimal window
    figure;
    hold on
    for idx = 1:length(optPercor_cmp(:,1))
    errorbar(nanmean(optPercor_cmp{idx,1},1),std(optPercor_cmp{idx,1},0,1,'omitnan'),'Marker','o','LineWidth',1,'Color',Colour(idx,:));
    end
    ylabel("Percent correct (%)",'FontSize',12,'FontWeight','bold');
    xlabel('Session','FontSize',12,'FontWeight','bold');
    ylim([0 100]);
    title('Percent correct','FontSize',14,'FontWeight','bold');
    hold off
    legend(input{2,1});
    
    
    %stats
    for stat = 1:length(optTime_cmp(:,1))
        optTime_cmp_all_st(:,stat) = nanmean(optTime_cmp{stat}(:,:));
    end
        p_optTime = friedman(optTime_cmp_all_st,2);
    
end


    
