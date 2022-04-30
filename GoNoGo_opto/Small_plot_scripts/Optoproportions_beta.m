%% Median proportions of success in control vs laser trials
% Compendium of boxplots of the data calculated
% by the GoNoGo_SingleSession_opto.m script
% Author: Victor Chamosa Pino
%


% optCR_proportions = [0.7895,0.8158; 1,1; 1,0.8; 0.8269,0.8235];
% optCR_meanprop = round(mean(optCR_proportions),3);
% optFA_meanprop = 1 - optCR_meanprop;
% optFA_meanprop = round(optFA_meanprop,3);
% optNOGO_meanprop = vertcat(optCR_meanprop,optFA_meanprop);
% trial_categories = categorical({'Control','Laser'});
% 
% figure
% hold on
% bar(trial_categories,optNOGO_meanprop','stacked');
% for plt = 1:length(optCR_proportions)
%     plot(optCR_proportions(plt,:),'Marker','.','Color','k')
% end
% hold off


% No-go data for C57 483, 485, 744, 745, 748, 749
%optCR_proportions = [1,0.56; 0.9754,0.7647; 0.875,0.4286; 0.76,0.82; 0.9305,0.9677; 0.7391,0.6392; 0.8,0.6154; 0.75,0.28];

% No-go data for the first somatic sessions of C57 483, 485, 748, 749, 1024
%optCR_proportions = [1,0.56; 0.9754,0.7647; 0.875,0.4286; 0.9305,0.9677; 0.8,0.6154; 0.9297,0.7974];
%optFA_proportions = 1 - optCR_proportions;
%optCR_hits = [0.888; 0.917; 1; 1; 1];

% No-go data for the summed first somatic sessions of C57 483, 485, 748, 749, 259, 1024
%optCR_proportions = [1,0.5833; 0.9444,0.7619; 0.931,0.9583; 0.8333,0.6; 0.9302,0.7647; 0.7857,0.6667];
%optFA_proportions = 1 - optCR_proportions;

% No-go data for the median first somatic sessions of C57 483, 485, 748 (1,1;), 749, 259, 1024 
optCR_proportions = [1,0.5833; 1,0.7778; 0.875,0.5833; 1,0.7083; 0.7708,0.625];
optFA_proportions = 1 - optCR_proportions;

% Go data for the summed somatic sessions of C57 745, 747, 748, 494, 495, 1025
optHit_proportions = [0.9444,0.9318; 0.4026,0.5606; 0.871,0.8302; 0.9744,1; 1,1; 0.8776,0.9111];

% Go data for the median somatic sessions of C57 745, 747, 748, 494, 495, 1025
%optHit_proportions = [1,1; 1,1; 1,1; 1,1; 1,1; 0.9375,1];

% Early push only data
%earlyHit_proportions = [32,33; 31,31; 36,18; 27,22; 0,0; 2,3]; % Raw number
% Removed datapoints: 0.2,0.9167; 0.2053,0.9394 (for reasons of abnormal shutter)
earlyHit_proportions = [0.3303,0.1698; 0.2596,0.2056; 0,0; 0.0172,0.0254]; % Proportions

% optCR_meanprop = round(mean(optCR_proportions),3);
% optFA_meanprop = 1 - optCR_meanprop;
% optFA_meanprop = round(optFA_meanprop,3);
% optNOGO_meanprop = vertcat(optCR_meanprop,optFA_meanprop);
trial_categories = categorical({'Control','Laser'});

% Go/no-go boxplot proportions
figure
hold on
for plt = 1:length(optCR_proportions)
    plot(optCR_proportions(plt,:),'Marker','.','MarkerSize',9,'Color',[0.75 0.75 0.75])
end
boxplot(optCR_proportions,trial_categories);
hold off
ylim([0 1]);
title('No-go','FontSize',14,'FontWeight','bold');

p = ranksum(optCR_proportions(:,1),optCR_proportions(:,2))


% Go boxplot proportions
figure
hold on
for plt = 1:length(optHit_proportions)
    plot(optHit_proportions(plt,:),'Marker','.','MarkerSize',9,'Color',[0.75 0.75 0.75])
end
boxplot(optHit_proportions,trial_categories);
hold off
ylim([0 1]);
title('Go','FontSize',14,'FontWeight','bold');


% Go-only boxplot proportions
figure
hold on
for plt = 1:length(earlyHit_proportions)
    plot(earlyHit_proportions(plt,:),'Marker','.','MarkerSize',9,'Color',[0.75 0.75 0.75])
end
boxplot(earlyHit_proportions,trial_categories);
hold off
ylim([0 1]);
title('Early pushes','FontSize',14,'FontWeight','bold');


% ROC space proportions
% cutoff = 1.5;
% ROC_proportions = figure('Renderer', 'painters'); % comment out the parentheses if figure generation is problematic
% axis([0, 1, 0, 1]);
% patch([0 0.5 0], [0 0.5 0.5], 'blue', 'EdgeAlpha',0, 'FaceAlpha', 0.2);
% patch([0 0.5 0.5 0], [0.5 0.5 1 1], 'green', 'EdgeAlpha',0, 'FaceAlpha', 0.2);
% patch([0.5 1 0.5], [0.5 1 1], 'red', 'EdgeAlpha',0, 'FaceAlpha', 0.2);
% patch([0 1 1], [0 0 1], 'black', 'EdgeAlpha',0, 'FaceAlpha', 0.2);
% hold on
% X = [0.01:0.01:0.99];
% Y = [0.01:0.01:0.99];
% [X,Y] = meshgrid(X,Y);
% dPrimeContour = norminv(Y)-norminv(X); % Y = Hit Rate, X = False Alarm Rate
% V1 = [cutoff,cutoff]; % set the desired cutoff value here
% contour(X,Y,dPrimeContour,V1,'LineColor','k','LineStyle','--','LineWidth',1);
% plotsize = 50;
% scatter(optFA_proportions(:,1),optCR_hits,plotsize,[0.8500 0.3250 0.0980],'filled','MarkerEdgeColor','black');
% scatter(optFA_proportions(:,2),optCR_hits,plotsize,[0 0.4470 0.7410],'filled','MarkerEdgeColor','black');
% hold off
% xlabel('False Alarm rate');
% ylabel('Hit Rate');
% title('Control/laser optimal rates comparison','FontSize',14,'FontWeight','bold');
% 



