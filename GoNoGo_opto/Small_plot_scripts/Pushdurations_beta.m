%% Push durations of various control and laser trials
%
% Compendium of push durations in push trials, calculated
% by the GoNoGo_SingleSession_opto.m script
% Author: Victor Chamosa Pino
%

clear all

seshes = load('push_durations.mat'); % Saved to Ardbeg
sessions = struct2cell(seshes);
% opt_pushdur_FA opt_pushdur_FA_laser opt_pushdur_hit

FA_control = [];
FA_laser = [];
hit_control = [];
for push = 1:length(sessions)
    FA_control = [FA_control sessions{push}{1,1}];
    FA_laser = [FA_laser sessions{push}{1,2}];
    hit_control = [hit_control sessions{push}{1,3}];
end


figure
hold on
histogram(FA_control,'BinWidth',0.1,'FaceColor',[0 0.68 0],'Normalization','probability')
histogram(FA_laser,'BinWidth',0.1,'FaceColor',[0.8500 0.3250 0.0980],'Normalization','probability')
histogram(hit_control,'BinWidth',0.1,'FaceColor',[0 0.3 1],'Normalization','probability')
title('Push duration','FontSize',14,'FontWeight','bold');
hold off
ylim([0 1]);
xlim([0 4]);
xlabel('Time (seconds)','FontSize',13,'FontWeight','bold');



