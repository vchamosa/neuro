%% Reaction times of various control and laser trials
%
% Compendium of reaction times in push trials, calculated
% by the GoNoGo_SingleSession_opto.m script
% Author: Victor Chamosa Pino
%

clear all

seshes = load('reaction_times.mat'); % Saved to Ardbeg
sessions = struct2cell(seshes);
% optRT_FA  optRT_FA_laser  optRT_hit


FA_control = [];
FA_laser = [];
hit_control = [];
for rt = 1:length(sessions)
    FA_control = [FA_control sessions{rt}{1,1}];
    FA_laser = [FA_laser sessions{rt}{1,2}];
    hit_control = [hit_control sessions{rt}{1,3}];
end


figure
hold on
histogram(FA_control,'BinWidth',0.05,'FaceColor',[0 0.68 0],'Normalization','probability')
histogram(FA_laser,'BinWidth',0.05,'FaceColor',[0.8500 0.3250 0.0980],'Normalization','probability')
histogram(hit_control,'BinWidth',0.05,'FaceColor',[0 0.3 1],'Normalization','probability')
title('Reaction times','FontSize',14,'FontWeight','bold');
hold off
ylim([0 1]);
xlim([0 2]);
xlabel('Time (seconds)','FontSize',13,'FontWeight','bold');



