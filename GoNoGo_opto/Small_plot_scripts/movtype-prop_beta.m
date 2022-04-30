%% To plot the proportions of movement types across opto trials
%
% Very temporary script, work in progress!!!! Please check for warning
% comments below!!!!!!!!!
%
% This plots proportions of qualitative forepaw movement types in direct
% response to various cues, as established by observing the videos of
% go/no-go opto sessions.
% The data can be found in a more qualitatively detailed/raw state in Ardbeg and at:
% https://docs.google.com/spreadsheets/d/1yQXitGKtN0I9zGdh9hVuiVvINgywoSeq0vG0hmD1WKs/edit?usp=sharing
% Author: Victor Chamosa Pino
%

trialTypes{1} = [7,5,5,5,2,5,1,5,5,4,5,1,2,5,1,5,5,5,1,5,5,1,4,1,5,5,5,1,7,7,5,7,5,5,5,5,5,7,5,5,5,5,1,1,5,1,1,5,2,5,5,7,5,2,5,5,1,1,5,4,5,5,5,1,5,5,1,5,2,4,4,5,5,2,5,2,1];
trialTypes{2} = [1,7,1,2,1,5,1,5,5,1,7,5,1,1,2,5,5,5,4,7,5,5,2,1,2,2,1,1,5,5,2,2,5,5,2,2,1,1,1,5,5,5,1,1,5,5,4,1,5,5,4,2,5,5,4,2,2,1,2,2,5,5,5,5,3,1,5,5,2,5,5,4,5,5,5,2];
trialTypes{3} = [5,5,5,4,5,1,5,2,5,5,7,5,2,5,1,1,1,1,1,1,5,5,1,1,1,1,2,2,1,5,5,2,1,5,4,2,1,5,1,2,1,5,5,3,7,5,4,2,5,1,5,1,5,5,7,5,5,5,1,5,4,1,5,5,5,5,5,1,5,1,5,2,1,5,2,2,5,5,5,2,1,5,5,2,1,5,1,5,5,1,5,7,1,1,4,5,5,5,4,5,5,1,2,1,4,5,5,5,5,1,5,2,5,5,5,5,5,2,2,1,1,3,5,1,5,1,3,5,5,5,2,2,5,4,2,5,5,2,2,2,1,2,2,5,5,5,5,1,2,5,1,5,1,2,5];
trialTypes{4} = [5,5,1,1,5,2,5,5,3,5,1,5,5,5,5,4,5,5,1,5,3,5,1,5,5,3,5,1,2,5,1,1,2,2,5,5,4,7,1,7,1,5,1,5,5];

moveType{1} = [13,9,9,9,13,9,12,9,9,9,9,11,12,9,13,9,9,9,11,9,9,13,10,11,10,9,9,11,13,13,9,13,9,9,9,9,9,12,9,9,9,9,11,13,9,12,11,9,12,9,9,12,9,11,9,9,11,11,9,9,9,9,9,11,9,9,13,9,13,10,10,9,9,13,9,11,12];
moveType{2} = [12,12,9,11,11,9,11,9,9,12,13,9,11,11,11,9,9,9,10,13,9,9,13,11,11,11,11,11,9,9,13,12,9,9,11,11,11,11,12,9,9,10,11,12,10,9,10,11,10,9,10,11,10,9,10,13,11,11,11,13,10,9,9,10,10,11,9,9,11,9,9,9,9,9,9,11];
moveType{3} = [9,9,9,9,9,11,9,13,9,9,11,9,12,9,13,11,13,12,11,13,9,9,12,11,12,12,12,12,11,9,9,11,12,9,10,11,13,9,12,11,11,9,9,9,12,9,9,11,9,11,10,13,10,10,13,9,10,9,11,9,9,11,9,9,9,9,9,11,9,13,9,12,11,9,11,11,9,9,9,12,11,9,9,11,11,10,13,9,9,13,10,11,11,12,9,9,9,9,9,9,9,13,12,11,12,9,9,9,9,13,9,13,9,9,10,9,9,12,13,12,11,10,9,12,9,13,13,9,9,9,13,11,9,9,13,9,9,13,13,13,11,12,12,9,9,9,9,13,13,9,13,9,13,13,9];
moveType{4} = [9,9,13,13,10,11,9,10,9,10,11,9,10,9,10,9,9,10,11,9,9,10,12,9,9,10,10,13,11,9,11,11,13,11,10,9,10,11,11,11,11,10,11,10,9];


ctrl_CR_still = 0;
ctrl_CR_length = 0;
ctrl_CR_repos = 0;
laser_CR_still = 0;
laser_CR_length = 0;
laser_CR_repos = 0;

ctrl_hit_repopush = 0;
ctrl_hit_simpush = 0;
laser_hit_repopush = 0;
laser_hit_simpush = 0;

ctrl_FA_repopush = 0;
ctrl_FA_simpush = 0;
laser_FA_repopush = 0;
laser_FA_simpush = 0;

ctrl_miss_still = 0;
ctrl_miss_length = 0;
ctrl_miss_repos = 0;
laser_miss_still = 0;
laser_miss_length = 0;
laser_miss_repos = 0;

for ses = 1:4
    ctrlCR_still(ses) = 0;
    for i = 1:length(trialTypes{ses})
        if trialTypes{ses}(i) == 1 && moveType{ses}(i) == 11 % Control CR & still
            ctrl_CR_still = ctrl_CR_still + 1;
            ctrlCR_still(ses) = ctrlCR_still(ses) + 1;
            
        elseif trialTypes{ses}(i) == 1 && moveType{ses}(i) == 12 % Control CR & lengthwise
            ctrl_CR_length = ctrl_CR_length + 1;
            
        elseif trialTypes{ses}(i) == 1 && moveType{ses}(i) == 13 % Control CR & reposition
            ctrl_CR_repos = ctrl_CR_repos + 1;
            
        elseif trialTypes{ses}(i) == 2 && moveType{ses}(i) == 11 % Laser CR & still
            laser_CR_still = ctrl_CR_still + 1;
            
        elseif trialTypes{ses}(i) == 2 && moveType{ses}(i) == 12 % Laser CR & lengthwise
            laser_CR_length = ctrl_CR_length + 1;
            
        elseif trialTypes{ses}(i) == 2 && moveType{ses}(i) == 13 % Laser CR & reposition
            laser_CR_repos = ctrl_CR_repos + 1;
            
        elseif trialTypes{ses}(i) == 3 && moveType{ses}(i) == 9 % Control FA & reposition + push
            ctrl_FA_repopush = ctrl_FA_repopush + 1;
            
        elseif trialTypes{ses}(i) == 3 && moveType{ses}(i) == 10 % Control FA & simple push
            ctrl_FA_simpush = ctrl_FA_simpush + 1;
            
        elseif trialTypes{ses}(i) == 4 && moveType{ses}(i) == 9 % Laser FA & reposition + push
            laser_FA_repopush = laser_FA_repopush + 1;
            
        elseif trialTypes{ses}(i) == 4 && moveType{ses}(i) == 10 % Laser FA & simple push
            laser_FA_simpush = laser_FA_simpush + 1;
            
        elseif trialTypes{ses}(i) == 5 && moveType{ses}(i) == 9 % Control hits & reposition + push
            ctrl_hit_repopush = ctrl_hit_repopush + 1;
            
        elseif trialTypes{ses}(i) == 5 && moveType{ses}(i) == 10 % Control hits & simple push
            ctrl_hit_simpush = ctrl_hit_simpush + 1;
            
            %         elseif trialTypesALL(i) == 6 && trialOutcomesALL(i) == 9 % Laser hits
            %             trials_hit = [trials_hit toneTimesALL(i)];
            %             laser_hit = [laser_hit toneTimesALL(i)];
            %
            %         elseif trialTypesALL(i) == 6 && trialOutcomesALL(i) == 10 % Laser hits
            %             trials_miss = [trials_miss toneTimesALL(i)];
            %             laser_miss = [laser_miss toneTimesALL(i)];
            
        elseif trialTypes{ses}(i) == 7 && moveType{ses}(i) == 11 % Control miss & still
            ctrl_miss_still = ctrl_miss_still + 1;
            
        elseif trialTypes{ses}(i) == 7 && moveType{ses}(i) == 12 % Control miss & lengthwise
            ctrl_miss_length = ctrl_miss_length + 1;
            
        elseif trialTypes{ses}(i) == 7 && moveType{ses}(i) == 13 % Control miss & reposition
            ctrl_miss_repos = ctrl_miss_repos + 1;
            
            %         elseif trialTypesALL(i) == 8 && trialOutcomesALL(i) == 9 % Laser miss
            %             trials_CR = [trials_CR toneTimesALL(i)];
            %             laser_CR = [laser_CR toneTimesALL(i)];
            %
            %         elseif trialTypesALL(i) == 8 && trialOutcomesALL(i) == 10 % Laser miss
            %             trials_FA = [trials_FA toneTimesALL(i)];
            %             laser_FA = [laser_FA toneTimesALL(i)];
        end
    end
end

total_ctrl_CR = ctrl_CR_still + ctrl_CR_length + ctrl_CR_repos;
prop_ctrl_CR_still = ctrl_CR_still/total_ctrl_CR;
prop_ctrl_CR_length = ctrl_CR_length/total_ctrl_CR;
prop_ctrl_CR_repos = ctrl_CR_repos/total_ctrl_CR;
prop_ctrl_CR = [prop_ctrl_CR_still prop_ctrl_CR_length prop_ctrl_CR_repos];

total_laser_CR = laser_CR_still + laser_CR_length + laser_CR_repos;
prop_laser_CR_still = laser_CR_still/total_laser_CR;
prop_laser_CR_length = laser_CR_length/total_laser_CR;
prop_laser_CR_repos = laser_CR_repos/total_laser_CR;
prop_laser_CR = [prop_laser_CR_still prop_laser_CR_length prop_laser_CR_repos];

total_ctrl_FA = ctrl_FA_repopush + ctrl_FA_simpush;
prop_ctrl_FA_repopush = ctrl_FA_repopush/total_ctrl_FA;
prop_ctrl_FA_simpush = ctrl_FA_simpush/total_ctrl_FA;
prop_ctrl_FA = [prop_ctrl_FA_repopush prop_ctrl_FA_simpush];

total_laser_FA = laser_FA_repopush + laser_FA_simpush;
prop_laser_FA_repopush = laser_FA_repopush/total_laser_FA;
prop_laser_FA_simpush = laser_FA_simpush/total_laser_FA;
prop_laser_FA = [prop_laser_FA_repopush prop_laser_FA_simpush];

total_ctrl_hit = ctrl_hit_repopush + ctrl_hit_simpush;
prop_ctrl_hit_repopush = ctrl_hit_repopush/total_ctrl_hit;
prop_ctrl_hit_simpush = ctrl_hit_simpush/total_ctrl_hit;
prop_ctrl_hit = [prop_ctrl_hit_repopush prop_ctrl_hit_simpush];

%total_laser_hit =

total_ctrl_miss = ctrl_miss_repos + ctrl_miss_length + ctrl_miss_still;
prop_ctrl_miss_repos = ctrl_miss_repos/total_ctrl_miss;
prop_ctrl_miss_length = ctrl_miss_length/total_ctrl_miss;
prop_ctrl_miss_still = ctrl_miss_still/total_ctrl_miss;
prop_ctrl_miss = [prop_ctrl_miss_still prop_ctrl_miss_length prop_ctrl_miss_repos];

%total_laser_miss =

prop_CR = [prop_ctrl_CR; prop_laser_CR]';
prop_FA = [prop_ctrl_FA; prop_laser_FA]';
%prop_hit = [prop_ctrl_hit; prop_laser_hit]';
%prop_miss = [prop_ctrl_miss; prop_laser_miss]';


%% Plots
% Watch out -- category names and bars don't seem to align!!!!
push_categories = categorical({'Reposition','Clean'});
nonpush_categories = categorical({'Reposition','Lengthwise','Still'});

figure
b = bar(nonpush_categories, prop_CR);
b(1).FaceColor = [0.8500 0.3250 0.0980];
b(2).FaceColor = [0.3010 0.7450 0.9330];
ylim([0 1])
title('Correct rejections')

figure
b = bar(push_categories, prop_FA);
b(1).FaceColor = [0.8500 0.3250 0.0980];
b(2).FaceColor = [0.3010 0.7450 0.9330];
ylim([0 1])
title('False alarms')

figure
b = bar(push_categories, prop_ctrl_hit);
b(1).FaceColor = [0.8500 0.3250 0.0980];
%b(2).FaceColor = [0.3010 0.7450 0.9330];
ylim([0 1])
title('Hits')

%figure
%bar(prop_miss)
%ylim([0 1])





