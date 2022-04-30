% Small script to plot a histogram of the depths of different morphological
% features of corticospinal neurons.
% Author: Victor Chamosa Pino
%

% Load data
if exist('\\ardbeg.mvm.ed.ac.uk\duguidlab\motor_choice')
    datapath = '\\ardbeg.mvm.ed.ac.uk\duguidlab\motor_choice\mapping\data\interim\CSNs\measurements\';
else
    datapath = 'D:\ALM_Ardbeg\Axioscan\measurements\';
end

depth_arbor = fullfile(datapath,'depth_arbor.mat');
depth_arbor = load(depth_arbor);
depth_arbor = depth_arbor.depth_arbor;

depth_nexus = fullfile(datapath,'depth_nexus.mat');
depth_nexus = load(depth_nexus);
depth_nexus = depth_nexus.depth_nexus;

depth_soma = fullfile(datapath,'depth_soma.mat');
depth_soma = load(depth_soma);
depth_soma = depth_soma.depth_soma;


%% Plots

% Plot raw
figure
hold on
histogram(depth_arbor,'BinWidth',50,'Orientation','horizontal')
histogram(depth_nexus,'BinWidth',50,'Orientation','horizontal')
histogram(depth_soma,'BinWidth',50,'Orientation','horizontal')
hold off
set(gca,'YDir','reverse')
ylim([0 1050])
xlim([0 18])
ylabel('Depth (um)')

% Plot normalised
figure
hold on
histogram(depth_arbor,'BinWidth',50,'Orientation','horizontal','Normalization','probability')
histogram(depth_nexus,'BinWidth',50,'Orientation','horizontal','Normalization','probability')
histogram(depth_soma,'BinWidth',50,'Orientation','horizontal','Normalization','probability')
hold off
set(gca,'YDir','reverse')
ylim([0 1050])
xlim([0 1])
ylabel('Depth (um)')




% old version of plot
% 
% % Calculate bin sizes
% [arbor_s, arbor_l] = bounds(depth_arbor);
% arbor_dist = arbor_l - arbor_s;
% arbor_bins = round(arbor_dist/50); % + 2;
% 
% [nexus_s, nexus_l] = bounds(depth_nexus);
% nexus_dist = nexus_l - nexus_s;
% nexus_bins = round(nexus_dist/50); % + 2;
% 
% [soma_s, soma_l] = bounds(depth_soma);
% soma_dist = soma_l - soma_s;
% soma_bins = round(soma_dist/50); % + 2;
%
% % figure 
% hold on
% histogram(depth_arbor,arbor_bins,'Orientation','horizontal')
% histogram(depth_nexus,nexus_bins,'Orientation','horizontal')
% histogram(depth_soma,soma_bins,'Orientation','horizontal')
% hold off
% set(gca, 'YDir','reverse')
% ylabel('Depth (um)')

