%% Analysis of go/no-go training sessions
%
% This script extracts behavioural data from go/no-go training sessions of
% mice selected by the user. For further analysis and plotting refer to
% GoNoGo_plotting.m
% Author: Victor Chamosa Pino
%

%% Select your mice
% selectanimals() will give you a cell array listing mouse IDs
% it is passed any number of key-value pairs and looks for them in animals.json

mice = selectanimals('task','GoNoGo','responsible_person','Victor');

% if you want to exclude mice you can do it here
mice_unwanted = horzcat(selectanimals('strain', 'Thy1G'),selectanimals('strain', 'Thy1GCaMP'),selectanimals('strain', 'Thy1GCaMP6s'));
mice = setdiff(mice, mice_unwanted);
for iMouse = {'VGATCre761' 'VGATCre764' 'MCos402' 'MCos403' 'MCos216' 'MCos217' 'MCos863' 'MCos864' 'MCos904' 'MCos905' 'MCos906'}
    mice = mice(~strcmp(mice, iMouse));
end
numMice = numel(mice); 

 
%% Download their data
% This only needs to be done when there is new training data
% mat files containing abf/daq file events will be downloaded locally for your mice

% it takes the list of mice, the local path to CueBehaviourAnalysis/Data/Raw,
% the path to the locally mounted ardbeg server's Raw, and a cell array of strings listing file types you want to download
syncdata(mice,'\\ardbeg.mvm.ed.ac.uk\duguidlab\motor_choice\behaviour','D:\\ALM_Ardbeg\\Behavioural_data\\Data_analysed',{'json','mat', 'daq', 'abf'});


%% Get metadata
% this gives you a structure containing the metadata contained in animals.json if you need it
metadata = getmousemetadata(mice);

tMetadata = cell(1, numMice);
for iMouse = 1 : numMice
    tMetadata{iMouse} = loadjson([mice{iMouse} '.json']);
end


%% Generate behaviour variables on training data

trainingdataanalysis(mice,'metadata',metadata,'noskip','dataDir','D:\\ALM_Ardbeg\\Behavioural_data\\Data_analysed');


%% Load the events saved in the mat files

Data = cell(1, numMice);
for iMouse = 1 : numMice
    Data{iMouse} = load([mice{iMouse} '_analysed.mat']);
end

