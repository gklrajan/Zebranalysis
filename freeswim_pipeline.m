%% MikroBirdDrishti Zebranalysis Script
% by Gokul Rajan, DEL-BENE Lab, Paris - Aug 2017

% adapted from Bonsai AQ55/A56/A60 analysis script written by Isaac Bianco, Aug 2015
% and modified by Chris, Oct 2015 - July 2016

%% some inits/ change here
clearvars; clc;

testfilename = '8_16_2017_AN012.bin';
csvfilename = '8_16_2017_AN012.csv';

root = '/Institut Curie/Lab/Projects/Reelin/Behavior/Swim_Kinematics/2017_08_16_RLN/';
foldername = fullfile(root,'SwimTap');
allFileNames =dir(fullfile(foldername,'*AN*.bin'));

%num_fish             = sum(structfun(@numel, free_swimming));
zebranalysis_folders = struct2cell( dir(root) )';


%% Read csv file



%% fish param cals

% inst velocity
INST_VEL = ((diff(posXfilt).^2+diff(posY).^2).^0.5)/(millisDiff);

fieldname = 'turn';

% delta orientation
DELTA_ORI = diff(angle);

% turn classification
if (DELTA_ORI > 0) | ((abs(DELTA_ORI) >= 250) & (DELTA_ORI < 0)) % assuming no turns >110 degree betn two frames
    p.(fieldname)= -1; %left
    
elseif (DELTA_ORI < 0)|((abs(DELTA_ORI) >= 250) & (DELTA_ORI > 0))
    p.(fieldname)= 1; %right
    
end

%% pooling



%% mean cals