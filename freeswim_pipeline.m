%% MikroBirdDrishti Zebranalysis Script
% by Gokul Rajan, DEL-BENE Lab, Paris - Aug 2017

% adapted from Bonsai AQ55/A56/A60 analysis script written by Isaac Bianco, Aug 2015
% and modified by Chris, Oct 2015 - July 2016

%% some inits/ chage here
clearvars; clc;

testfilename = '8_16_2017_AN012.bin';
csvfilename = '8_16_2017_AN012.csv';

root = '/Institut Curie/Lab/Projects/Reelin/Behavior/Swim_Kinematics/2017_08_16_RLN/';
foldername = fullfile(root,'SwimTap');
allFileNames =dir(fullfile(foldername,'*AN*.bin'));

%num_fish             = sum(structfun(@numel, free_swimming));
zebranalysis_folders = struct2cell( dir(root) )';

%% Read the bin file
IFTms = 1000/949; % change this if acq was not at 500 Hz
h = fopen(testfilename);
test = fread(h,inf,'float');
NumVar = 7;

frame = test(1:NumVar:end-1);
posX = test(3:NumVar:end);
posY = test(4:NumVar:end);
angle = test(5:NumVar:end);
angle = rad2deg(angle);

ledStatus = test(6:NumVar:end);
ledRaw = test(7:NumVar:end);
CamTime = test(2:NumVar:end);

%% cleaning data

if (isnan(angle)) % replaces all frames where fish track was lost into NaNs. 
    posX = 'NaN';
    posY = 'NaN';
end

%% compute time series and other params and write to csv

frameNorm = frame-(frame(1)-1);

timeMillis = frame*IFTms;
timeNormMillis = timeMillis-(timeMillis(1)-1);
millisDiff = diff(timeNormMillis);
totalMillis = timeMillis(end)-timeMillis(1);

timeSec = ((frame*IFTms)/1000);
totalSec = timeSec(end)-timeSec(1);
timeNormSecs = timeSec - (timeSec(1)-1);
totalMin = totalSec/60;

posXfilt=filtfilt(ones(1,5)/5,1,posX);
posYfilt=filtfilt(ones(1,5)/5,1,posY);

angleDiff = (diff(angle));
angleDiff(end+1) = NaN;
blank = nan(size(posX));
Reye = blank;
Leye = blank;
loomX = blank;
loomY = blank;
            
mydata = [angle,frame,timeMillis,posX,posY,ledStatus,ledRaw,CamTime,loomX,loomY,Leye,Reye,angleDiff];
%dt = datestr(now,'mm-dd-yyyy_HH-MM');

csvwrite(csvfilename,mydata);
%comet3(posX,posY,timeMillis);
%plot(posX,posY);


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

%% File save stuff

figure(1); plot(posXfilt,posYfilt);
figure(2); plot(ledStatus*100,'r');
hold on; plot(INST_VEL,'b');

filename1 = sprintf('/Institut Curie/Lab/Projects/Reelin/Behavior/preTests/track_%s.jpg',testfilename);
filename2 = sprintf('/Institut Curie/Lab/Projects/Reelin/Behavior/preTests/stim_%s.jpg',testfilename);
filename3 = sprintf('/Institut Curie/Lab/Projects/Reelin/Behavior/preTests/allVars_%s.mat',testfilename);
saveas(figure(2),filename1);
saveas(figure(1),filename2);
save(filename3);