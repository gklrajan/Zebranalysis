%to read bin files generated from MikroBirdDrishti_Grab&Stim
%by Gokul Rajan, DEL-BENE Lab, 2017-08-10.
%added more parameters on 2017-08-11.

%% some inits/ change here before run
clearvars; clc;

testfilename = '8_11_2017_S1.bin';
csvfilename = '8_11_2017_S1.csv';

IFTms = 1000/949; % change this if acq freq changes
IFTs = IFTms/1000;
NumVar = 7;

%% Read the bin file
h = fopen(testfilename);
test = fread(h,inf,'float');

frame = test(1:NumVar:end);
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

%% compute time series and other params

frameNorm = frame-(frame(1));
frameNormDiff = diff(frameNorm);
frameNormDiff = [0;frameNormDiff];
timeMillis = frameNormDiff*IFTms;
timeNormMillis = timeMillis-(timeMillis(1));
millisDiff = diff(timeNormMillis);
totalMillis = timeMillis(end)-timeMillis(1);

timeSec = ((frame*IFTms)/1000);
totalSec = timeSec(end)-timeSec(1);
timeNormSecs = timeSec - (timeSec(1));
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


%% some basic plots and writing to csv
csvwrite(csvfilename,mydata);

figure(1);plot(posXfilt,posYfilt);
figure(2);plot3(posXfilt,posYfilt,timeNormMillis);

INST_VEL = (diff(posXfilt).^2+diff(posY).^2).^0.5;
figure(3);plot(ledStatus*100,'r');
hold on;plot(INST_VEL,'b');


%% Other saving stuff
filename1 = sprintf('/Institut Curie/Lab/Projects/Reelin/Behavior/preTests/track_%s.jpg',testfilename);
filename2 = sprintf('/Institut Curie/Lab/Projects/Reelin/Behavior/preTests/trqck3d_%s.jpg',testfilename);
filename3 = sprintf('/Institut Curie/Lab/Projects/Reelin/Behavior/preTests/stim_%s.jpg',testfilename);
filename4 = sprintf('/Institut Curie/Lab/Projects/Reelin/Behavior/preTests/allVars_%s.mat',testfilename);
saveas(figure(1),filename1);
saveas(figure(2),filename2);
saveas(figure(3),filename3);
save(filename4);
