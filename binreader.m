%to read bin files generated from MikroBirdDrishti_Grab&Stim
%by Gokul Rajan, DEL-BENE Lab, 2017-08-10.
%added more parameters on 2017-08-11.

IFTms = 1000/949; % change this if acq was not at 500 Hz
h = fopen('8_11_2017_1.bin');
test = fread(h,inf,'float');
NumVar = 7;

frame = test(1:NumVar:end-1);
frameNorm = frame-(frame(1)-1);

timeMillis = frame*IFTms;
timeNormMillis = timeMillis-(timeMillis(1)-1);
millisDiff = diff(timeMillis);
totalMillis = timeMillis(end)-timeMillis(1);

timeSec = ((frame*IFTms)/1000);
totalSec = timeSec(end)-timeSec(1);
timeNormSecs = timeSec - (timeSec(1)-1);
totalMin = totalSec/60;

posX = test(3:NumVar:end);
posXfilt=filtfilt(ones(1,5)/5,1,posX);

posY = test(4:NumVar:end);
posYfilt=filtfilt(ones(1,5)/5,1,posY);

angle = test(5:NumVar:end);
%angle = rad2deg(test(5:4:end));

ledStatus = test(6:NumVar:end);

ledRaw = test(7:NumVar:end);

CamTime = test(2:NumVar:end);

angleDiff = (diff(angle));
angleDiff(end+1) = NaN;
blank = nan(size(posX));
Reye = blank;
Leye = blank;
loomX = blank;
loomY = blank;

INST_VEL = (diff(posXfilt).^2+diff(posY).^2).^0.5;
            
mydata = [angle,frame,timeMillis,posX,posY,ledStatus,ledRaw,CamTime,loomX,loomY,Leye,Reye,angleDiff];
csvwrite('/Institut Curie/Lab/Projects/Reelin/Behavior/preTests/2017-08-07_FS01/trial_AN100.csv',mydata);
%comet3(posX,posY,timeMillis);
%plot(posX,posY);
figure(1); plot(posXfilt,posYfilt);

figure(2); plot(ledStatus*100,'r');
hold on; plot(INST_VEL,'b');
