%% MikroBirdDrishti Zebranalysis Script
% by Gokul Rajan, DEL-BENE Lab, Paris - Aug 2017

% adapted from Bonsai AQ55/A56/A60 analysis script written by Isaac Bianco, Aug 2015
% and modified by Chris, Oct 2015 - July 2016

%% some inits/ change here
clearvars; clc;

%testfilename = '8_11_2017_S1.bin';
csvfilename = '8_11_2017_S1.csv';

IFTms = 1000/949; % change this if acq freq changes
IFTs = IFTms/1000;
px_per_mm = 20.6;

root = '/Institut Curie/Lab/Projects/Reelin/Behavior/Swim_Kinematics/2017_08_16_RLN/';
foldername = fullfile(root,'SwimTap');
allFileNames =dir(fullfile(foldername,'*AN*.bin'));

%num_fish             = sum(structfun(@numel, free_swimming));
zebranalysis_folders = struct2cell( dir(root) )';


%% Read csv file

% angle,
% frame,timeMillis,
% posX,posY,
% ledStatus,ledRaw,
% CamTime,
% loomX,loomY,
% Leye,Reye,
% angleDiff

fid = fopen('8_11_2017_S1.csv');
tmp_data = textscan(fid,'%n%n%n%n%n%n%n%n%n%n%n%n%n', 'delimiter', ',');
fclose(fid);

tmp_data = cell2mat(tmp_data);

%% correction for delta angle

%tmp_angle = tmp_data(:,1);
tmp_delta_ori_raw = tmp_data(:,13);% CAVEAT: also remaps the values such that the turn is always between the shortest delta orientation

% correction of discontinuties when the fish turns from 0 to 360 or vice-versa
  for kk = 1: length(tmp_delta_ori_raw)
                
      if (abs(tmp_delta_ori_raw(kk)) >= 250) && (tmp_delta_ori_raw(kk) < 0)  % left turn
         tmp_delta_ori_raw(kk) =  tmp_delta_ori_raw(kk) + 360;
                    
      elseif (abs(tmp_delta_ori_raw(kk)) >= 250) && (tmp_delta_ori_raw(kk) > 0) % right turn
         tmp_delta_ori_raw(kk) = tmp_delta_ori_raw(kk) - 360;
      end
                
  end
                  
% save delta-orientation derived from the raw, not interpolated values
  tmp_data(:, 13) = tmp_delta_ori_raw;
  

%%  No interpolation

        % extract time vector
        tmp_time = (tmp_data(:,3)); % time in ms
        
%         %LINEAR INTERPOLATION of data points for the frames that were dropped by the camera 
% 
%         % fixed new timebase
%         % tmp_timebase = 0:(1/free_swimming.(identifier).brightfield.(trial_field).datarate_Hz):tmp_t_vect(end); 
% 
%         for ll = 1:10 % 10 columns of data from csv file
% 
%            tmp_datainterp(:, ll) = interp1(tmp_t_vect, tmp_data(:, ll), tmp_timebase);
% 
%         end
  
 %% Copy all data for cals 
 %No selection of dish center as data already cleaned for missing points 
 %in the binreader script
 
 tmp_posX = tmp_data(:,4); %posX
 tmp_posY = tmp_data(:,5); %posY
 tmp_ori = tmp_data(:,1); %angle
 tmp_deltaOri = tmp_data(:,13);%deltaAngle
 
 
%% free swimming param cals

dx = [0; diff(tmp_posX)]; % distance between two consecutive x-coordinates
dy = [0; diff(tmp_posY)]; % distance between two consecutive y-coordinates
tmp_dist_unfilt =  sqrt(dx.^2 + dy.^2);

tmp_dist_unfilt = tmp_dist_unfilt/px_per_mm;  % convert to distance in mm
tmp_vel_unfilt  = tmp_dist_unfilt./tmp_time;  % convert to velocity in mm/ms

    filterB = ones(1,10)./10; % broad filter for event detection, better detection of bout_off
    filterV = ones(1, 3)./ 3; % narrow filter for peak velocity/timing detection, better detection of bout_on
     
            dxB                       = filtfilt(filterB, 1, dx);
            dyB                       = filtfilt(filterB, 1, dy);
            tmp_fdistB                = sqrt(dxB.^2 + dyB.^2);     % distance moved between iterations, in pixels               
            tmp_fdistB                = tmp_fdistB./px_per_mm;  % convert to distance in mm
            tmp_fvelB                 = tmp_fdistB./tmp_time;  % convert to velocity in mm/s           
            
            dxV                       = filtfilt(filterV, 1, dx);
            dyV                       = filtfilt(filterV, 1, dy);
            tmp_fdistV                = sqrt(dxV.^2 + dyV.^2); % distance moved between iterations, in pixels
            tmp_fdistV                = tmp_fdistV./px_per_mm; % convert to distance in mm
            tmp_fvelV                 = tmp_fdistV./tmp_time;  % convert to velocity in mm/s
            
figure(1);plot(tmp_vel_unfilt,'k'); hold on; plot(tmp_fvelB,'b'); hold on; plot(tmp_fvelV,'r');
figure(2); polarhistogram(tmp_ori,30);

% inst velocity
%INST_VEL = ((diff(posX).^2+diff(posY).^2).^0.5)/(millisDiff);

% turn classification
DELTA_ORI = tmp_deltaOri;
left = (DELTA_ORI > 0); % assuming no turns >110 degree betn two frames
right = (DELTA_ORI < 0);
TR= sum(right);
TL = sum(left);
rightProp = TR/(TL+TR);
leftProp = -(TL/(TL+TR));
turns = [rightProp,leftProp];
figure(3); barh(turns);
