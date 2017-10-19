%% Figures on/off

set(0, 'DefaultFigureVisible', 'on');

%%

% Data key
% 1 frame #
% 2 camera time in microsec (for data transfer purposes microseconds instead of nanoseconds used)
% 3 system time in ms
% 4 xpos, start at 0
% 5 ypos
% 6 angle in rad

clearvars; clc;

[FileName, PathName] = uigetfile('*.bin');
cd(PathName);

num_data_categories = 6;
camscale            = 20.6; % px/mm
datarate_acquis     = 750;  % Hertz

% Read and reorganize the bin file

h        = fopen([PathName, FileName]);
tmp_data = fread(h, inf, 'float');
fclose(h);

tmp_data = tmp_data(1:(end-mod(size(tmp_data,1), num_data_categories)), 1); % cuts when parts of the data categories are missing at the end
tmp_data = reshape(tmp_data, [num_data_categories, size(tmp_data, 1)/num_data_categories])';

% camera timecounter and datarate
% linearize the saw-tooth function timestamps

time_diff      = [0; diff(tmp_data(:, 2))];                     % in microseconds
idx            = time_diff <= -2^32/1000 + 2*median(time_diff); % find the frames when timecounter was reset to zero
time_diff(idx) = median(time_diff);                             % reset the time difference between these frames to the median in microseconds

% camera frame length in microseconds based on the camera timestamps
frame_length   = median(time_diff);

% aquisition datarate in Hertz based on the camera timestamps
datarate_calc  = (1/(median(time_diff)).*10^6);  

idx_time       = abs(time_diff) >= 1.1*frame_length; % searches for time differences between frames that are +-10% of the expected frame duration

% camera frame counters and missing frames
% fill in nan values when frames were lost (most likely due to fish being lost during tracking ) 

frame_diff = [0; diff(tmp_data(:, 1))]; 
idx_frame  = frame_diff > 1;                             % index of missing frames

idx_lost        = find(idx_frame == 1); % first frame in the block of missed frames
num_frames_lost = sum(frame_diff(idx_frame) - 1);

fprintf('first frame in the block of missed frames : number of frames lost\n\n');
fprintf('%d: %d \n',  [idx_lost, frame_diff(idx_frame)-1].');

% checks if missed timestamps coincide with missed frames, is 0 if inconsistent timestamps outside of missed frames 
isTime = isequal(idx_time, idx_frame);

if isTime ==0
    fprintf('timing counter shows flawed interframe timing independent of missing frames!!!');
end

% insert nans for lost frames...

insert_blocks = @(insert_block, vector, n) cat(1,  vector(1:n-1,:), insert_block, vector(n:end,:) );

raw_data = tmp_data;

for ii = nnz(idx_frame):-1:1
 
    nan_block       = nan(frame_diff(idx_lost(ii)) - 1, num_data_categories);
    nan_block(:, 1) = tmp_data(idx_lost(ii)-1, 1)+1: tmp_data(idx_lost(ii)-1, 1) + frame_diff(idx_lost(ii))-1; % fill the first column of the Nan blocks with frame numbers that were missing
    
    tmp_data = insert_blocks(nan_block, tmp_data, idx_lost(ii));
    
end

tmp_data(:,1) = tmp_data(:,1) - tmp_data(1,1) + 1; % framecounter starts at 1

%-------------------------------------------------------------------------------------
%
% IDENTIFICATION OF SWIM BOUTS (bout speeds, IBIs etc)
%
%-------------------------------------------------------------------------------------

xpos = tmp_data(:, 4);
ypos = tmp_data(:, 5);

% plot fish position 
xx = [xpos'; xpos']; yy = [ypos'; ypos'];
z  = 1:size(xpos', 2); zz = [z; z]; % frame-vector

figure
hs = surf(xx, yy, zz, zz, 'EdgeColor', 'interp', 'LineWidth', 2);
colormap('parula');
view([0 90 0]); 
xlabel('x-position');
ylabel('y-position');
zlabel('frames');

dx   = [0; diff(xpos)]; % distance between two consecutive x-coordinates
dy   = [0; diff(ypos)]; % distance between two consecutive y-coordinates

tmp_dist_unfilt           = sqrt(dx.^2 + dy.^2);
tmp_dist_unfilt           = tmp_dist_unfilt./camscale;  % convert to distance in mm
tmp_vel_unfilt            = tmp_dist_unfilt.*datarate_acquis;  % convert to velocity in mm/s

idx_nan     = isnan(dx);
dx(idx_nan) = 0; % for filtering nan values need to be removed 
dy(idx_nan) = 0;

filter    = ones(1, 101)./101; % filter tbd

dxf        = filtfilt(filter, 1, dx);
dyf        = filtfilt(filter, 1, dy);

tmp_fdist = sqrt(dxf.^2 + dyf.^2);     % distance moved between iterations, in pixels
tmp_fdist = tmp_fdist./camscale;  % convert to distance in mm
tmp_fvel  = tmp_fdist.*datarate_acquis;  % convert to velocity in mm/s

tmp_fvel(idx_nan) = nan; % re-insert the nan values

% orientation

tmp_delta_ori = [nan; diff(tmp_data(:, 6))];

% correction of discontinuties when fish turns from 0 to 2*pi or vice versa

for kk = 1: length(tmp_delta_ori)
    
    if tmp_delta_ori(kk) > pi % right turn
        tmp_delta_ori(kk) =  tmp_delta_ori(kk) - 2*pi;
        
    elseif tmp_delta_ori(kk) < -pi % left turn
        tmp_delta_ori(kk) =  2*pi + tmp_delta_ori(kk);
    end
    
end   

figure, 
hold on; 
plot(tmp_vel_unfilt); 
%plot(tmp_fvelB,'LineWidth',5); 
%plot(tmp_data(:,6));
plot(50*tmp_delta_ori-30);
hold off;

