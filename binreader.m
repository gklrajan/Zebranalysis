
clearvars; clc;

fname = uigetfile('*.bin');

% data key
% 1 frame #
% 2 camera time in microsec (for data transfer purposes microseconds instead of nanoseconds used)
% 3 system time in ms
% 4 xpos, start at 0
% 5 ypos
% 6 angle in rad

num_data_categories = 6;

% Read and reorganize the bin file

h        = fopen(fname);
tmp_data = fread(h, inf, 'float');
fclose(h);

tmp_data = tmp_data(1:(end-mod(size(tmp_data,1), num_data_categories)), 1); % cuts when parts of the data categories are missing at the end
tmp_data = reshape(tmp_data, [num_data_categories, size(tmp_data, 1)/num_data_categories])';

%%

% camera timecounter and datarate
% linearize the saw-tooth function timestamps

time_diff           = [0; diff(tmp_data(:, 2))];                     % in microseconds
idx_time            = time_diff <= -2^32/1000 + 2*median(time_diff); % find the frames when timecounter was reset to zero
time_diff(idx_time) = median(time_diff);                             % reset the time difference between these frames to the median in microseconds
    
timebase = (cumsum(time_diff) + 1) ./10^6;                           % new timebase in seconds

% camera frame length in microsecond
frame_length = median(time_diff);

%
datarate     = round(1/(median(time_diff)).*10^6);                   % aquisition datarate in Hertz based on the frame timestamps

% camera frame counters and missing frames
% fill in nan values when frames were lost (most likely due to fish being lost during tracking ) 

frame_diff = [0; diff(tmp_data(:, 1))]; 
idx_frame  = frame_diff > 1;                             % find missing frames

num_frames_lost = sum(frame_diff(idx_frame) - 1);      
fprintf('total number of frames when fish was most likely lost during tracking: %d\n ', num_frames_lost);

% insert nans for lost frames...

insert_blocks = @(insert_block, vector, n) cat(1,  vector(1:n-1,:), insert_block, vector(n:end,:) );
idx_lost      = find(idx_frame == 1);

for ii = nnz(idx_frame):-1:1
 
    nan_block       = nan(frame_diff(idx_lost(ii)) - 1, num_data_categories);
    nan_block(:, 1) = tmp_data(idx_lost(ii)-1, 1)+1: tmp_data(idx_lost(ii)-1, 1) + frame_diff(idx_lost(ii))-1; % fill the first column of the Nan blocks with frame numbers that were missing
    
    
    tmp_data = insert_blocks(nan_block, tmp_data, idx_lost(ii));
    
end

tmp_data(:,1) = tmp_data(:,1) - tmp_data(1,1)+1;


%%



xpos = tmp_data(:, 4);
ypos = tmp_data(:, 5);

dx            = [0; diff(xpos)]; % distance between two consecutive x-coordinates
dy            = [0; diff(ypos)]; % distance between two consecutive y-coordinates

tmp_dist_unfilt           =  sqrt(dx.^2 + dy.^2);
tmp_dist_unfilt           = tmp_dist_unfilt./20;  % convert to distance in mm
tmp_vel_unfilt            = tmp_dist_unfilt.*datarate;  % convert to velocity in mm/s

filterB = ones(1, 50)./50; % broad filter for event detection, better detection of bout_off

dxB                       = filtfilt(filterB, 1, dx);
dyB                       = filtfilt(filterB, 1, dy);

tmp_fdistB                = sqrt(dxB.^2 + dyB.^2);     % distance moved between iterations, in pixels
tmp_fdistB                = tmp_fdistB./20;  % convert to distance in mm
tmp_fvelB                 = tmp_fdistB.*datarate;  % convert to velocity in mm/s



figure, plot(tmp_vel_unfilt); hold on; plot(tmp_fvelB,'LineWidth',5); hold on; plot(50*diff([0; tmp_data(:,6)])-100);

