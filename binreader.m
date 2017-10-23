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

% extract parameters from filename
tmp_str = strsplit(FileName, '_');

% save parameters in strings
acquis_date = tmp_str{1, 1}; acquis_time = tmp_str{1, 2}; exp_type    = tmp_str{1, 3}; fish_num    = tmp_str{1, 4}; trial_num   = tmp_str{1, 6};

% acquisition parameters
num_data_categories = 6;
camscale_px_per_mm  = 20.6; % px/mm
datarate_Hz         = 750;  % Hertz

% Read and reorganize the bin file
h        = fopen([PathName, FileName]);
tmp_data = fread(h, inf, 'float');
fclose(h);

tmp_data = tmp_data(1:(end-mod(size(tmp_data,1), num_data_categories)), 1); % cuts when parts of the data categories are missing at the end
tmp_data = reshape(tmp_data, [num_data_categories, size(tmp_data, 1)/num_data_categories])';

%%
% CHECK FOR TIMING PROBLEMS AND LOST FRAMES

% TIMER COUNTER: 
% time difference between frames in microseconds, based on the cameras 32bit time stamp counter (25Mhz, 40ns)

time_diff      = [0; diff(tmp_data(:, 2))];                      

% linearize the "saw-tooth function" for the timecounter
idx            = time_diff <= -2^32/1000 + 1.5*median(time_diff); % find the frames when 32bit timecounter was reset to zero
time_diff(idx) = median(time_diff);                               % reset the time difference between these frames to the median in microseconds

% camera frame length in microseconds calculated from the camera timecounter
frame_length_calc_ms = median(time_diff);

% aquisition datarate in Hertz calculated from the camera timecounter
datarate_calc_Hz = (1/(median(time_diff)).*10^6);  

% CHECK for timing problems (e.g. frames that took unusually long or are
% unusually shorter than what is expected from the used datarate)

idx_time       = abs(time_diff) >= 1.01*frame_length_calc_ms;  % searches for time differences between frames that are +-1% of the expected frame duration


% FRAME COUNTER:
% index difference between frames, based on the cameras 24bit frame counter

frame_diff = [0; diff(tmp_data(:, 1))]; 

% linearize the "saw-tooth function" for the frame counter (should not
% happen at low datarates) 
idx             = frame_diff <= -2^24 + 1.5*median(frame_diff); % find the frames when 24bit framecounter was reset to zero
frame_diff(idx) = 1;                                            % reset the frame difference between these frames to 1

% CHECK for missing frames
idx_frame  = frame_diff > 1;                         % index of missing frames
idx_lost   = find(idx_frame == 1);                   % first frame in the block of missed frames

% checks if missed timestamps coincide with missed frames, is 0 if inconsistent timestamps outside of missed frames 
isTime = isequal(idx_time, idx_frame);

% prints the above calculated values

fprintf('\nacquisition time: %s %s', acquis_date, acquis_time); 
fprintf('\nexperiment type: %s', exp_type); 
fprintf('\nfish number: %s', fish_num); 
fprintf('\ntrial number: %s', trial_num);
fprintf('\n\nfirst frame in the block of missed frames : number of frames lost\n');
fprintf('\n %d: %d',  [idx_lost, frame_diff(idx_frame)-1].');
fprintf('\n\ntiming flawed (outside of lost frames):  %d  \n', ~isTime  );


%%
% INSERT nans for lost frames...

% define anonymous function that inserts (nxm) blocks into (oxm) matrices
insert_blocks = @(insert_block, matrix, n) cat(1,  matrix(1:n-1,:), insert_block, matrix(n:end,:) );

data_raw = tmp_data;

for ii = nnz(idx_frame):-1:1 % starts from the last row in the matrix to keep the indizes for block-insertion
 
    nan_block       = nan(frame_diff(idx_lost(ii)) - 1, num_data_categories);
    nan_block(:, 1) = tmp_data(idx_lost(ii)-1, 1)+1: tmp_data(idx_lost(ii)-1, 1) + frame_diff(idx_lost(ii))-1; % fill the first column of the Nan blocks with frame numbers that were missing
    
    tmp_data        = insert_blocks(nan_block, tmp_data, idx_lost(ii));
    
end

tmp_data(:,1) = tmp_data(:,1) - tmp_data(1,1) + 1; % framecounter starts at 1

%%

%-------------------------------------------------------------------------------------
%
% IDENTIFICATION OF SWIM BOUTS (bout speeds, IBIs etc)
%
%-------------------------------------------------------------------------------------

xpos = tmp_data(:, 4);
ypos = tmp_data(:, 5);

% plot fish position 
xx = [xpos'; xpos']; yy = [ypos'; ypos'];
z  = 1:size(xpos', 2); zz = [z; z];       % create frame-vector

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
tmp_dist_unfilt           = tmp_dist_unfilt./camscale_px_per_mm;  % convert to distance in mm
tmp_vel_unfilt            = tmp_dist_unfilt.*datarate_Hz;         % convert to velocity in mm/s

%%

idx_nan     = isnan(dx);
dx(idx_nan) = 0; % for filtering nan values need to be removed 
dy(idx_nan) = 0;

filter     = ones(1, 55)./55; % filter tbd

dxf        = filtfilt(filter, 1, dx);
dyf        = filtfilt(filter, 1, dy);

dxff       = fbtrim(dx, [0.5, 1, 0.1, 15, 5]); % [0.5, 0.6, 0.01, 5, 5]

tmp_fdist = sqrt(dxf.^2 + dyf.^2);     % distance moved between iterations, in pixels
tmp_fdist = tmp_fdist./camscale_px_per_mm;  % convert to distance in mm
tmp_fvel  = tmp_fdist.*datarate_Hz;  % convert to velocity in mm/s

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

%%
fig1 = figure;
hold on; 
plot(tmp_vel_unfilt); 
plot(tmp_fvel,'LineWidth', 3);

plot(dx-10);
plot(20*dxf-10);
%plot(20*dxff-10);

plot(20*dyf-20);
plot(dy-20);
plot(50*tmp_delta_ori-40);
hold off;

dcm1 = datacursormode(fig1);
set(dcm1, 'UpdateFcn', @Data_Cursor_precision, 'Enable', 'on');

