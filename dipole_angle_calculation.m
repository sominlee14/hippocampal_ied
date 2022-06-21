% Script for calculating the dipole angle 
% Script assumes that scalp STAs have been saved by group with the naming
% convention: 'pt_id_scalp_sta_group_x'

% load channel locations
load('channel_locs.mat')
channel_locs = struct2table(channel_locs);

epoch_length = 512;
sample_rate = 1024;

% Find all groups for patient scalp STAs
pt_dir = struct2table(dir('pt_id_scalp_sta.mat'));
pt_filenames = pt_dir.name;
n_groups = length(pt_filenames)
    
time_points = 1:512;
dipole_threshold = 3.5;

% Allocate for various results
max_channel = nan(n_groups,1);
max_abs_potential = nan(n_groups, 1);
i_max = nan(n_groups, 1);

sta_time_course = cell(n_groups, 1);
dipole_diff = nan(length(time_points), n_groups);

for i = 1:n_groups
    
    % load scalp STA for given group
    load(pt_filenames{i})
    figure; hold('on')
    title(['Group ', num2str(i)])
    offset = 5;
    x = (1:epoch_length)/sample_rate*1000;
    
    % Plot to look at scalp STAs
    plot_offset_signals2(x, sta(:, 1:epoch_length), offset, 'k')
    plot_channel_labels(channel_locs.labels, offset)
    
    fprintf(' Group %d results', i)

    % Find channel and time of max amplitude to determine t = 0
    [ch_peak, i_ch_peak] = max(abs(sta(:, 1:epoch_length)), [], 2);
    
    [~, max_channel(i)] = max(ch_peak);
    
    max_channel_label = channel_locs.labels(max_channel(i))
    m = i_ch_peak(max_channel(i))
    i_max(i) = m;
    
    m_pot = sta(max_channel(i), i_max(i))
    
    max_abs_potential(i) = m_pot;
    t = time_points;
    sta_time_course{i} = sta(:, t);
    
end


significant_dipole_angles = cell(n_groups,1);
std_dipole_angles = nan(n_groups,1);

for i = 1:n_groups
    
    % calculate dipole angles for each time point using scalp coordinates 
    time_course = sta_time_course{i};
    max_values = nan(size(time_course,2),1);
    min_values = nan(size(time_course,2),1);
    max_channel = nan(size(max_values));
    min_channel = nan(size(min_values));
    dipole_angle = nan(size(time_course,2),1);

    for j = 1:size(time_course,2)
        
        [max_values(j), max_channel(j)] = max(time_course(:,j));
        [min_values(j), min_channel(j)] = min(time_course(:,j));

        x_max = channel_locs.X(max_channel(j));
        y_max = channel_locs.Y(max_channel(j));
        x_min = channel_locs.X(min_channel(j));
        y_min = channel_locs.Y(min_channel(j));
        
        x_diff = x_max - x_min;
        y_diff = y_max - y_min;
        angle = atan(abs(y_diff)/abs(x_diff)) * 180 / pi;

        dipole_angle(j) = angle;
        
    end
    
    dipole_diff(:,i) = max_values - min_values;
    
    % Include only time points where there is a significant dipole 
    significant_dipole_angles{i} = dipole_angle(dipole_diff(:,i) >= dipole_threshold);

    % standrad dev of significant dipole angles 
    std_dipole_angles(i) = std(significant_dipole_angles{i});
end
















