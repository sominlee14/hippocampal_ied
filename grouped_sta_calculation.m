% This script assumes that peak detection has already been performed 
% using peak_detection.m 

% This script performs the following analyses: 
% (1) Grouping of hippocampal interictal discharges (HC-IED) using
% connected components
% (2) RHD and scalp STA calculations by group
% (3) Cross-correlation analysis between RHD and scalp STA


% Load peaks and metadata for patient
% prior to this analysis, all hippocampal discharges were detected and
% saved into a file with the naming convention: 'pt_id_signal_clips.mat'
% metadata (channel list and sample rate) stored in mat file
% 'pt_id_metadata.mat'

addpath('signal_proc_functions/')

% Set pt id to load relevant data
pt_file_ID = 'pt_id'; 

load([pt_file_ID,'_metadata.mat'])
load([pt_file_ID,'_signal_clips.mat'])
%%
mat_factor = 250; % Conversion into uV
signal_clips = signal_clips / mat_factor; 
%%
n_clips = size(signal_clips, 3);
clip_length = size(signal_clips, 2);
clip_trigger = clip_length/2 + 1; % Peak location (i.e. center of clip)

% Find relevant channels (RHD & scalp EEG)
rhd_channels = find(ismember(channel_list.channel_type, 'RHD'));
n_rhd_channels = length(rhd_channels);

scalp_channels = find(ismember(channel_list.channel_type, 'scalp_eeg'));
n_scalp_channels = length(scalp_channels);
scalp_channel_labels = channel_list.channel_label(scalp_channels); % Names of scalp EEG channels
%%

% **** ANALYSIS PARAMETERS ******

filt_bp_freq = [10,50] % Bandpass filter frequency range
filt_order = 2
filt_rep = 0 % How many additional times to run through filter (alternatively, just increase filter order)

% Get channel # and labels for scalp EEG channels used in analysis
scalp_exclude = {'FP1', 'FP2', 'M1', 'M2'}; % Scalp EEG channels not included in analysis
excluded_channels = find(ismember(scalp_channel_labels, scalp_exclude));
included_channels = find(~ismember(scalp_channel_labels, scalp_exclude));
included_channel_labels = scalp_channel_labels(included_channels);

% Length of epoch to use for analysis
epoch_length = 256; % n samples for epoch
epoch_range = ...
    (clip_trigger-(epoch_length/2)):...
    (clip_trigger+(epoch_length/2)-1); % indices to grab epoch from larger clip

% Define significant SNR 
sig_snr_dB = 5
% *********************************
% Filter signals using Butterworth bandpass filter 
% Function uses filtfilt to avoid phase shift 

rhd_filt = filter_butter(...
    signal_clips(rhd_channels,:,:), ...
    filt_bp_freq, sample_rate, filt_order);

scalp_filt = filter_butter(...
    signal_clips(scalp_channels,:,:),...
    filt_bp_freq, sample_rate, filt_order);

% Run through filter more times if needed (alternatively, can use higher
% order for Butterworth filter)

if filt_rep > 0
    
    for i = 1:filt_rep
        rhd_filt = filter_butter(rhd_filt, filt_bp_freq, sample_rate, filt_order);
        scalp_filt = filter_butter(scalp_filt, filt_bp_freq, sample_rate, filt_order);
    end
    
end

% Make average reference montage 
% Both RHD and scalp signals referenced to average of scalp signals
% FP1, FP2, M1, M2 channels not included

rhd_avg_ref_mont = montage_average_ref(...
    rhd_filt,...
    scalp_filt(included_channels,:,:));

scalp_avg_ref_mont = montage_average_ref(...
    scalp_filt,...
    scalp_filt(included_channels,:,:));

% ***** Classify peaks into groups based on peak morphology *****

% Concatenate RHD channels into one vector for each epoch
rhd_epochs = rhd_filt(:, epoch_range,:);
rhd_epochs_concat = string_channels(rhd_epochs);

% Correlate concatenated RHD signals with each other
rhd_epochs_concat_coeffs = corrcoef(rhd_epochs_concat');

% Determine threshold correlation coefficient to find connected components
% Find connected components using a range of thresholds and choose
% threshold that results in the greatest number of "large" groups 

thresholds_range = 0.3:0.01:1; % range of thresholds to be tested
n_group_min = 20; % minimum number of peaks to constitute a group

% Test # of groups formed for each threshold 
threshold_test = test_conn_comp_thresholds(rhd_epochs_concat_coeffs, thresholds_range, n_group_min);

% Plot results of threshold test 
all_figures(1) = figure;
subplot(2,1,1); hold('on')
plot(thresholds_range, threshold_test.n_groups, '-ko')
title(['# of bins with minimum ', num2str(n_group_min), ' peaks'])
%%
% ***** Form groups of peaks using connected components *****

% Threshold is corr coeff that results in the most number of groups 
threshold = 0.56 % determined by threshold testing code above 

[bin_assignment, n_bins, bin_sizes, conncomp_graph] = ...
    bin_by_conn_components(rhd_epochs_concat_coeffs, threshold);

groups = find(bin_sizes >= n_group_min); % Find bin number for groups that contain minimum # of discharges 
n_groups = length(groups);
group_sizes = bin_sizes(groups); % # of discharges in each group

% Find discharges for each group
grouped_discharges = cell(n_groups, 1); % ID# for discharges in a group
for i = 1:n_groups
    grouped_discharges{i} = find(bin_assignment == groups(i));
end

n_discharges = sum(group_sizes); % number of peaks that was sorted into a group


% ***** Confirm classification procedure worked as intended *****
% -Demonstrate that discharges are more similar to other discharges in own
% group vs. other groups 


% Allocate space to store inter-group comparisons (i.e. corr coeffs)
inter_group_coeffs_avg = zeros(n_groups, n_groups);
inter_group_coeffs_indiv = cell(n_groups, n_groups);

for i = 1:n_groups
    
    for j = 1:n_groups
        
        if i == j % If comparing discharges within a group
            
            % Take top half of coeff matrix containing corr coeffs of
            % discharges within group (diagonal not included)
            comparison_coeffs = triu(...
                rhd_epochs_concat_coeffs(grouped_discharges{i}, grouped_discharges{j}), 1);
            
            comparison_coeffs(comparison_coeffs == 0) = NaN; % Replace zeros with NaN
            
            inter_group_coeffs_avg(i,j) = ...
                mean(comparison_coeffs(:), 'omitnan'); % Store average corr coeff, excluded NaN
            
            % Collapse matrix, remove NaN
            comparison_coeffs = comparison_coeffs(:);
            comparison_coeffs(isnan(comparison_coeffs)) = [];
            inter_group_coeffs_indiv{i,j} = comparison_coeffs;
            
        else
            
            % If comparing across groups, get corr coeffs of discharges in
            % group i vs. group j
            comparison_coeffs = rhd_epochs_concat_coeffs(grouped_discharges{i}, grouped_discharges{j});
            inter_group_coeffs_avg(i,j) = mean(comparison_coeffs(:));
            inter_group_coeffs_indiv{i,j} = comparison_coeffs(:);
            
        end   
    end
end

% Calculate noise estimates & SNR for RHD signals

group_rhd_avg = cell(1, n_groups);
group_rhd_noise = cell(1, n_groups);
group_rhd_snr = cell(1, n_groups);

for i = 1:n_groups
    
    group_epochs = rhd_epochs(:,:,grouped_discharges{i});
    avg_signal = mean(group_epochs, 3);
%     half = ceil(size(group_epochs,3)/2);
%     
%     plus_min_avg = mean(group_epochs(:,:,1:half),3) - ...
%         mean(group_epochs(:,:,half+1:end),3);
    plus_min_avg = ....
        mean(group_epochs(:, :, 1:2:end), 3) - mean(group_epochs(:, :, 2:2:end), 3);
    rms_noise = rms(plus_min_avg');
    rms_signal = rms(avg_signal');
    
    group_rhd_avg{i} = avg_signal;
    group_rhd_noise{i} = plus_min_avg;
    group_rhd_snr{i} = 20*log10(rms_signal./rms_noise);

end


all_figures(2) = figure;
ic_offset = 900;
epoch_trigger = epoch_length/2 + 1;

for i = 1:n_groups
   
   group_epochs = rhd_epochs(:, :, grouped_discharges{i});
   
   subplot(2,n_groups, i); hold('on')
   subplot_title = {['Group ', num2str(i)]; ['n = ', num2str(group_sizes(i))]};
   title(cellstr(subplot_title))
   
   % Plot individual traces
   for j = 1:group_sizes(i)
       
       plot_offset_signals(...
           group_epochs(:,:,j), ic_offset, ...
           'Color', my_color("blue_gray"))
   end
   
   % Plot average signal
   plot_offset_signals(group_rhd_avg{i}, ic_offset, 'k', ...
       'LineWidth', 1.5)
   
   ylim([(n_rhd_channels + 1)*-1*ic_offset, 0])
   xlim([0, epoch_length])
   
   % Plot line for trigger
   plot([epoch_trigger, epoch_trigger], ylim, 'r')
   
   % Plot noise (plus/min average)
   plot_offset_signals(group_rhd_noise{i}, ic_offset, ...
       'Color', my_color("forest"), 'LineWidth', 1)
   
   % Boxplots of inter-group correlations 
   subplot(2,n_groups, i+n_groups); hold('on')
   
   boxplot_data = inter_group_coeffs_indiv(i,:);
   group_labeled_boxplot_data = [];
   
   % Make an extra column with group numbers to use as grouping variable in
   % boxplot function
   for j = 1:n_groups    
       x = [(boxplot_data{j}), (repelem(j, length(boxplot_data{j})))'];
       group_labeled_boxplot_data = [group_labeled_boxplot_data; x];
   end
   
   % Set color for group being plotted to 1 so that it's a different color
   % from other groups it is being compared to 
   group_color = zeros(n_groups, 1);
   group_color(i) = 1;
   
   boxplot(group_labeled_boxplot_data(:,1), group_labeled_boxplot_data(:,2),...
       'Color', 'kb', 'ColorGroup', group_color, 'Symbol', '')
   xlim([0,n_groups+1])
   ylim([-1, 1])
end


% *** Analysing scalp signals *****

% Trim scalp clips into epochs 
% Don't include channels FP1, FP2, M1, and M2 
scalp_epochs = scalp_avg_ref_mont(included_channels, epoch_range, :);

% For each group, make average of each scalp channel 
% Calculate plus/minus average & SNR
group_scalp_avg = cell(1,n_groups);
group_scalp_noise = cell(1, n_groups);
group_scalp_snr = cell(1, n_groups);

for i = 1:n_groups
    
    group_epochs = scalp_epochs(:,:,grouped_discharges{i});
    avg_signal = mean(group_epochs, 3);
    
    plus_min_avg = ...
        mean(group_epochs(:, :, 1:2:end), 3) - mean(group_epochs(:, :, 2:2:end), 3);
    rms_noise = rms(plus_min_avg');
    rms_signal = rms(avg_signal');
    
    group_scalp_avg{i} = avg_signal;
    group_scalp_noise{i} = plus_min_avg;
    group_scalp_snr{i} = 20*log10(rms_signal./rms_noise);
    
end


% Plot average scalp signals for each group 
scalp_offset = 25; % Space between traces 
scalp_ylim = [(size(scalp_epochs,1)+1)*-1*scalp_offset, 0];
channel_label_x = -75; % location on x axis to write channel label
snr_label_x = epoch_length + 25;


all_figures(3) = figure;

for i = 1:n_groups
    
    group_epochs = scalp_epochs(:,:, grouped_discharges{i});   
    snr = group_scalp_snr{i};
    
    subplot(1, n_groups, i); hold('on')
    subplot_title = {['Group ', num2str(i)]; ['n = ', num2str(group_sizes(i))]};
    title(cellstr(subplot_title))
    
    for j = 1:group_sizes(i)
        plot_offset_signals(group_epochs(:,:,j), scalp_offset,...
            'Color', my_color("blue_gray"))
    end
    
    for j = 1:size(scalp_epochs,1)
        text(channel_label_x, -j*scalp_offset,...
            scalp_channel_labels{included_channels(j)},...
            'Color', 'k')
        
        % If SNR is significant, write text in red
        if snr(j) >= sig_snr_dB
            text(snr_label_x, -j*scalp_offset,...
                num2str(round(snr(j),2)),...
                'Color', 'r')
        else
            text(snr_label_x, -j*scalp_offset,...
                num2str(round(snr(j),2)),...
                'Color', 'k')
        end
    end
    plot_offset_signals(group_scalp_avg{i}, scalp_offset, 'k', 'LineWidth', 1.5)
    plot_offset_signals(group_scalp_noise{i}, scalp_offset, 'Color', my_color("forest"),...
        'LineWidth', 1)
    
    xlim([0, epoch_length])
    ylim(scalp_ylim)
    ax = gca;
    ax.YAxisLocation = 'left';
    plot([epoch_trigger, epoch_trigger], ylim, 'r')

end


%%

% Cross-corr with inidividual RHD channels 
group_rhd_avg_by_channel = cell(1, n_groups);
rhd_x_scalp_correlograms = cell(1, n_groups);
xcorr_analysis_info = cell(1, n_groups);
xcorr_analysis_results = cell(1, n_groups);

for i = 1:n_groups 
    
    % Average signal of 1 RHD channel in group
    group_rhd_avg_by_channel{i} = mean(rhd_epochs(:,:, grouped_discharges{i}),3);
    
    % Allocate space
    group_correlograms = cell(1, size(scalp_epochs,1));
    group_max_abs_coeffs = cell(1, size(scalp_epochs,1));
    max_channels = zeros(size(scalp_epochs,1),1);
    max_channel_coeffs = zeros(size(max_channels));
    max_channel_lags = zeros(size(max_channels));
   
    for j = 1:size(scalp_epochs, 1)
        
        % Correlogram against each RHD channel
        single_channel_correlograms = zeros(n_rhd_channels,epoch_length*2-1);
        single_channel_max_coeffs = zeros(n_rhd_channels,1); % Max coeff for each RHD channel correlogram
        single_channel_max_coeff_locs = zeros(n_rhd_channels,1); % index for max coeff
        single_channel_max_coeff_lags = zeros(n_rhd_channels,1); % corresponding lag for max coeff
        
        for k = 1:n_rhd_channels
            
            % Cross-correlate STA from 1 RHD channel with STA from 1 scalp
            % EEG channel
            [coeffs, xcorr_lags] = xcorr(...
                group_rhd_avg_by_channel{i}(k,:), group_scalp_avg{i}(j,:), 'coeff');
            
            single_channel_correlograms(k,:) = coeffs; % Store raw correlogram
            i_max_coeff = find(abs(coeffs) == max(abs(coeffs))); % Find index for greatest absolute corr coeff
            
            single_channel_max_coeffs(k) = coeffs(i_max_coeff); % Store max coeff 
            single_channel_max_coeff_locs(k) = i_max_coeff; % Store index for max coeff
            single_channel_max_coeff_lags(k) = xcorr_lags(i_max_coeff); % Store lag for max coeff
        end
        
        group_correlograms{j} = single_channel_correlograms; 
        
        group_max_abs_coeffs{j} = table(single_channel_max_coeffs, single_channel_max_coeff_locs, single_channel_max_coeff_lags);
        
        % Find RHD channel with greatest correlation with scalp
        max_channel = find(abs(single_channel_max_coeffs) == max(abs(single_channel_max_coeffs)));
        max_channels(j) = max_channel;
        max_channel_coeffs(j) = single_channel_max_coeffs(max_channel);
        max_channel_lags(j) = single_channel_max_coeff_lags(max_channel);
        max_coeff_results = table(included_channel_labels, max_channels, max_channel_coeffs, max_channel_lags);
        
    end
    
    rhd_x_scalp_correlograms{i} = group_correlograms;
    xcorr_analysis_info{i} = group_max_abs_coeffs; % Tables containing all xcorr info for all RHD and scalp channels
    xcorr_analysis_results{i} = max_coeff_results; % Only data relevant to maximum xcorr 
    
end




% Plot cross-correlation data
xcorr_offset = 2.5;
text_offset = 500;
xcorr_threshold = 0.7;

all_figures(4) = figure;
for i = 1:n_groups
    
    snr = group_scalp_snr{i};
    subplot(1, n_groups, i); hold('on')
    for j = 1:size(max_coeff_results, 1)
        
        plot_results = xcorr_analysis_results{i};
        
        % Plot correlogram with the greatest RHD x scalp EEG coeff
        plot(xcorr_lags,rhd_x_scalp_correlograms{i}{j}(plot_results.max_channels(j),:)-j*xcorr_offset, 'b')
        plot([0,0], ylim, 'r')
        
        % Print max correlation coefficient & corresponding lags
        % Colors depend on SNR
        
        % If scalp EEG had sufficient SNR, print correlation coefficient in
        % green, lag in red. Otherwise, print numbers in black 
        if snr(j) >= sig_snr_dB
            text(0-text_offset/2, -(j-1)*xcorr_offset-1, ...
                num2str(round(plot_results.max_channel_coeffs(j), 2)), 'Color', 'g')
            
            text(0+text_offset-100, -(j-1)*xcorr_offset-1,...
                num2str(plot_results.max_channel_lags(j)), 'Color', 'r')
        else
            text(0-text_offset/2, -(j-1)*xcorr_offset-1, ...
                num2str(round(plot_results.max_channel_coeffs(j), 2)), 'Color', 'k')
            
            text(0+text_offset-100, -(j-1)*xcorr_offset-1,...
                num2str(plot_results.max_channel_lags(j)), 'Color', 'k')
        end


        % Print RHD channel number
        text(0+text_offset/2, -(j-1)*xcorr_offset-1, num2str(plot_results.max_channels(j)), 'Color', 'b')
        
        % Print scalp channel labels
        text(-600, -(j-1)*xcorr_offset-1,...
            scalp_channel_labels{included_channels(j)},...
            'Color', 'k')
    end
    ax = gca;
    axis off
end




