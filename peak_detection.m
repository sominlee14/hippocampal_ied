% Try analysis with absolute value peaks


% Run peak detection on patient XLTEK files and save clipped signals into
% matfile

clear; clc;

pt = 'PX' % file ID for XLTEK .raw files 


% File names for saving

diary_name = [pt, '_log_', datestr(now, 'yyyy_mm_dd')]
diary(diary_name) % initiate diary to keep track command term

% *** Peak Detection Parameters ****
bp_freq = [10,50] %#ok<*NOPTS>
filter_order = 2
min_peak_dist = .5 % in seconds v(flank around peaks)
peak_zscore = 5
clip_length = 2 % Length of clip in seconds (centered around peak)

% Run peak detection 
[signal_clips, channel_list, n_rhd_channels, rhd_channels, sample_rate] = ...
    detect_abs_peaks_rhd(...
    bp_freq, filter_order, min_peak_dist, peak_zscore, clip_length);

%Save all relevant variables for future analysis 
    

save([pt, '_metadata_', datestr(now, 'yyyy_mm_dd'), '.mat'], 'channel_list', 'sample_rate')


save([pt, '_signal_clips_', datestr(now, 'yyyy_mm_dd'), '.mat'], 'signal_clips', '-v7.3')

diary off


% Internal peak detection function 
function [signal_clips, ...
    channel_list, ...
    n_rhd_channels, ...
    rhd_channels,...
    sample_rate] = detect_abs_peaks_rhd(bp_freq, filter_order, min_peak_dist, peak_zscore, ...
    clip_length)

% Given a set of xltek files, runs peak detection and saves results &
% signal clips 

% (1) Read signals from XLTEK files & run peak detection to find epochs

% Read in channel information file (.xls)
channel_list = read_channel_file();

% Find channel numbers for RHD electrodes 
rhd_channels = find(ismember(channel_list.channel_type, 'RHD'));
n_rhd_channels = length(rhd_channels);

% Choose XLTEK files to use for analysis (user input) 
[xltek_filename, xltek_filepath] = uigetfile('*.raw', ...
    'Select .raw file(s)',...
    'MultiSelect', 'on');
% Exit if dialogue box is cancelled
if isequal(xltek_filename, 0)
    fprintf('File selection cancelled')
    return
end

n_files = length(xltek_filename);
fprintf('%d files selected. \n \n', n_files)
full_path = fullfile(xltek_filepath, xltek_filename);

signal_clips = [];
total_length_seconds = 0;

% (2) Peak detection & signal clipping

for i = 1:n_files
    
    % Read in single xltek file
    filename = full_path{i};
    
    % Extract signals & sample rate
    [raw_signals, ~, sample_rate] = readXLTEK(filename);
    raw_signals = subtract_mean(raw_signals); % De-mean raw recording
    
    % Keep track of total length of all files 
    file_length_seconds = round(size(raw_signals, 2)/sample_rate);
    total_length_seconds = total_length_seconds + file_length_seconds;
    
    fprintf('Information for %s \n \n', filename)
    fprintf('Length of recording (seconds): %d \n', file_length_seconds)
    fprintf('Sample rate: %d \n', sample_rate)
    
    
    % Make RHD signal to use for peak detection
    rhd_signals = raw_signals(rhd_channels,:); % Get only RHD channels
    
    % Filter and make average RHD signal 
    peak_detection_signal = mean(...
        filter_butter(rhd_signals, bp_freq, sample_rate, filter_order), 1);
    
    % Find peak using cutoff
    % Peak using absolute value (so both up and downward peaks)
    min_peak_height = std(peak_detection_signal) * peak_zscore
%     min_peak_height = median(abs(peak_detection_signal))/0.6745 * peak_zscore
    
    [~, peak_locs] = findpeaks(...
        abs(peak_detection_signal), sample_rate,...
        'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_dist);
    i_peaks = peak_locs * sample_rate + 1; % findpeaks returns time of peak. Convert back to # samples
    
    fprintf('# of peaks found: %d \n \n \n', length(i_peaks))
    
    % Don't include peaks that were found too close to beginning or end of
    % recording and can't be clipped in full 
    n_samples = sample_rate * clip_length;
    i_remove = find(i_peaks <= n_samples/2 | ...
        i_peaks > size(rhd_signals, 2)-n_samples/2);
    i_peaks(i_remove) = [];
    
    fprintf('# peaks removed: %d \n\n\n', length(i_remove))
    
    % Save all clips found for each file
    signal_clips = cat(3, signal_clips, ...
        get_clips(raw_signals, i_peaks - (n_samples/2), n_samples));
    
end

fprintf('Peak detection parameters: \n')
fprintf('zscore = %d \n', peak_zscore)
fprintf('Min Peak Distance = %d seconds \n', min_peak_dist)
fprintf('Total # peaks: %d \n\n', size(signal_clips, 3))
fprintf('Total time all files: %d seconds \n\n', total_length_seconds)
  
end

    
