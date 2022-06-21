% Runs simulation to figure out jitter for "0" lag 

sample_rate = 1024; 

% Frequency for sine signal
freq = 10;
% number of trials for monte carlo siulation 
n_trials = 10000;

% noise levels to test 
noise_levels = 0.1:0.1:.9;

% Allocate for saving standard deviations and SNRs for results 
stdevs = nan(length(noise_levels), 1);
snrs = nan(length(noise_levels), 1);

% Loop to run for each noise level 
for i = 1:length(noise_levels)
    
    % allocate to save lag for each trial
    max_lags = nan(n_trials, 1);

    % run each trial 
    for j = 1:n_trials
        
        % Generate two sine waves 
        [sig_1, ~] = generate_toy_signal(freq, sample_rate, 0.5, 1, noise_levels(i), 's');
        [sig_2, ~] = generate_toy_signal(freq, sample_rate, 0.5, 1, noise_levels(i), 's');
        
        % Cross-correlate and find lag of maximum correlation 
        [xcorrelogram, lags] = xcorr(sig_1, sig_2, 'coeff');
        [max_xcorr, i_max] = max(xcorrelogram);
        max_lags(j) = lags(i_max);
        
    end
    
    pd = fitdist(max_lags, 'Normal');
    x = lags;
    y = pdf(pd, x);

    stdevs(i) = pd.std;
    
    s = nan(n_trials, 1);
    
    % calculate SNR 
    for k = 1:n_trials
        
        [sig_1, time_axis] = generate_toy_signal(freq, sample_rate, 0.5, 1, 0, 's');
        noise = randn(1, length(sig_1))* noise_levels(i);

        sig_w_noise = sig_1 + noise;

        rms_sig = rms(sig_1);
        rms_noise = rms(noise);
        s(k) = rms_sig/rms_noise;
    end
    
    snrs(i) = 20*log10(mean(s));
end

%%
figure
plot(snrs, 3*stdevs, 'ko-')
xlabel('SNR')
ylabel('jitter (+/- ms)')
ylim([0, 12])