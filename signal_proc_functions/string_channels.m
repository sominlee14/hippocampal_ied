function [signal_vector] = string_channels(signals)
  % String channels of a multichannel signal to make 1D signal
  
  n_channels = size(signals, 1);
  signal_length = size(signals, 2); 
  n_signals = size(signals,3);
  
  signal_vector = zeros(n_signals, n_channels*signal_length);
  
  for i = 1:n_signals
      
      signal_vector(i,:) = reshape(signals(:,:,i)', [1, n_channels*signal_length]);
      
  end
end