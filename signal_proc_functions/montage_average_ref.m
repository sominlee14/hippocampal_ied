function montaged_signal = montage_average_ref(signals, reference_signals)

    montaged_signal = signals - mean(reference_signals, 1);

end
