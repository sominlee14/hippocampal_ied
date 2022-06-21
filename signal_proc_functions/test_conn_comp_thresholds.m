  
  function [bin_info] = test_conn_comp_thresholds(corr_matrix, thresholds, n_minimum)
    
    n_bins = zeros(length(thresholds), 1);
    bin_sizes = cell(length(thresholds), 1);
    n_groups = zeros(length(thresholds), 1);
    
    for i = 1:length(thresholds)
        
        corr_matrix_binary = corr_matrix - diag(diag(corr_matrix));
        corr_matrix_binary(find(corr_matrix_binary >= thresholds(i))) = 1;
        corr_matrix_binary(find(corr_matrix_binary < thresholds(i))) = 0;
        
        corr_graph = graph(corr_matrix_binary, 'upper');
        [~, bin_size] = conncomp(corr_graph);
        
        n_bins(i) = length(bin_size);
        n_groups(i) = length(find(bin_size >= n_minimum));
        bin_sizes{i} = sort(bin_size, 'descend');
        
    end
    
    thresholds = thresholds';
    bin_info = table(thresholds, n_bins, bin_sizes, n_groups);
    
  end