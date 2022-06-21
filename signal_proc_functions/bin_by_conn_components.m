  function [bin_assignment, n_bins, bin_sizes, conncomp_graph] = bin_by_conn_components(corr_matrix, threshold)
    
        corr_matrix = corr_matrix - diag(diag(corr_matrix));
        corr_matrix(find(corr_matrix >= threshold)) = 1;
        corr_matrix(find(corr_matrix < threshold)) = 0;
        
        conncomp_graph = graph(corr_matrix, 'upper');
        [bin_assignment, bin_sizes] = conncomp(conncomp_graph);
        n_bins = length(bin_sizes);

  end