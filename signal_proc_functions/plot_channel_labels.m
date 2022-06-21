% Adds channel labels to a plot 
function plot_channel_labels(channel_labels, y_offset)
    
    n_channels = numel(channel_labels);
    
    
    if y_offset == 'default'
        
        y_offset = diff(ylim)/n_channels;
        
    end
    
    yticks(-1*n_channels*y_offset:y_offset:-y_offset)
    yticklabels(flipud(channel_labels))
    
    
end
