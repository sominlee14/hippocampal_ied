% Plots signals in single plot using an offset 
% Signals is an array/matrix containing a separate signals or channel in
% each row
% Offset is the amount of space between each trace. 'default' calculates
% offset by using the biggest range found in signal
% varargin accepts parameters for plot function 

function plot_offset_signals(signals, offset, varargin)

    
    if offset == 'default'
        
        plot_offset = max(range(signals,2)) * 0.75;
        
    else
        
        plot_offset = offset;
        
    end
    
    n_channels = size(signals, 1);
    
    for i = 1:n_channels
        
        plot((signals(i,:) - (i*plot_offset)), varargin{:})
        ax = gca;
        ax.ColorOrderIndex = 1;
        
    end
    
    xlim([0,size(signals,2)])
    

end
