function avg_spikes(time, spikes, idx, num_clusters, rows, cols, channel)
    figure()
    for i=1:num_clusters
        avg_neurons(i,:) = mean(spikes(idx==i, :), 1);
        std_neurons(i,:) = std(spikes(idx==i, :), 1);
    
        subplot(rows,cols,i)
        plotEB(time, avg_neurons(i,:), std_neurons(i,:))
        if (i>cols)
            xlabel('Time (ms)')
        end
    
        if (i==1)||(i==cols+1)
            ylabel('Voltage (uV)')
        end
        title(['Cluster ', num2str(i)])
    end
    sgtitle(['Average Neuons Based on Clustering, Channel ', num2str(channel)])
end