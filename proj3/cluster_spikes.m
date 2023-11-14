function spikes_coeff = cluster_spikes(spikes, k_max, channel)
    
    spikes_coeff = pca(zscore(spikes)');

    figure()
    scatter(spikes_coeff(:, 1), spikes_coeff(:, 2), 'filled')
    xlabel('PCA 1')
    ylabel('PCA 2')
    zlabel('PCA 3')
    title(['Unclustered Neurons, Channel ', num2str(channel)])
    
    inertias = zeros([1 k_max]);
    
    for i=1:k_max
        [~, ~, sumD] = kmeans(spikes_coeff(:, 1:3), i, "MaxIter", 500);
        inertias(i) = sum(sumD);
    end
    
    figure()
    plot(1:k_max, inertias)
    xlabel('Number of Clusters')
    ylabel('Inertia')
    title(['Inertia vs. Number of Clusters, Channel ', num2str(channel)])
end