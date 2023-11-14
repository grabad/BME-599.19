function plotEB(x, y, error)
    % x: x-axis values
    % y: y-axis values
    % error: error values for each data point
    % lineColor: color of the line
    % faceColor: color of the shaded error bars
    
    % Plot the line
    hold on;

    % Plot shaded error bars
    xError = [x, fliplr(x)];
    yError = [y + error, fliplr(y - error)];
    fill(xError, yError, [0 0.4470 0.7410], 'FaceAlpha', 0.3, 'EdgeAlpha', 0);
    
    plot(x, y, 'LineWidth', 2);

    hold off;
end