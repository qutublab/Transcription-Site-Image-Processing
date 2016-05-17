function h = figureInitPlot(x, y, marker, cmap, legendText)

h = figure;
plot(x, y, marker{1}, 'MarkerFaceColor', cmap(1, :), 'MarkerEdgeColor', cmap(1, :));
hold on;
for i = 2:size(cmap, 1)
    plot(x, y, marker{i}, 'MarkerFaceColor', cmap(i, :), 'MarkerEdgeColor', cmap(i, :));
end
legend(legendText, 'Location', 'northeastoutside');
end

