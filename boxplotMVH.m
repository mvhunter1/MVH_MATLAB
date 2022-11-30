
% ************************************************************************
% (c) 2022:
%       Miranda Hunter, White Lab, MSKCC
%           hunterm@mskcc.org | mirandavhunter@gmail.com
%
% Boxplot function.
% data_matrix = data you want to plot as matrix with each condition = column.
% colours = optional, colours for each group.
% xlabels = text labels for each group.
% dot_size = self explanatory, default = 75.
% thickness = linewidth of box, default = 3.
% ************************************************************************



function boxplotMVH(data_matrix, colours, xlabels, dot_size, thickness)

[~,n_conditions] = size(data_matrix);

if nargin == 1
    dot_size = 75;
    thickness = 3;
    colours = parula(3);
    xlabels = [];
elseif nargin == 2
    dot_size = 75;
    thickness = 3;
    xlabels = [];
elseif nargin == 3
    dot_size = 75;
    thickness = 3;
end

for ii = 1:n_conditions
    if ii == 1
        xrange = [1 2];
    else
        xrange = [1+(1.5*(ii-1)) 2+(1.5*(ii-1))];
    end

    data = rmmissing(data_matrix(:,ii));
    colour = colours(ii,:);

    % plot individual points
    x = linspace(xrange(1), xrange(2), length(data)); % generating evenly spaced points along the x axis for dots to be plotted
    scatter(x,data,dot_size,colour,'filled');
    hold on

    % plot box
    mean_data = mean(data);
    max_data = mean_data + std(data);
    min_data = mean_data - std(data);
    firstQ = mean_data - std(data)/sqrt(length(data)); % standard error
    thirdQ = mean_data + std(data)/sqrt(length(data)); % standard error

    % plot vertical lines
    plot([mean(xrange), mean(xrange)], [thirdQ, max_data], 'k', 'LineWidth', thickness);
    plot([mean(xrange), mean(xrange)], [firstQ, min_data], 'k', 'LineWidth', thickness);
    plot([xrange(1), xrange(1)], [firstQ, thirdQ], 'k', 'LineWidth', thickness);
    plot([xrange(2), xrange(2)], [firstQ, thirdQ], 'k', 'LineWidth', thickness);

    % plot horizontal lines
    plot([xrange(1), xrange(2)], [thirdQ, thirdQ], 'k', 'LineWidth', thickness);
    plot([xrange(1), xrange(2)], [firstQ, firstQ], 'k', 'LineWidth', thickness);
    plot([xrange(1), xrange(2)], [mean_data, mean_data], 'k', 'LineWidth', thickness);

end
set(gca, 'Box', 'off', 'FontSize', 26, 'LineWidth', 3);
set(gca, 'XTick', (1.5:1.5:100), 'XTickLabel', xlabels);
xlim([0.5 n_conditions+2]);
ylim([0 Inf]);

end