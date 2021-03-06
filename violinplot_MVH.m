% 
% ************************************************************************
% (c) 2022:
%       Miranda Hunter, White Lab, MSKCC
%           hunterm@mskcc.org | mirandavhunter@gmail.com
% 
% Create violin plot with data points plotted on top of violins.
% ************************************************************************

% Requires (from MATLAB File Exchange): 
%   violin.m
%   plotSpread.m

function violinplot_MVH(data, x_labels, cols)

[~,c] = size(data);

% set colours for 2 and 3 component violins
if nargin < 3 && c == 2
    cols = [227 91 114; 91 227 204] ./256;
    markersize = 12;
elseif nargin < 3 && c == 3
    cols = [227 91 114; 91 227 204; 255 172 53] ./256;
    markersize = 10;
elseif nargin < 3
    cols = parula(c);
    markersize = 10;
else
    markersize = 10;
end

% plot
violin(data, 'edgecolor', 'none', 'facecolor', cols);
plotSpread(data, 'distributionColors', 'k', 'xNames', x_labels, 'spreadWidth', 0.5)
box off
if c == 2
    xlim([0.5 2.5]);
end
ax = gca;
ax.FontSize = 24;
ax.LineWidth = 2;
h = get(gca, 'children');
for ii = 1:length(h)
    if strcmp(get(h(ii), 'type'), 'line')
        set(h(ii), 'markersize', markersize);
    end
end
end

