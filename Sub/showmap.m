function hf = showmap(Map, meterPerPixel, fig)
% showmap(Map, meterPerPixel, fig)
% 
% Caution: for the matrix Map(x, y), the "origin" is at the top-left corner,
% while for the image to display, the "origin" is at the bottom-left corner.
% In other words, the x-axis in the image corresponds to the rows of the
% matrix Map(x, y), and the y-axis in the image corresponds to the columns.

if nargin == 3
    hf = figure(max(1, round(fig)));
else
    hf = [];
end
    
[lenX, lenY] = size(Map);

[Xmeter, Ymeter] = meshgrid((0:lenX - 1) * meterPerPixel, ...
                            (0:lenY - 1) * meterPerPixel);
surf(Xmeter, Ymeter, Map.', 'edgecolor', 'none');
view(0, 90);
axis square
xlim([0, (lenX - 1) * meterPerPixel]);
ylim([0, (lenY - 1) * meterPerPixel]);
% set(gca, 'FontSize', 14);
xlabel('Longitude (meter)');
ylabel('Latitude (meter)');

% XLableDistance_y = - 0.09;     % Units = 'normalized'
% YLableDistance_x = - 0.10;      % Units = 'normalized'
% 
% xlab = get(gca, 'XLabel');
% set(xlab, 'Units', 'normalized');
% set(xlab, 'Position', [0.5 XLableDistance_y 0]);
% 
% ylab = get(gca, 'YLabel');
% set(ylab, 'Units', 'normalized');
% set(ylab, 'Position', [YLableDistance_x 0.5 0]);