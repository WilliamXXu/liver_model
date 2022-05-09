load('/Users/williamxu/MATLAB-Drive/final/test1/data/mito6_365.mat')
res=res(5:17,1:137)*100;
res=flip(res,1);
hold on
contourf(res,10);
[C,h]=contour(res,[25 25],'--r');
colorbar
clabel(C,h)
ax = gca;
ax.FontSize = 17;
set(gca,'XTick',[], 'YTick', [])
title('Simulated Cell Density (% of Basal Value)')
xlabel('Axial');
ylabel('Radial');
% h.Title = 'Simulated Cell Density (% of Basal Value)';
% h.YLabel = 'Radial';
% h.XLabel = 'Axial';

% h=heatmap(res);
% colormap parula
% set(gca,'FontSize',18)
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% h.Title = 'Simulated Cell Density (% of Basal Value)';
% h.YLabel = 'Radial';
% h.XLabel = 'Axial';
% grid off


% res=oxy;
% res=res(5:17,1:137);
% h=heatmap(res);
% colormap winter
% set(gca,'FontSize',18)
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% h.Title = 'Simulated Hepatocyte Oxygen Concentration (mmHg)';
% h.YLabel = 'Radial';
% h.XLabel = 'Axial';
% grid off



