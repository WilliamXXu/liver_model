% res=load('result.mat').res;
% res=res(1:3,1:137);
% h=heatmap(res);
% colormap copper
% set(gca,'FontSize',18)
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% h.Title = 'Plasma Oxygen Concentration (mmHg)';
% h.YLabel = 'Radial';
% h.XLabel = 'Axial';
% grid off


% res=load('ocr.mat').ocr;
% res=res(5:17,1:137);
% h=heatmap(res);
% colormap summer
% set(gca,'FontSize',18)
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% h.Title = 'Hepatocyte Oxygen Consumption Rate (mmHg/s)';
% h.YLabel = 'Radial';
% h.XLabel = 'Axial';
% grid off

% 
% res=load('result.mat').res;
% res=res(5:17,1:137);
% h=heatmap(res);
% colormap jet
% set(gca,'FontSize',18)
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% h.Title = 'Hepatocyte Oxygen Concentration (mmHg)';
% h.YLabel = 'Radial';
% h.XLabel = 'Axial';
% grid off




% 
% res=load('hemo.mat').hemo;
% res=res(1:3,1:137);
% h=heatmap(res);
% colormap gray
% set(gca,'FontSize',18)
% Ax = gca;
% Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% h.Title = 'Haemoglobin Oxygen Concentration (mmHg)';
% h.YLabel = 'Radial';
% h.XLabel = 'Axial';
% grid off



% res=load('result.mat').res;
% res=res(5:17,1:137);
% 
% res=flip(res,1);
% contourf(res,10);
% colorbar
% ax = gca;
% ax.FontSize = 17;
% set(gca,'XTick',[], 'YTick', [])
% title('Hepatocyte Oxygen Concentration (mmHg)')
% xlabel('Axial');
% ylabel('Radial');

% res=load('ocr.mat').ocr;
% res=res(5:17,1:137);
% 
% colormap cool
% 
% res=flip(res,1);
% contourf(res,10);
% colorbar
% ax = gca;
% ax.FontSize = 17;
% set(gca,'XTick',[], 'YTick', [])
% title('Hepatocyte Oxygen Consumption Rate (mmHg/s)')
% xlabel('Axial');
% ylabel('Radial');
