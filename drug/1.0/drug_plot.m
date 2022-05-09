% res=dru;
% h=heatmap(res);
% colormap jet
% set(gca,'FontSize',18)
% Ax = gca;
%  Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% %ax.XDisplayLabels=linspace(0,275,10);
% 
% h.Title = 'High Permeability Drug Concentration (mol/m^3)';
% h.YLabel = 'Radial';
% h.XLabel = 'Axial (\mu m)';
% grid off
% 



res=dru;
res=flip(res,1);
contourf(res,10);
colorbar
ax = gca;
ax.FontSize = 18;
set(gca,'XTick',[], 'YTick', [])
title('Low Permeability Drug Concentration (mol/m^3)')
xlabel('Axial');
ylabel('Radial');