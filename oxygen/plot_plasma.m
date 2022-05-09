clear
clc
res=load('result.mat').res;
x=linspace(0,275,137);
plot(x,res(1,1:137))
ax = gca;
ax.FontSize = 18;
title('Simulated Plasma Oxygen Concentration In Sinusoid')
xline(275/3,'--r')
xline(275/3*2,'--r')
axis tight
xlabel('Distance From Hepatic Triad (\mum)');
ylabel('Plasma Oxygen Concentration (mmHg)');