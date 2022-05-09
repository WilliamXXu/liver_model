



% plot([0 4.5],[0 1])
% hold on
% ax = gca;
% ax.FontSize = 17;
% title('Simulated Cellular ATP After Drug Intake (DMO)')
% xlabel('Minutes') 
% ylabel('% of Basal Cellular ATP') 

x=linspace(0,800,100);
y=arrayfun(@dr,x);
plot(x,y*100);
hold on
xline(120,'--r')
xline(390,'--r')
yline(50,'--r')
hold off
title('Drug Concentration in Option 3 of Drug Submodel')
yticks([0 25 50 75 100])
ylabel('% of C_{max}')
set(gca,'FontSize',18)
set(gca,'xtick',[])
legend('Drug concentration')
function [conc]=dr(t)   %ez drug concentration
    t=t/60;
    t_max=2;
    c_max=1;
    grad=c_max/t_max;
    half_life=4.5;
    k=log(2)/half_life;

    conc=(t<=t_max)*t*grad+(t>t_max)*c_max*exp(-k*(t-t_max));

end