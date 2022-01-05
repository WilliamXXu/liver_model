load('q5.mat');
y=res1;
step=3
x=(1:step:16);
plot(x,y,'-o')
xlabel('Wall thickness, $\frac{\epsilon}{4}$','Interpreter','latex');
% xlabel('$$','interpreter','latex');
 ylabel('Q values');