ind=5;
load('results.mat');
x=linspace(0,60,30);
y=deval(sol,x,ind);
plot(x,y);
