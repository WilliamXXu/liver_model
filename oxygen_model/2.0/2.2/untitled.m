tspan = [0 3000];
y0 = [2 0];
tic
sol = ode15s(@vdp1000,tspan,y0);
toc
x = linspace(0,3000,2500);
tic
y = deval(sol,5.2452,1)
toc