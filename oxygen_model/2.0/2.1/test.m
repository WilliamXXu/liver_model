
k1=9200;
k2=26^2.73;
k3=5700;

tic
syms x
digits(32)
res=vpasolve(x^3.73+(k1-k3)*x^2.73+k2*x-k2*k3==0,30)
toc