k1=9200/1.34;
k2=26^2.73;
k3=5700;
f=@(x) x^3.73+(k1-k3)*x^2.73+k2*x-k2*k3;
df=@(x) 3.73*x^2.73+2.73*(k1-k3)*x^1.73+k2;
df2=@(x) 2.73*3.73*x^1.73+1.73*2.73*(k1-k3)*x^0.73;
%plot(arrayfun(df,[20:5:65]))
tic
a=halley(f,df,df2,42,1)
toc

function [sol] = halley(f,df,df2,x0,epsilon)
ling=f(x0);
yi=df(x0);
er=df2(x0);
sol = x0 - 2*ling*yi/(2*yi^2-ling*er);


count=0;
while abs(f(sol)) > epsilon
    sol = sol - f(sol)/df(sol);
    count=count+1
end
end

