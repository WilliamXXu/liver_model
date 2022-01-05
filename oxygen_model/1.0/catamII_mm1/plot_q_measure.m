k=1;
a=0.8;
delta_li=(8:8:96);
res=[];
N=128;
for delta=delta_li
    ep=round(k*(delta.^a));
    delta
    ep
    matr=poisson_insulate_inf(N,ep,delta);
    res=[res,q_measure(matr)];
end
k=0.5;
a=1;
res1=[];
for delta=delta_li
    ep=round(k*(delta.^a));
    delta
    ep
    matr=poisson_insulate_inf(N,ep,delta);
    res1=[res1,q_measure(matr)];
end

plot(delta_li,res,'Displayname','k=1,a=0.8');
hold on
plot(delta_li,res1,'Displayname','k=0.5,a=1');
xlabel('$\delta$','Interpreter','latex');
ylabel('Q value');
legend('Location','northwest');
hold off
