delta_list=(1:0.1:2);
iter=[];
for delta=delta_list
    delta
    [x,counter]=poisson(128,delta);
    iter=[iter,counter];
end
plot(delta_list,iter);
xlabel('$\delta$','Interpreter','latex');
ylabel('The number of iterations');
