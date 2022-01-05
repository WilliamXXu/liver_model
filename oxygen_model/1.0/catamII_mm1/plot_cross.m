
    s=50;
    e=80;
    k=65;
    x=poisson_insulate(128,4,8);
    part=x(s:e,k);
    part2=x(s:e,k+2);
    part3=x(s:e,k+4);
    
    plot((s-1:e-1),part.','Displayname','y_0=64/128');
    hold on
    plot((s-1:e-1),part2.','Displayname','y_0=66/128');
    plot((s-1:e-1),part3.','Displayname','y_0=68/128');
   xlabel('$\frac{x}{\Delta x}$','Interpreter','latex');
   ylabel('$T(x,y_0)$','Interpreter','latex');
    hold off
    legend('Location','northwest');
