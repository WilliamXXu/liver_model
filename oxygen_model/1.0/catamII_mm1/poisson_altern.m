function [res,counter]=poisson_altern(N,delta)  %Gauss seidel, fast
    if nargin<2
        delta=1.9;
        if nargin<1
          N=128;
        end
    end
    h=1/N;
    err_ord=3*log10(h);
    x=zeros(N+1,N+1);
    x(N+1,:)=1;
    zuida=1;
    function [result]=normal(x,x_old,delta,i,j,N)
        up=j+1;
        down=j-1;
        if j==N+1
            up=N-1;
        end
        if j==1
            down=2;
        end
        result=(1-delta)*x_old(i,j)+delta*0.25*(x(i-1,j)+x(i+1,j)+x(i,down)+x(i,up));
    end

    function [res]=evaluate(x,N)
        res=zeros(N+1,N+1);
        for j=1:N+1
            for i=2:N
                res(i,j)=normal(x,x,1,i,j,N)-x(i,j);
            end
        end
        res=4*res;
    end
    
    
    counter=0;
    while log10(zuida)>err_ord
        x_old=x;
        
        for j=1:N+1
            for i=2:N
                x(i,j)=normal(x,x_old,delta,i,j,N);

            end
        end
        if mod(counter,100)==0
            err=evaluate(x,N);
            zuida=max(max(err))
        end
        counter=counter+1;
    end
    %contour(x.');
    %xlabel('$\frac{x}{\Delta x}$','Interpreter','latex');
    %ylabel('$\frac{y}{\Delta y}$','Interpreter','latex');
    %err_ord=4*log10(h)
    %counter
    res=x;
end