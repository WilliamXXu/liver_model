function [res]=rebalance(c,h)
    k1=9200/1.34;
    k2=26^2.73;
    k3=c+h;
    f=@(x) x^3.73+(k1-k3)*x^2.73+k2*x-k2*k3;
    df=@(x) 3.73*x^2.73+2.73*(k1-k3)*x^1.73+k2;
    
    res=Newton(f,df,c,1);
    
    function [sol] = Newton(f,df,x0,epsilon)
        sol = x0 - f(x0)/df(x0);
        
        
        %count=0;
        while abs(f(sol)) > epsilon
            sol = sol - f(sol)/df(sol);
         %   count=count+1
        end
    end
end

