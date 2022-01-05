function [res]=poisson_insulate_cold(N,size,apart)
     function [first,indic]=get_first(N,size,apart)
        indic=1;
    insu=apart-size;
    remain=mod(N,apart);
    if remain<=insu
       first=ceil((insu+remain)/2);
    elseif remain>size
         first=ceil((remain-size)/2);
    else
       indic=-1;
       [first,~]=get_first(N,insu,apart);
    end
    
     end

    function [flag,tick]=classify(j,size,apart,N,first,indicator)
            tick=5;      
            ind=mod(j-1,apart);
           if ind==mod(first+size,apart)    
                  flag=2; tick=-indicator;
           elseif ind==first
                flag=2;tick=indicator;
           elseif or(ind<first,ind>first+size)
               flag=1*indicator;
           else
                flag=-1*indicator;
           end
           
           if and(j==2,indicator<0)
               flag=2;tick=1;
           end
           if and(j==N,indicator<0)
               flag=2;tick=-1;
           end
           
           
           if j==1
              flag=3;tick=1;
           end
           if j==N+1
              flag=3;tick=-1;
           end 
    end
    
    
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
    
    function [result]=horizontal(x,x_old,delta,i,j,direct)
        result=(1-delta)*x_old(i,j)+delta*0.25*(2*x(i+direct,j)+x(i,j+1)+x(i,j-1));
    end
    function [result]=cornor(x,x_old,delta,i,j,direct_hor,direct_ver)
         result=(1-delta)*x_old(i,j)+delta*0.25*(2*x(i+direct_hor,j+direct_ver)+x(i+direct_hor,j-direct_ver)+x(i-direct_hor,j+direct_ver));
    end
    function [result]=cornor_new(x,x_old,delta,i,j,direct_horizontal,direct_vertical)
        result=(1-delta)*x_old(i,j)+delta*0.5*(x(i+direct_horizontal,j)+x(i,j+direct_vertical));
    end

    h=1/N;
    err_ord=3*log10(h);
    wall=floor((N+1)/2);
    [first,indicator]=get_first(N,size,apart);
    if indicator<0
       size=apart-size; 
    end
    class=zeros(N+1,2);
    for j=1:N+1
        [class(j,1),class(j,2)]=classify(j,size,apart,N,first,indicator);
    end
    
    
   
    
    for delta=1.9:-0.1:1
        qi=1;
            x=zeros(N+1,N+1);
            err=zeros(N+1,N+1);
             x(N+1,:)=1;
             err(N,:)=1;
             zuida=1;
         while log10(zuida)>err_ord
        x_old=x;
        err_old=err;
        for j=1:N+1
            for i=2:wall-1
                
                x(i,j)=normal(x,x_old,delta,i,j,N);
                err(i,j)=normal(err,err_old,delta,i,j,N);
                
            end
            
            flag=class(j,1);tick=class(j,2);
           switch flag
                case -1   %blank
               x(wall,j)=normal(x,x_old,delta,wall,j,N);
               x(wall+1,j)=normal(x,x_old,delta,wall+1,j,N);
                err(wall,j)=normal(err,err_old,delta,wall,j,N);
                err(wall+1,j)=normal(err,err_old,delta,wall+1,j,N);
               case 1 % block
                    x(wall,j)=horizontal(x,x_old,delta,wall,j,-1);
                 x(wall+1,j)=horizontal(x,x_old,delta,wall+1,j,1);
                 err(wall,j)=horizontal(err,err_old,delta,wall,j,-1);
                    err(wall+1,j)=horizontal(err,err_old,delta,wall+1,j,1);
               case 2  %up
                 x(wall,j)=cornor(x,x_old,delta,wall,j,-1,tick);
                 x(wall+1,j)=cornor(x,x_old,delta,wall+1,j,1,tick);
                 err(wall,j)=cornor(err,err_old,delta,wall,j,-1,tick);
                 err(wall+1,j)=cornor(err,err_old,delta,wall+1,j,1,tick);
              
               case 3 % cornor up
                 x(wall,j)=cornor_new(x,x_old,delta,wall,j,-1,tick);
                 x(wall+1,j)=cornor_new(x,x_old,delta,wall+1,j,1,tick);
                 err(wall,j)=cornor_new(err,err_old,delta,wall,j,-1,tick);
                 err(wall+1,j)=cornor_new(err,err_old,delta,wall+1,j,1,tick);
              

           end
            
            for i=wall+2:N
                
                x(i,j)=normal(x,x_old,delta,i,j,N);
                err(i,j)=normal(err,err_old,delta,i,j,N);
                
            end

        end 
%        % x
        delta
        zuida=max(max(abs(err)))
        
        if zuida>2
            qi=-1;
            break
        end
        
        end
         if qi==1
             break
         end
    end
%     contour(x.');
%         xlabel('$\frac{x}{\Delta x}$','Interpreter','latex');
%     ylabel('$\frac{y}{\Delta y}$','Interpreter','latex');
    res=x;
end