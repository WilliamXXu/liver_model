function [res]=poisson_insulate_inf(N,size,apart)

    function [result]=normal(x,x_old,delta,i,j,apart)
        result=(1-delta)*x_old(i,j)+delta*0.25*(x(i-1,j)+x(i+1,j)+x(i,mod(j-2,apart)+1)+x(i,mod(j,apart)+1));
    end
    
   function [result]=horizontal(x,x_old,delta,i,j,direct,apart)
        result=(1-delta)*x_old(i,j)+delta*0.25*(2*x(i+direct,j)+x(i,mod(j-2,apart)+1)+x(i,mod(j,apart)+1));
    end
    function [result]=cornor(x,x_old,delta,i,j,direct_hor,direct_ver,apart)
         result=(1-delta)*x_old(i,j)+delta*0.25*(2*x(i+direct_hor,mod(j+direct_ver-1,apart)+1)+x(i+direct_hor,mod(j-direct_ver-1,apart)+1)+x(i-direct_hor,mod(j+direct_ver-1,apart)+1));
    end

function [flag,tick]=classify(j,size,apart,first)
tick=5;
    ind=j-1;
            switch ind    
                  case mod(first+size,apart)
                      flag=2;tick=-1;
                  case first
                      flag=2;tick=1;
                  otherwise
                       if  and(ind>first,ind<first+size)
                           flag=-1;
                       else
                           flag=1;
                       end
            end
end
function [res]= evaluate (x,apart,N,wall,class)
    res=zeros(N+1,apart);
    for j=1:apart
            for i=2:(wall-1)                
                res(i,j)=normal(x,x,1,i,j,apart)-x(i,j);
            end
                   flag=class(j,1);tick=class(j,2);
           switch flag
                case -1   %blank
                  res(wall,j)=normal(x,x,1,wall,j,apart)-x(wall,j);
                  res(wall+1,j)=normal(x,x,1,wall+1,j,apart)-x(wall+1,j);
                    
               case 1 % block
                    res(wall,j)=horizontal(x,x,1,wall,j,-1,apart)-x(wall,j);
                 res(wall+1,j)=horizontal(x,x,1,wall+1,j,1,apart)-x(wall+1,j);
                
               case 2  %up
                 res(wall,j)=cornor(x,x,1,wall,j,-1,tick,apart)-x(wall,j);
                 res(wall+1,j)=cornor(x,x,1,wall+1,j,1,tick,apart)-x(wall+1,j);
           end
            for i=(wall+2):N                
                res(i,j)=normal(x,x,1,i,j,apart)-x(i,j);
            end 
    end
     res=4*res;
end


    h=1/N;
    err_ord=3*log10(h);
    wall=floor((N+1)/2);
    first=ceil((apart-size)/2);
        class=zeros(apart,2);
    for j=1:apart
        [class(j,1),class(j,2)]=classify(j,size,apart,first);
    end

    for delta=1.9:-0.1:1
        qi=1;
    x=zeros(N+1,apart);
    x(N+1,:)=1;
    zuida=1;
    counter=0;
         while log10(zuida)>err_ord
        x_old=x;
        for j=1:apart
            for i=2:wall-1              
                x(i,j)=normal(x,x_old,delta,i,j,apart);
         
            end
            
          flag=class(j,1);tick=class(j,2);

           switch flag
                case -1   %blank
               x(wall,j)=normal(x,x_old,delta,wall,j,apart);
               x(wall+1,j)=normal(x,x_old,delta,wall+1,j,apart);
               case 1 % block
                    x(wall,j)=horizontal(x,x_old,delta,wall,j,-1,apart);
                 x(wall+1,j)=horizontal(x,x_old,delta,wall+1,j,1,apart);      
               case 2  %up
                 x(wall,j)=cornor(x,x_old,delta,wall,j,-1,tick,apart);
                 x(wall+1,j)=cornor(x,x_old,delta,wall+1,j,1,tick,apart);
           end
            
            for i=wall+2:N          
                x(i,j)=normal(x,x_old,delta,i,j,apart);
            end

        end 
       
        if mod(counter,100)==0
            err=evaluate (x,apart,N,wall,class);
            zuida=max(max(err)) ;       
        if zuida>2
            qi=-1;
            break
        end
        end
            counter=counter+1;
        
        end
         if qi==1
             break
         end
    end
    contour(x.');
        xlabel('$\frac{x}{\Delta x}$','Interpreter','latex');
    ylabel('$\frac{y}{\Delta y}$','Interpreter','latex');
    res=x;
end