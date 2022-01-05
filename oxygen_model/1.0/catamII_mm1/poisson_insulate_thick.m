function [res,res2]=poisson_insulate_thick(N)
%res is the list of Q values, res2 is the 3D matrix of all solutions
    function [result]=normal(x,x_old,delta,i,j,apart)
        result=(1-delta)*x_old(i,j)+delta*0.25*(x(i-1,j)+x(i+1,j)+x(i,mod(j-2,apart)+1)+x(i,mod(j,apart)+1));
    end
    
   function [result]=horizontal(x,x_old,delta,i,j,direct,apart)
        result=(1-delta)*x_old(i,j)+delta*0.25*(2*x(i+direct,j)+x(i,mod(j-2,apart)+1)+x(i,mod(j,apart)+1));
   end

   function [result]=vertical(x,x_old,delta,i,j,direct,apart)
        result=(1-delta)*x_old(i,j)+delta*0.25*(2*x(i,mod(j+direct-1,apart)+1)+x(i-1,j)+x(i+1,j));
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

function [res]= evaluate (x,apart,N,wall,thick,class)
    res=zeros(N+1,apart);
    for j=1:apart
            for i=2:(wall-thick-1)                
                res(i,j)=normal(x,x,1,i,j,apart)-x(i,j);
            end
                   flag=class(j,1);tick=class(j,2);
           switch flag
                case -1   %blank
                    for k=-thick:thick
                  res(wall+k,j)=normal(x,x,1,wall+k,j,apart)-x(wall+k,j);
           
                    end
               case 1 % block
                    res(wall-thick,j)=horizontal(x,x,1,wall-thick,j,-1,apart)-x(wall-thick,j);
                 res(wall+thick,j)=horizontal(x,x,1,wall+thick,j,1,apart)-x(wall+thick,j);
                
               case 2  %up
                 res(wall-thick,j)=cornor(x,x,1,wall-thick,j,-1,tick,apart)-x(wall-thick,j);
                 for k=(-thick+1):(thick-1)
                     res(wall+k,j)=vertical(x,x,1,wall+k,j,tick,apart)-x(wall+k,j);
                 end
                 res(wall+thick,j)=cornor(x,x,1,wall+thick,j,1,tick,apart)-x(wall+thick,j);
           end
            for i=(wall+thick+1):N                
                res(i,j)=normal(x,x,1,i,j,apart)-x(i,j);
            end 
    end
     res=4*res;
end
    function [res]=interpo(x,apart,N,wall,thick,step,class)
        res=zeros(N+1,apart);
        for j=1:apart
            tem=class(j);
            if tem(1)~=1
                res(wall-thick-step:wall+thick+step,j)=interp1((-thick:thick),x(wall-thick:wall+thick,j),(-thick:(thick)/(thick+step):thick));
            end
            res(1:wall-thick-step,j)=interp1((0:wall-thick-1),x(1:wall-thick,j),(0:(wall-thick-1)/(wall-thick-step-1):wall-thick-1));
             res(wall+thick+step:N+1,j)=interp1((0:wall-thick-1),x(wall+thick:N+1,j),(0:(wall-thick-1)/(wall-thick-step-1):wall-thick-1));
        end
        
    end

    h=1/N;
    err_ord=3*log10(h);
    wall=ceil((N+1)/2);
    apart=64;
    size=8;
    first=ceil((apart-size)/2);
     class=zeros(apart,2);
    for j=1:apart
        [class(j,1),class(j,2)]=classify(j,size,apart,first);
    end
    res=[];
    step=3;
    
     x_static=zeros(N+1,apart);
     x_static(N+1,:)=1;
     li=(1:step:16);
     res2=zeros(N+1,apart,numel(li));

    for thick=li
        thick
    for delta=1.9:-0.1:1
        qi=1;
    x=x_static;
    zuida=1;
    counter=0;

         while log10(zuida)>err_ord
        x_old=x;
        for j=1:apart
            for i=2:(wall-thick-1)                
                x(i,j)=normal(x,x_old,delta,i,j,apart);        
            end         
           
       flag=class(j,1);tick=class(j,2);
           switch flag
                case -1   %blank
                    for k=-thick:thick
                  x(wall+k,j)=normal(x,x_old,delta,wall+k,j,apart);
                    end
               case 1 % block
                    x(wall-thick,j)=horizontal(x,x_old,delta,wall-thick,j,-1,apart);
                 x(wall+thick,j)=horizontal(x,x_old,delta,wall+thick,j,1,apart);
               case 2  %up
                 x(wall-thick,j)=cornor(x,x_old,delta,wall-thick,j,-1,tick,apart);
                 for k=(-thick+1):(thick-1)
                     x(wall+k,j)=vertical(x,x_old,delta,wall+k,j,tick,apart); 
                 end
                 x(wall+thick,j)=cornor(x,x_old,delta,wall+thick,j,1,tick,apart);
           end
            
            for i=(wall+thick+1):N           
                x(i,j)=normal(x,x_old,delta,i,j,apart);  
            end

        end 
       
        
         if mod(counter,1000)==0
            err=evaluate(x,apart,N,wall,thick,class);
            zuida=max(max(err))       
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
    
         
    res=[res,q_measure(x)];
    res2(:,:,(thick-1+step)/step)=x;

     x_static=interpo(x,apart,N,wall,thick,step,class);
    end
end