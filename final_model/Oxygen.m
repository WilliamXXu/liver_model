classdef Oxygen < handle
    properties (Access=private)
        par
        method
        A
        B
        D
        b
        d
        option
        list_sinu
        listHepato
    end
    properties
                C
                H
                b_current
    end
    
    methods 
        function obj=Oxygen(par,option,simple,listHepato)
            if simple==1
                obj.C=[par.oxy_pp par.oxy_pp*0.7 par.oxy_pp*0.49].';
            else

            obj.par=par;
         obj.method=Method(par);
         obj.option=option;
            up=@(x,y)(obj.method.up(x,y));
            obj.C=load('result_predefined.mat').C;  
            obj.H=load('hemo_predefined.mat').H;  

                obj.listHepato=listHepato;
            list_sinu=[];
        for i=2:(par.length+1)    
             for j=2:(par.sinusoid+1)
                n=up(i,j);
                list_sinu=[list_sinu n];  
             end
         end
            obj.list_sinu=list_sinu;





timestep=par.time(1); %second, time step
unit=par.unit;
velo=par.velo;
sinusoid=par.sinusoid;
radius=par.radius;
length=par.length;
total=par.total;

miu1=par.miu1;
miu2=par.miu2;

b=zeros(total,1);
d=zeros(total,1);

f=timestep/unit;   %check courant for upwind method'
if abs(f*velo)>1
    error('courant number>1. upwind method diverges. quit')
end

%boundary conditions: oxygen concentration and flux at entry


sA=[]; %sparse matrix for iteration matrices A and B
sB=[]; 
sD=[];

%construction of A and B, D
corners=[1 up(1,sinusoid+2) up(1,radius+2) up(length+2,1) up(length+2,sinusoid+2) up(length+2,radius+2)];  %4 corners and 2 ends of the interface set to 0
for k=1:6
    n=corners(k);
    sA=[sA;n n 1];
end

for j=2:(sinusoid+1)   %sinusoid lower/upper boundary
     l=up(1,j); 
     m=up(2,j);
     n=up(3,j);

     p=up(length+2,j);  
     q=up(length+1,j);
     r=up(length,j);

     sA=[sA;l l 1;m m 1;p p 1;q q 1;];
     sB=[sB;q q 1-f*velo;q r f*velo ];
     b(m)=par.oxy_pp;
     d(m)=par.hemo_pp;
end

for j=(sinusoid+3):(radius+1)   %hepatocyte lower/upper boundary, zero flux
     n=up(1,j);   %lower
     m=up(3,j);          
     p=up(length+2,j);  %higher
     q=up(length,j);

     sA=[sA;n n 1;n m -1;p p 1;p q -1];
end


for i=2:(length+1)    %inside
    n=up(i,1);  %sinusoid left boundary, symmetry
    sA=[sA;n n 1;n up(i,3) -1];
if and(i>2,i<length+1)
    j=2;
    n=up(i,j);
    sA=[sA;n up(i+1,j) -0.5*miu1;n up(i-1,j) -0.5*miu1;n up(i,j+1) -0.5*miu1;n up(i,j-1) -0.5*miu1;n n 1+2*miu1];
    sB=[sB;n up(i+1,j) 0.5*miu1;n up(i-1,j) 0.5*miu1+f*velo;n up(i,j+1) 0.5*miu1;n up(i,j-1) 0.5*miu1;n n 1-2*miu1-f*velo];
    sD=[sD;n up(i-1,j) f*velo;n n 1-f*velo];        

    for j=3:(sinusoid+1)  %sinusoid
       
        n=up(i,j);
        sA=[sA;n up(i+1,j) -0.5*miu1;n up(i-1,j) -0.5*miu1;n up(i,j+1) -0.5*miu1-0.25*miu1/(j-2);n up(i,j-1) -0.5*miu1+0.25*miu1/(j-2);n n 1+2*miu1];
        sB=[sB;n up(i+1,j) 0.5*miu1;n up(i-1,j) 0.5*miu1+f*velo;n up(i,j+1) 0.5*miu1+0.25*miu1/(j-2);n up(i,j-1) 0.5*miu1-0.25*miu1/(j-2);n n 1-2*miu1-f*velo];
  
        sD=[sD;n up(i-1,j) f*velo;n n 1-f*velo];
      
    end
end
    j=sinusoid+2;    %interface
    n=up(i,j);
    sA=[sA; n n -miu1-miu2;n up(i,j+1) miu2;n up(i,j-1) miu1];

    
    for j=(sinusoid+3):(radius+1)  %hepatocytes
        n=up(i,j);
        sA=[sA;n up(i+1,j) -0.5*miu2;n up(i-1,j) -0.5*miu2;n up(i,j+1) -0.5*miu2-0.25*miu2/(j-2);n up(i,j-1) -0.5*miu2+0.25*miu2/(j-2);n n 1+2*miu2];
        sB=[sB;n up(i+1,j) 0.5*miu2;n up(i-1,j) 0.5*miu2;n up(i,j+1) 0.5*miu2+0.25*miu2/(j-2);n up(i,j-1) 0.5*miu2-0.25*miu2/(j-2);n n 1-2*miu2];
    end 
    n=up(i,radius+2);       %hepatocytes right boundary, zero flux
    sA=[sA;n n 1;n up(i,radius) -1];
end

%sparse matrix defined
obj.b=b;
obj.d=d;
obj.A=sparse(sA(:,1),sA(:,2),sA(:,3),total,total);
obj.B=sparse(sB(:,1),sB(:,2),sB(:,3),total,total);
obj.D=sparse(sD(:,1),sD(:,2),sD(:,3),total,total);
            end
        end










        function iterate(obj,mito,drug)

c_t=1/0;
while c_t>obj.par.error(1)
%while count<iteration
    b_cur=obj.b;  
    %count
    list_mito=obj.listHepato;
    b_cur(list_mito)=-obj.par.time(1)*respiration(obj.C(list_mito),mito(list_mito,:),drug); % update the oxygen consumption term

    obj.H=obj.D*obj.H+obj.d;  %iterate hemoglobin
    result= obj.A\(obj.B*obj.C+b_cur); %iterate disolved oxygen
    

    p=obj.method.parti(obj.list_sinu,obj.par.batch_number(1));
            for q=1:numel(p)
                activ=p{q};
                % update the hemoglobin oxygen release
                [result(activ),obj.H(activ)]=rebalance(result(activ),obj.H(activ));
            end

    C_previous=obj.C;   %update
    obj.C=result;
    c_t=norm(obj.C-C_previous)/norm(obj.C)/obj.par.time(1);
end

obj.b_current=b_cur;

        function res=respiration(C,mito,drug)
          para=obj.par;
          MM=@(x,y,z)obj.method.MM(x,y,z);
                drug_etc=MM(drug.coeff{1},1,drug.dose(1));
                if obj.option==1
                    res=MM(C,para.v_respire,para.k_respire)*drug_etc/para.basal(4).*mito(:,4).*(0.5255*mito(:,7)+0.4573);
                else
                    res=MM(C,para.v_respire,para.k_respire)*drug_etc/para.basal(4).*mito(:,4);
                end
            end

function [r_c,r_h]=rebalance(c,h)  %oxygen
    k1=9200/1.34;
    k2=26^2.73;
    k3=c+h;
    f=@(x) x.^3.73+(k1-k3).*x.^2.73+k2*x-k2.*k3;
    df=@(x) 3.73*x.^2.73+2.73*(k1-k3).*x.^1.73+k2;
 
    r_c=Newton(f,df,c,1);
    r_h=c+h-r_c;

    function [sol] = Newton(f,df,x0,epsilon)
        sol = x0 - f(x0)./df(x0);
        
        
       ct=0;
        while max(abs(f(sol)))> epsilon
            sol = sol - f(sol)./df(sol);
            ct=ct+1;
        end

    end
end



       end


    end

end