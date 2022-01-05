%----------------------  parameters  ----------------------
unit=2;
hepatocyte_diamter=27.5;   
sinusoid_diameter=12.5;
sinusoid_length=275;

sinusoid=round(sinusoid_diameter/2/unit);
global radius
radius=round(sinusoid+hepatocyte_diamter/unit);
global length
length=round(sinusoid_length/unit);
global total
total=(length+2)*(radius+2)    ;
conversion=1.34e-3;
global v_max
v_max=0.044/conversion;  %mmHg/sec  
global k_max
k_max=6.24e-3/conversion; %mmHg

r=1581.3;
%options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
tic
%r=fminunc(@optimise,r);
%tic
%r=fminsearch(@optimise,r);
toc
optimise(r)
function [error]=optimise(para)
    error=0;
    temp=para;
    %para=[temp(1:6)  1576.2 temp(7)];
    para=[0.4014 0.3000 0.2000 3.1001 0 32.3000 temp 0.1094];
    data1=load('data/training.mat');
    data2=load('data/training_graphic.mat');
    drug=[data1.drug data2.drug];
    target=[data1.target data2.target];

    num_data=0;
    for x=1:size(drug,2)
        res=iterate(para,drug(:,x));
        for y=1:size(target,1)
            if ~(target(y,x)==0)               
               adder=abs(res(y)/target(y,x)-1);               
                %if adder<1
                   num_data=num_data+1;
                   error=error+adder;
                    
                    %{
                    x
                    y
                    drug(:,x)
                    target(y,x)
                    adder
                    %}
               % end
            end
        end

    end
    error=100*error/num_data;
end


function [res]=iterate(para,drug)
%----------------------  ODE parameters  ----------------------

tspan=[0 40];
initial=[0.7 0.2 0.1722 199.6*5 4.0];
    M=zeros(5);
    M(2,2)=1;
    M(3,3)=1;
    M(5,5)=1;
options = odeset('Mass',M,'RelTol',1e-4,'AbsTol',1e-6);



%---------------------- Bio parameters  ----------------------
ocr=199.6/60/9*100/1.34;
oxygen_basal=-(312*ocr)/(67*ocr - 2200);
c=oxygen_basal;

    global v_max
    global k_max


    drug_etc=1/(1+drug(1)/para(1));

    ocr=199.6/60/9*100/1.34;
    ox_basal=0.3;
    etc2ocr=ocr/0.6; 
    ox2ocr=(drug_etc/ox_basal)*c*v_max/(c+k_max);
    ox2etc=ox2ocr/etc2ocr;

    conversion1=198.650/0.6;
    conversion2=198.650/199.6/5;
    conversion3=1/90;
    consumption=199.6*5*conversion3/3;

    

    ecar_basal=66.2;%mPH/min
    ecar2_py2ox=4.28/90*16/84;

%----------------------  post-modifications  ----------------------

[t,y]=ode15s(@(t,y)mito_sys(t,y,oxygen_basal,para,drug),tspan,initial,options);

last=y(size(y,1),:);
res=[last(2)*ox2ocr/ocr  last(1)/ecar2_py2ox/ecar_basal last(3)/initial(3) last(5)/initial(5)];




%----------------------  ODE definition  ----------------------

function out=mito_sys(t,y,c,para,drug)  %para=[k_etc   v_max_uncp k_uncp    k_atpase      k_mito  n_mito  v_max_atp k_atp]  drug=[etc uncp atpase]
    

    drug_uncp=y(3)*para(2)*drug(2)/(para(3)+drug(2));
    drug_atpase=1/(1+drug(3)/para(4));

    k_mitoox=para(5);
    n1=para(6);
    v_max_atp=para(7);%199.6*5*1.25
    k_atp=para(8);%0.15;


    out=[y(1)-0.6*(k_mitoox^n1+ox_basal^n1)/(k_mitoox^n1+y(2)^n1)*(1+(0.1722/0.1657)^30)/(1+(y(3)/0.1657)^30)   %pry2ox flux
    -y(2)*ox2etc+y(1) %ox  dynamics
    conversion1*y(2)*ox2etc-conversion2*y(4)-drug_uncp      %mmp dynamics
    y(4)-y(3)*v_max_atp/(k_atp+y(3))*8/(4+y(5))*drug_atpase %ATP_mitocondria flux
    4/3*conversion3*y(4)-consumption*y(5)];   %ATP_cellular dynamics

end

end