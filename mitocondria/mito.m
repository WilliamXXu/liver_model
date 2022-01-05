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


tspan=[0 10];
initial=[0.7 0.2 0.1722 199.6*5 4.0];
    M=zeros(5);
    M(2,2)=1;
    M(3,3)=1;
    M(5,5)=1;
options = odeset('Mass',M);



ocr=199.6/60/9*100/1.34;

oxygen_basal=-(312*ocr)/(67*ocr - 2200);
drug=[00 0 0];
%para=[0.2225    0.2987    0.2350    3.0432   0  32.2711 1576   0.0842];
para=[0.4014 0.3000 0.2000 3.1001 0 32.3000 1581.3 0.1094];

%tic;
%sol=ode15s(@(t,y)mito_sys(t,y,oxygen_basal),tspan,initial,options);
%toc;
%save('results.mat','sol')


[t,y]=ode15s(@(t,y)mito_sys(t,y,oxygen_basal,para,drug),tspan,initial,options);
siz=size(y);
last=y(siz(1),:)
plot(t,y(:,5),'-o');


function out=mito_sys(t,y,c,para,drug)  %para=[k_etc   v_max_uncp k_uncp    k_atpase      k_mito  n_mito  v_max_atp k_atp]  drug=[etc uncp atpase]
    
    global v_max
    global k_max


    drug_etc=1/(1+drug(1)/para(1));
    drug_uncp=y(3)*para(2)*drug(2)/(para(3)+drug(2));
    drug_atpase=1/(1+drug(3)/para(4));

    ocr=199.6/60/9*100/1.34;
    ox_basal=0.3;
    etc2ocr=ocr/0.6; 
    ox2ocr=(drug_etc/ox_basal)*c*v_max/(c+k_max);
    ox2etc=ox2ocr/etc2ocr;

    conversion1=198.650/0.6;
    conversion2=198.650/199.6/5;
    conversion3=1/90;
    consumption=199.6*5*conversion3/3;

    
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