function [res]=drug_test(eff_k_m,c_max,time_max)  %0.03250 for k_m=0.35, c_max=5.1     0.06503 for k_m=0.75, c_max=1.2


para_ind=1;
drug_type=1;
tspan=1;

para=load('r_9_ms.mat').r;
para(para_ind)=eff_k_m;


drug=[0 0 0];
%---------------------- Bio parameters  ----------------------
initial=[0.6 0.3 0.1722 199.6*5 4.0];
    M=zeros(5);
    M(2,2)=1;
    M(3,3)=1;
    M(5,5)=1;
options = odeset('Mass',M);


ocr=199.6/60/9*100/1.34;
oxygen_basal=-(312*ocr)/(67*ocr - 2200);
c=oxygen_basal;

conversion=1.34e-3;
v_max=0.044/conversion;  %mmHg/sec  
k_max=6.24e-3/conversion; %mmHg



last_prev=initial;
res=[];
for t=1:tspan:time_max
t
drug(drug_type)=dr(t/60);

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

[t,y]=ode15s(@(t,y)mito_sys_9(t,y,oxygen_basal,para,drug),[0 tspan],last_prev,options);
last=y(size(y,1),:);
last_prev=last;
res=[res;[last(2)*ox2ocr/ocr last(3)/initial(3) last(5)/initial(5)]];
end






function [conc]=dr(t)   %ez drug concentration
    t_max=2;
    grad=c_max/t_max;
    half_life=4.5;
    k=log(2)/half_life;

    conc=(t<=t_max)*t*grad+(t>t_max)*c_max*exp(-k*(t-t_max));

end
end


