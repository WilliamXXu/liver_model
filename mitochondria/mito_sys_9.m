
function out=mito_sys_9(t,y,c,para,drug)  %para=[k_etc   v_max_uncp k_uncp    k_atpase      k_mito  n_mito  v_max_atp  ]  drug=[etc uncp atpase]
    
conversion=1.34e-3;

v_max=0.044/conversion;  %mmHg/sec  
k_max=6.24e-3/conversion; %mmHg


    drug_etc=1/(1+drug(1)/para(1));
    drug_uncp=y(3)*para(2)*drug(2)/(para(3)+drug(2));
    drug_atpase=1/(1+drug(3)/para(4));
    
    
    ocr=199.6/60/9*100/1.34;
    ox_basal=0.3;
    atp_basal=199.6*5;
    mmp_basal=0.1722;
    etc2ocr=ocr/0.6; 
    ox2ocr=(drug_etc/ox_basal)*c*v_max/(c+k_max);
    ox2etc=ox2ocr/etc2ocr;

    conversion1=198.650/0.6;
    conversion2=198.650/atp_basal;
    conversion3=1/90;
    consumption=199.6*5*conversion3/3;
    

    k_mitoox=para(5);
    n1=para(6);
    v_max_atp=para(7);%199.6*5*1.25
    k_atp=mmp_basal*(para(7)-atp_basal)/atp_basal;%0.15;


    out=[y(1)-0.6*(k_mitoox^n1+ox_basal^n1)/(k_mitoox^n1+y(2)^n1)*(1+(mmp_basal/para(8))^para(9))/(1+(y(3)/para(8))^para(9))   %pry2ox flux
    -y(2)*ox2etc+y(1) %ox  dynamics
    conversion1*y(2)*ox2etc-conversion2*y(4)-drug_uncp      %mmp dynamics
    y(4)-y(3)*v_max_atp/(k_atp+y(3))*8/(4+y(5))*drug_atpase %ATP_mitocondria flux
    4/3*conversion3*y(4)-consumption*y(5)];   %ATP_cellular dynamics
end