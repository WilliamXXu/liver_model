function holder=validation(in)   %in=[drug_type drug_conc return_option]
    %----------------------  parameters  ----------------------
error=1e-4;
timestep=1;

initial=[0.6 0.3 0.1722 199.6*5 4.0];
    M=zeros(5);
    M(2,2)=1;
    M(3,3)=1;
    M(5,5)=1;
options = odeset('Mass',M);

drug=[0 0 0];
drug(in(1))=in(2);

para=load('r_9_ms.mat').r;
%---------------------- Bio parameters  ----------------------
ocr=199.6/60/9*100/1.34;
oxygen_basal=-(312*ocr)/(67*ocr - 2200);
c=oxygen_basal;

conversion=1.34e-3;
v_max=0.044/conversion;  %mmHg/sec  
k_max=6.24e-3/conversion; %mmHg


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





last_prev=initial;
m_t=1/0;
count=0;
tic
while m_t>error

[t,y]=ode15s(@(t,y)mito_sys_9(t,y,oxygen_basal,para,drug),[0 timestep/60],last_prev,options);
last=y(size(y,1),:);
m_t=max(abs((last-last_prev)./last))/timestep;
count=count+1;
last_prev=last;
end

res=[last(2)*ox2ocr/ocr  last(1)/ecar2_py2ox/ecar_basal last(3)/initial(3) last(5)/initial(5)];
%holder=last(2)/0.6;
%mmp_basal=0.1722;
%holder=(1+(mmp_basal/0.1657)^30)/(1+(last(3)/0.1657)^30);


%  k_mitoox=para(5);
%     n1=para(6);
% 
%    holder= 0.6*(k_mitoox^n1+ox_basal^n1)/(k_mitoox^n1+y(2)^n1);


%holder=last(2)*ox2ocr/ocr;
holder=res(in(3));


end