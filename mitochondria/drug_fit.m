function [r]=drug_fit(k_m)
k_m
%--------------------  data prep ---------
para_ind=1;
drug_type=1;
 dec= @(b,x)(b(1)./(b(1)+x));
 beta=k_m;
 in=[0.0005 0.001 0.01 0.015 0.03 0.06 0.12 0.25 0.3 0.5 0.6 0.75 0.8 1 1.2 1.5 2 2.5 3 3.5 4 4.5 5 5.5];
 out=arrayfun((@(x)dec(beta,x)),in);

 error=1e-4;
timestep=1;

initial=[0.6 0.3 0.1722 199.6*5 4.0];
    M=zeros(5);
    M(2,2)=1;
    M(3,3)=1;
    M(5,5)=1;
options = odeset('Mass',M);
%-------------------- optimisation ------------
original=load('r_9_ms.mat').r;
r=original(para_ind);
A=-eye(numel(r)); 
b=zeros(numel(r),1);
%options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');
tic
opt = optimoptions('patternsearch','MaxFunctionEvaluations',300,'FunctionTolerance',1e-4);%,'Display','iter');
r=patternsearch(@objectiv,r,A,b,[],[],[],[],[],opt);
toc

%objectiv(r)

% beta_fit = nlinfit(in,valida(r),dec,beta)
% 
 semilogx(in,out,'-x')
 hold on
 semilogx(in,valida(r),'-x')
 hold off
 legend('data','simu');





%----------- objective function --------
    function[res]=objectiv(par)
        res=norm(valida(par)-out);
    end

function [re]=valida(par)
re=[];
 for x=1:size(in,2)

drug=[0 0 0];
drug(drug_type)=in(x);
para=load('r_9_ms.mat').r;
para(para_ind)=par;
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

re=[re last(2)*ox2ocr/ocr];
 end
end




    %----------------------  parameters  ----------------------

%holder=last(2)/0.6;
%mmp_basal=0.1722;
%holder=(1+(mmp_basal/0.1657)^30)/(1+(last(3)/0.1657)^30);


%  k_mitoox=para(5);
%     n1=para(6);
% 
%    holder= 0.6*(k_mitoox^n1+ox_basal^n1)/(k_mitoox^n1+y(2)^n1);


%holder=last(2)*ox2ocr/ocr;


% inc=@(b,x)(x*b(1)/(x+b(2)));   %MM kinetics
%  beta0=[0.6];
%  beta = nlinfit(t.Var1,atp,dec,beta0);





end
