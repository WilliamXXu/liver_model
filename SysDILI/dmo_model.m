%with everything,mito_sys_9
function holder=dashboard(in)   % in=[drug_type drug_conc], holder=average

%----------------------  canvas -----------------------------------------



%duration=1800;
c_max=1.3;
eff=0.0338;
para_ind=1;
drug_type=1;



tmax=180;   %minute
time=[4.5e-3 60];   % seconds timestep=[oxygen dynamics, mitochondria]
pinlv=5;


error=[1e-4 1e-4 1e-5];      %[oxygen cellular ATP drug] mmHg/second  unit/second


unit=2; % miu m, space step
hepatocyte_diamter=27.5;   
sinusoid_diameter=12.5;
sinusoid_length=275;

sinusoid=round(sinusoid_diameter/2/unit);
radius=round(sinusoid+hepatocyte_diamter/unit);
length=round(sinusoid_length/unit);
total=(length+2)*(radius+2)  ;


% A c^(n+1)= B c^n+ b
% h^(n+1)= D h^n+d
C=load('result_predefined.mat').C;     %oxygen
%C=zeros(total,1);
b=zeros(total,1);
%H=zeros(total,1);
H=load('hemo_predefined.mat').H;    %haemoglobin
d=zeros(total,1);
M=zeros(total,6);  %mitochondria states
DRUG=[0 0 0];



%----------------------  common parameters ---------------------------------
conversion=1.34e-3;    %mmHg to mol/m3, Henry's law

ocr=199.6/60/9*100/1.34;

para_mito=load('r_9_ms.mat').r;
para_mito(para_ind)=eff;
para_den=[0.94281 99.9289 0.0405821 99.95526 0.39389455];
    %0.001137430834279   9.981686756318517   0.000227374467196   0.047024090440809  14.777696658267223];

v_respire=0.044/conversion;  %mmHg/sec  
k_respire=6.24e-3/conversion; %mmHg

ox_basal=0.3;
%den_basal=2.667;   % e5 cells per mm3

%----------------------  common functions -------------------------------
    function res=HILL(c,k,n)
        res=1/(1+(k/c)^n);
    end
%Michalis-Menton
function [res]=MM(c,v,k)
    res=c*v/(c+k);
end

function [res]=up(i,j)  %matrix index to vector index
    res=(i-1)*(radius+2)+j;
end


function [c]=matrix_assemble(C)
    c=zeros(radius+2,length+2);
    for n=1:total
        c(n)=C(n);
    end
end

    function res=respiration(n,C,mito,drug)
                
                drug_etc=1/(1+drug(1)/para_mito(1));
                res=MM(C(n),v_respire,k_respire)*drug_etc/ox_basal*mito(n,2)*mito(n,6);
            end

%---------------------- oxygen dynamics setup----------------------------
for w=1:1
timestep=time(1); %second, time step

%iteration=100;

%advection and diffusion   
velo=400;   %miu/sec?
f=timestep/unit;   %check courant for upwind method'
if abs(f*velo)>1
    'courant number>1. upwind method diverges. quit'
    return
end

d1=1620;     % diffusion coeff of sinusoid    miu m^2/sec
d2=1600 ;    %diffusion coeff of hepatocyte     3600(not sure)

miu1=timestep*d1/(unit*unit);
miu2=timestep*d2/(unit*unit);

%boundary conditions: oxygen concentration and flux at entry
oxy_pp=70;     %mmHg 
hemo_pp=4*2300/1.34/(1+(26/oxy_pp)^2.73); %mmHg



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
     b(m)=oxy_pp;
     d(m)=hemo_pp;
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
A=sparse(sA(:,1),sA(:,2),sA(:,3),total,total);
B=sparse(sB(:,1),sB(:,2),sB(:,3),total,total);
D=sparse(sD(:,1),sD(:,2),sD(:,3),total,total);
end





%---------------------- mito/drug/density setup  ------------------------------
initial=[0.6 0.3 0.1722 199.6*5 4.0 1];  %basal state

    for i=2:(length+1)    % initialise mitochondria
        for j=(sinusoid+3):(radius+1)
            n=up(i,j);
            M(n,:)=initial;
        end
    end




%---------------------- iterate ---------------------------------------





m_t=1/0;
count=0;
tim=0;

%while m_t>error(2)
while tim<=tmax
    tic

    tim=count*time(2)/60


    



    DRUG(drug_type)=dr(tim);
    [C,H,b_current]=oxygen(C,H,M,DRUG);
    M=mito(C,M,DRUG);


if mod(count,pinlv)==0
        oxy=matrix_assemble(C);
hemo=matrix_assemble(H);
resp=matrix_assemble(-b_current/timestep);
%remove boundry units
oxy=oxy(2:radius+1,2:length+1);
hemo=hemo(2:radius+1,2:length+1);
resp=resp(2:radius+1,2:length+1);

for ind=1:6
    res=matrix_assemble(M(:,ind));
    res=res(2:radius+1,2:length+1);
    res=res/initial(ind);
    save(append('data/mito',num2str(ind),'_',num2str(count),'.mat'),"res");
end

%save('drug.mat','sol');
save(append('data/oxy','_',num2str(count),'.mat'),"oxy");
save(append('data/hemo','_',num2str(count),'.mat'),'hemo');
save(append('data/resp','_',num2str(count),'.mat'),'resp');
%save('mito.mat','M');
%plot
%heatmap(res);
    end

%     M_prev=M;
%     m_t= mito_error(M,M_prev)
    count=count+1;
    toc
end

%count*time(2)



%{

for iter=1:time(2):duration
    iter
    tic
    [C,H,b_current]=oxygen(C,H,M,DRUG);

    M=mito(C,M,DRUG);
    toc
end
%}


%---------------------- post-modification---------------------------------------
%map back to 2-d








%---------------------- components-------------------------------------

function [conc]=dr(t)   %ez drug concentration
    t=t/60;
    t_max=2;
    grad=c_max/t_max;
    half_life=4.5;
    k=log(2)/half_life;

    conc=(t<=t_max)*t*grad+(t>t_max)*c_max*exp(-k*(t-t_max));

end


    function [C,H,b_current]=oxygen(C,H,mito,drug)

c_t=1/0;
while c_t>error(1)
%while count<iteration
    b_current=b;  
    %count
    for i=2:(length+1)    % update the oxygen consumption term
        for j=(sinusoid+3):(radius+1)
            n=up(i,j);
            b_current(n)=-timestep*respiration(n,C,mito,drug);
        end
    end
    H=D*H+d;  %iterate hemoglobin
    result= A\(B*C+b_current); %iterate disolved oxygen

    for i=2:(length+1)    % update the hemoglobin oxygen release
        for j=2:(sinusoid+1)
            n=up(i,j);
            [result(n),H(n)]=rebalance(result(n),H(n));
        end
    end
    C_previous=C;   %update
    C=result;

    c_t=norm(C-C_previous)/norm(C)/timestep;
end



function [r_c,r_h]=rebalance(c,h)  %oxygen
    k1=9200/1.34;
    k2=26^2.73;
    k3=c+h;
    f=@(x) x^3.73+(k1-k3)*x^2.73+k2*x-k2*k3;
    df=@(x) 3.73*x^2.73+2.73*(k1-k3)*x^1.73+k2;
    
    r_c=Newton(f,df,c,1);
    r_h=c+h-r_c;
    function [sol] = Newton(f,df,x0,epsilon)
        sol = x0 - f(x0)/df(x0);
        
        
        %count=0;
        while abs(f(sol)) > epsilon
            sol = sol - f(sol)/df(sol);
         %   count=count+1
        end
    end
end



end
function mitocon=mito(c,mitocon,drug)

    Mt=zeros(6);
    Mt(2,2)=1;
    Mt(3,3)=1;
    Mt(5,5)=1;
    Mt(6,6)=1;
options = odeset('Mass',Mt);

%count=0;
    for i=2:(length+1)    
        for j=(sinusoid+3):(radius+1)
            n=up(i,j);  
       
            [t,y]=ode15s(@(t,y)mito_sys(t,y,c(n,:),drug),[0 time(2)/60],mitocon(n,:),options);
             mitocon(n,:)=y(size(y,1),:);
            %count=count+1
        end
    end

function out=mito_sys(t,y,c,drug)  %para=[k_etc   v_max_uncp k_uncp    k_atpase      k_mito  n_mito  v_max_atp ]  drug=[etc uncp atpase]


    para=para_mito;
    para2=para_den;
    
    drug_etc=1/(1+drug(1)/para(1));
    drug_uncp=y(3)*para(2)*drug(2)/(para(3)+drug(2));
    drug_atpase=1/(1+drug(3)/para(4));
   
    densi=y(6);
     apop=5e-4/3;
    logistic=[0.000499253014622   1.555726168588230];
    
    atp_basal=199.6*5;
    mmp_basal=0.1722;
    etc2ocr=ocr/0.6; 
    ox2ocr=(drug_etc/ox_basal)*MM(c,v_respire,k_respire);
    ox2etc=ox2ocr/etc2ocr;

    conversion1=198.650/0.6;
    conversion2=198.650/atp_basal;
    conversion3=1/90;
    consumption=199.6*5*conversion3/3;
    

    k_mitoox=para(5);
    n1=para(6);
    v_max_atp=para(7);%199.6*5*1.25
    k_atp=mmp_basal*(para(7)-atp_basal)/atp_basal;%0.15;
    atp_dec=(1-(y(5)/4));
    atp_dec=(atp_dec>=0)*atp_dec;


    out=[y(1)-0.6*(k_mitoox^n1+ox_basal^n1)/(k_mitoox^n1+y(2)^n1)*(1+(mmp_basal/0.1657)^30)/(1+(y(3)/0.1657)^30)   %pry2ox flux
    -y(2)*ox2etc+y(1) %ox  dynamics
    conversion1*y(2)*ox2etc-conversion2*y(4)-drug_uncp      %mmp dynamics
    y(4)-y(3)*v_max_atp/(k_atp+y(3))*8/(4+y(5))*drug_atpase %ATP_mitocondria flux
    4/3*conversion3*y(4)-consumption*y(5)   %ATP_cellular dynamics
    0];
    %(1-HILL(atp_dec,para2(1),para2(2)))*logistic(1)*densi*(1-densi/logistic(2))*(densi>0.25)-apop*densi-para2(3)*HILL(atp_dec,para2(4),para2(5))*densi]; % cell density


end
end


end